import os
import pathlib
import pickle
import re
import lithops
from cuckooHash import HashTab


class FastaIndex:

    def __init__(self, storage, bucket, key, workers):
        self.storage = storage
        self.bucket = bucket

        self.__generate_fasta_index(key, workers)

    def __get_length(self, min_range, hastab, key, data, start_base, end_base):
        if key != '':
            start_base -= min_range
            end_base -= min_range
            len_base = len(data[start_base:end_base].replace('\n', ''))
            seq = hastab.find(key)
            if seq is not None:
                # name_id num_chunks_has_divided offset_head offset_bases ->
                # name_id num_chunks_has_divided offset_head offset_bases len_bases
                hastab.delete(key)
                hastab.insert(key, f'{seq} {len_base}')
                return False
            return True

    def __update_number_of_partitions(self, split, i, dicts):
        for x in range(i - split, i):  # Update previous sequences
            seqs_prev = dicts[x]['sequences']
            prev = dicts[x]['name_sequence'][-1]
            seq = seqs_prev.find(prev)
            if seq is not None:
                # name_id num_chunks_has_divided offset_head offset_bases ->
                # name_id num_chunks_has_divided offset_head offset_bases len_bases
                self.__update_data_hastab_with_raplace(seqs_prev, prev, seq, f' {split} ',
                                                       f' {split + 1} ')  # num_chunks_has_divided + 1 (i+1: total of current partitions of sequence))

    def __update_data_hastab_with_raplace(self, hastab, key, sequence, before, after):
        hastab.delete(key)
        seq = sequence.replace(before, after)
        hastab.insert(key, seq)
        return seq

        # Generate metadata from fasta file

    def __generate_chunks(self, id, key, chunk_size, obj_size, partitions):
        min_range = id * chunk_size
        max_range = int(obj_size) if id == partitions - 1 else (id + 1) * chunk_size
        data = self.storage.get_object(bucket=self.bucket, key=key,
                                       extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
        ini_heads = list(
            re.finditer(r"\n>", data))  # If it were '>' it would also find the ones inside the head information
        heads = list(re.finditer(r">.+\n", data))
        len_head = len(heads)
        ht = HashTab(len_head//2)
        content = []
        miss_prev_seq = False
        if ini_heads or data[
            0] == '>':  # If the list is not empty or there is > in the first byte (if is empty it will return an empty list)
            first_sequence = True
            seq_prev = ''
            prev = -1
            for m in heads:
                start = min_range + m.start()
                end = min_range + m.end()
                if first_sequence:
                    first_sequence = False
                    if id > 0 and start - 1 > min_range:  # If it is not the worker of the first part of the file and in addition it
                        # turns out that the partition begins in the middle of the base of a sequence.
                        # (start-1): avoid having a split sequence in the index that only has '\n'
                        match_text = list(re.finditer('.*\n', data[0:m.start()]))
                        if match_text:
                            text = match_text[0].group().split(' ')[0]
                            length_0 = len(data[match_text[0].start():m.start()].replace('\n', ''))
                            offset_0 = match_text[0].start() + min_range
                            if len(match_text) > 1:
                                offset_1 = match_text[1].start() + min_range
                                length_1 = len(data[match_text[1].start():m.start()].replace('\n', ''))
                                length_base = f"{length_0}-{length_1}"
                                offset = f"{offset_0}-{offset_1}"
                            else:
                                length_base = f"{length_0}"
                                offset = f'{offset_0}'
                            # >> num_chunks_has_divided offset_head offset_bases_split length/s first_line_before_space_or_\n
                            ht.insert('first_split_sequence', f">> <X> <Y> {str(offset)} {length_base} ^{text}^")  # Split sequences
                            content.append('first_split_sequence')
                        else:  # When the first header found is false, when in a split stream there is a split header that has a '>' inside (ex: >tr|...o-alpha-(1->5)-L-e...\n)
                            first_sequence = True
                            start = end = -1  # Avoid entering the following condition
                if prev != start:  # When if the current sequence base is not empty
                    if prev != -1:
                        miss_prev_seq = self.__get_length(min_range, ht, seq_prev, data, prev, start)
                    id_seq = m.group().split(' ')[0].replace('>', '')
                    # name_id num_chunks_has_divided offset_head offset_bases
                    ht.insert(id_seq, f"{id_seq} 1 {str(start)} {str(end)}")
                    content.append(id_seq)
                    seq_prev = id_seq
                if miss_prev_seq:
                    break
                prev = end

            if not miss_prev_seq and ini_heads[-1].start() + 1 > heads[
                -1].start():  # Check if the last head of the current one is cut. (ini_heads[-1].start() + 1): ignore '\n'
                last_seq_start = ini_heads[-1].start() + min_range + 1  # (... + 1): ignore '\n'
                if len_head != 0:  # Add length of bases to last sequence
                    miss_prev_seq = self.__get_length(min_range, ht, seq_prev, data, prev, last_seq_start)
                if not miss_prev_seq:
                    text = data[last_seq_start - min_range::]
                    # [<->|<_>]name_id_split offset_head
                    ht.insert('all_id_sequence' if ' ' in text else 'half_id_sequence',
                              f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")  # if '<->' there is all id
            elif not miss_prev_seq and len_head != 0:  # Add length of bases to last sequence
                self.__get_length(min_range, ht, seq_prev, data, prev, max_range)

        return {'min_range': min_range,
                'max_range': max_range,
                'sequences': ht,  # if not miss_prev_seq else None
                'name_sequence': content}  # if not miss_prev_seq else None

    def __reduce_generate_chunks(self, results):
        '''if len(results) > 1:
            for i, dict in enumerate(results):
                if dict['sequences'] is not None:
                    dictio = dict['sequences']
                    sequence = dictio.find('first_split_sequence')
                    if sequence is not None:
                        dict_prev = results[i - 1]
                        param = sequence.split(' ')
                        seq_all_id = dictio.find('all_id_sequence')
                        seq_half_id = dictio.find('half_id_sequence')
                        if seq_all_id is not None or seq_half_id is not None:
                            if seq_all_id is not None:  # If the split was after a space, then there is all id
                                param_seq_prev = seq_all_id.split(' ')
                                name_id = param_seq_prev[0].replace('<->', '')
                            else:
                                param_seq_prev = seq_half_id.split(' ')
                                name_id = param_seq_prev[0].replace('<_>', '') + param[5].replace('^', '')
                            length = param[4].split('-')[1]
                            offset_head = param_seq_prev[1]
                            offset_base = param[3].split('-')[1]
                            dict['name_sequence'][0] = name_id  # Update data
                            split = 0
                            # Update ranges
                            dict['min_range'] = int(offset_head)
                            dict_prev['max_range'] = int(offset_head)
                        else:
                            seq_prev = dictio.find(dict_prev['name_sequence'][-1])
                            if seq_prev is not None:
                                param_seq_prev = seq_prev.split(' ')
                                length = param[4].split('-')[0]
                                name_id = param_seq_prev[0]
                                offset_head = param_seq_prev[2]
                                offset_base = param[3].split('-')[0]
                                split = int(param_seq_prev[1])
                                # Update number of partitions of the sequence
                                self.__update_number_of_partitions(split, i, results)
                        sequence = self.__update_data_hastab_with_raplace(dictio, name_id, sequence, f' {param[5]}', '')  # Remove 4rt param
                        sequence = self.__update_data_hastab_with_raplace(dictio, name_id, sequence, f' {param[3]} ', f' {offset_base} ')  # [offset_base_0-offset_base_1|offset_base] -> offset_base
                        sequence = self.__update_data_hastab_with_raplace(dictio, name_id, sequence, f' {param[4]}', f' {length}')  # [length_0-length_1|length] -> length
                        sequence = self.__update_data_hastab_with_raplace(dictio, name_id, sequence, ' <X> ', f' {str(split + 1)} ')  # X --> num_chunks_has_divided
                        sequence = self.__update_data_hastab_with_raplace(dictio, name_id, sequence, ' <Y> ', f' {offset_head} ')  # Y --> offset_head
                        self.__update_data_hastab_with_raplace(dictio, name_id, sequence, '>> ', f'{name_id} ')  # '>>' -> name_id'''
        return results

    def __generate_index_file(self, data, file_name):
            output_path = './output_data/'
            if not os.path.exists(output_path):
                os.mkdir(output_path)
            with open(f'{output_path}{file_name}_index.pkl', 'wb') as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

    def __generate_fasta_index(self, key, workers):
        fexec = lithops.FunctionExecutor(max_workers=2000, runtime_memory=4096)  # log_level='DEBUG

        # location = obj.split('//')[1].split('/')  # obj.split('/') #
        # for_head = location[1:]
        # for_head = '/'.join(for_head)
        # data_bucket_name = location[0]
        fasta = self.storage.head_object(self.bucket, key)
        chunk_size = int(int(fasta['content-length']) / workers)
        # seq_name = location[-1].split('.')[0]

        print('===================================================================================')
        print('metadata chunks: ' + str(chunk_size))
        print('bucket to access data: ' + str(self.bucket))
        print('reference file name: ' + pathlib.Path(key).stem)
        print('fasta size: ' + str(fasta['content-length']) + ' bytes')
        print('===================================================================================')

        map_iterdata = [{'key': key} for _ in range(workers)]
        extra_args = {'chunk_size': chunk_size, 'obj_size': fasta['content-length'], 'partitions': workers}
        fexec.map_reduce(map_function=self.__generate_chunks, map_iterdata=map_iterdata, extra_args=extra_args,
                         reduce_function=self.__reduce_generate_chunks)
        fexec.wait()
        results = fexec.get_result()

        self.__generate_index_file(results, f'./output_data/{pathlib.Path(key).stem}')

        print('... Done, generated index')

        # fexec.plot()
        fexec.clean()


class FunctionsFastaIndex:
    def __init__(self, path_index_file):
        with open(path_index_file, "rB") as f:
            self.data = pickle.load(f)

    def get_info_sequence(self, identifier):
        length = offset_head = offset = -1
        found = False
        if identifier != '':
            for i, dict in enumerate(self.data):
                for j, sequence in enumerate(dict['sequences']):
                    if identifier in sequence:
                        found = True
                        param_seq = sequence.split(' ')
                        split = int(param_seq[1])
                        if split > 1:
                            length = 0
                            for x in range(i + 1, i + split):
                                length += int(self.data[x]['sequences'][0].split(' ')[4])
                            length += int(param_seq[4])
                        else:
                            length = int(param_seq[4])
                        offset_head = int(param_seq[2])
                        offset = int(param_seq[3])
                        break
                if found:
                    break

        return {'length': length, 'offset_head': offset_head, 'offset': offset}

    def __get_sequences(self, dict_seqs, previous_dict_seqs, split, sequences, index):
        if split > 1:
            seq_prev = previous_dict_seqs[-1]
            sequences.append({'identifier': seq_prev.split(' ')[0], 'offset': seq_prev.split(' ')[3]})
        else:
            sequences.append(dict_seqs[index].split(' ')[0])

    def get_sequences_of_range(self, min_range, max_range):
        sequences = []
        if min_range < max_range:
            i_min_range = -1
            i_max_range = -1
            for i, dict in enumerate(self.data):
                if dict['min_range'] == min_range and dict['max_range'] == max_range:  # If it is the default ranges
                    split = int(dict['sequences'][0].split(' ')[1])
                    self.__get_sequences(dict['sequences'], self.data[i - split + 1]['sequences'], split, sequences, 0)
                    for e in dict['sequences'][1::]:
                        sequences.append(e.split(' ')[0])
                    break
                else:
                    if dict['min_range'] <= min_range <= dict['max_range']:
                        i_min_range = i
                    if dict['min_range'] <= max_range <= dict['max_range']:
                        i_max_range = i
                    if i_min_range != -1 and i_max_range != -1:
                        break

            if i_min_range != -1 and i_max_range != -1: # If it is not the default ranges
                first_pos = self.__binary_search_modified(self.data[i_min_range], min_range, 'min')
                last_pos = self.__binary_search_modified(self.data[i_max_range], max_range, 'max')
                dict_min = self.data[i_min_range]
                if i_min_range == i_max_range:    # If it is the same dictionary
                    split = int(dict_min['sequences'][first_pos].split(' ')[1])
                    self.__get_sequences(dict_min['sequences'], self.data[i - split + 1]['sequences'], split, sequences, first_pos)
                    for e in self.data[i_min_range]['sequences'][first_pos + 1:last_pos]:
                        sequences.append(e.split(' ')[0])
                else:
                    for i in range(i_min_range, i_max_range + 1):
                        if i == i_min_range:
                            start = first_pos
                            end = len(self.data[i]['sequences'])
                        elif i == i_max_range:
                            start = 0
                            end = last_pos
                        else:
                            start = 0
                            end = len(self.data[i]['sequences'])
                        split = int(self.data[i]['sequences'][start].split(' ')[1])
                        if i == i_min_range and start == 0 or split == 1:
                            self.__get_sequences(self.data[i]['sequences'], self.data[i - split + 1]['sequences'], split, sequences, start)
                        elif i == i_min_range and start == end - 1:  # If it is the last sequence in the first dictionary in the range and is split
                            sequences.append(self.data[i]['sequences'][start].split(' ')[0])
                        for e in self.data[i]['sequences'][start + 1:end]:
                            sequences.append(e.split(' ')[0])
        return sequences

    def __binary_search_modified(self, dict, x, side):
        if side == 'min' or side == 'max':
            arr = dict['sequences']
            length = len(arr) - 1
            low = 0
            high = length
            mid = 0

            while low <= high:
                mid = (high + low) // 2
                # If x is greater, ignore left half
                param = arr[mid].split(' ')
                offset = int(param[3] if param[1] != '1' and mid == 0 else param[2])
                if offset < x:
                    low = mid + 1
                # If x is smaller, ignore right half
                elif offset > x:
                    high = mid - 1
                # means x is present at mid
                else:
                    return mid

            if mid == length and offset < x:
                return dict['max_range'] if side == 'min' else length
            elif mid == 0 and x < offset:
                return 0 if side == 'min' else dict['min_range']
            else:   # 0 <= mid <= length
                if mid != 0:
                    param = arr[mid - 1].split(' ')
                    offset2 = int(param[3] if param[1] != '1' and mid == 0 else param[2])
                    if offset2 < x < offset:
                        return mid if side == 'min' else mid - 1
                if mid != length:
                    param = arr[mid + 1].split(' ')
                    offset2 = int(param[3] if param[1] != '1' and mid == 0 else param[2])
                    if offset < x < offset2:
                        return mid + 1 if side == 'min' else (mid if mid != 0 else dict['min_range'])
        return -1
