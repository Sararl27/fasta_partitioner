import json
import pathlib
import re
import lithops


class FastaPartitioner:

    def __init__(self, storage, bucket, key, workers):
        self.storage = storage
        self.bucket = bucket

        self.__generate_fasta_index(key, workers)

    def __get_length(self, min_range, content, data, start_base, end_base):
        start_base -= min_range
        end_base -= min_range
        len_base = len(data[start_base:end_base].replace('\n', ''))
        # name_id num_chunks_has_divided offset_head offset_bases ->
        # name_id num_chunks_has_divided offset_head offset_bases len_bases
        content[-1] = f'{content[-1]} {len_base}'

    # Generate metadata from fasta file
    def __generate_chunks(self, id, key, chunk_size, obj_size, partitions):

        min_range = id * chunk_size
        max_range = int(obj_size) if id == partitions - 1 else (id + 1) * chunk_size
        data = self.storage.get_object(bucket=self.bucket, key=key,
                                       extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
        content = []
        ini_heads = list(
            re.finditer(r"\n>", data))  # If it were '>' it would also find the ones inside the head information
        heads = list(re.finditer(r">.+\n", data))

        if ini_heads or data[0] == '>':  # If the list is not empty or there is > in the first byte (if is empty it will return an empty list)
            first_sequence = True
            prev = -1
            for m in heads:
                start = min_range + m.start()
                end = min_range + m.end()
                if first_sequence:
                    first_sequence = False
                    if id >= 1 and start - 1 > min_range:  # If it is not the worker of the first part of the file and in addition it
                        # turns out that the partition begins in the middle of the base of a sequence.
                        # (start-1): avoid having a split sequence in the index that only has '\n'
                        match_text = list(re.finditer('.*\n', data[0:m.start()]))
                        length_base = '0'
                        text = ''
                        if match_text:
                            text = match_text[0].group().split(' ')[0]
                            length_0 = len(data[match_text[0].start():m.start()].replace('\n', ''))
                            offset_0 = match_text[0].start()+min_range
                            if len(match_text) > 1:
                                offset_1 = match_text[1].start()+min_range
                                length_1 = len(data[match_text[1].start():m.start()].replace('\n', ''))
                                length_base = f"{length_0}-{length_1}"
                                offset = f"{offset_0}-{offset_1}"
                            else:
                                length_base = f"{length_0}"
                                offset = f'{offset_0}'

                        # >> num_chunks_has_divided offset_head offset_bases_split length/s first_line_before_space_or_\n
                        content.append(f">> <X> <Y> {str(offset)} {length_base} ^{text}^")  # Split sequences
                if prev != start:  # When if the current sequence base is not empty
                    if prev != -1:
                        self.__get_length(min_range, content, data, prev, start)
                    # name_id num_chunks_has_divided offset_head offset_bases
                    content.append(f"{m.group().split(' ')[0].replace('>', '')} 1 {str(start)} {str(end)}")
                prev = end

            len_head = len(heads)
            if ini_heads[-1].start() + 1 > heads[-1].start():  # Check if the last head of the current one is cut. (ini_heads[-1].start() + 1): ignore '\n'
                last_seq_start = ini_heads[-1].start() + min_range + 1 # (... + 1): ignore '\n'
                if len_head !=0: # Add length of bases to last sequence
                    self.__get_length(min_range, content, data, prev, last_seq_start)
                text = data[last_seq_start - min_range::]
                # [<->|<_>]name_id_split offset_head
                content.append(
                    f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")  # if '<->' there is all id
            elif len_head !=0: # Add length of bases to last sequence
                self.__get_length(min_range, content, data, prev, max_range)

        return {'min_range': min_range,
                'max_range': max_range,
                'sequences': content}

    def __reduce_generate_chunks(self, results):
        if len(results) > 1:
            for i, dict in enumerate(results):
                dictio = dict['sequences']
                dict_prev = results[i - 1]
                seq_range_prev = dict_prev['sequences']
                if i > 0 and seq_range_prev and dictio and '>>' in dictio[
                    0]:  # If i > 0 and not empty the current and previous dictionary and the first sequence is split
                    param = dictio[0].split(' ')
                    seq_prev = seq_range_prev[-1]
                    param_seq_prev = seq_prev.split(' ')
                    if '<->' in seq_prev or '<_>' in seq_prev:
                        if '<->' in seq_range_prev[-1]:  # If the split was after a space, then there is all id
                            name_id = param_seq_prev[0].replace('<->', '')
                            length = param[4].split('-')[1]
                        else:
                            name_id = param_seq_prev[0].replace('<_>', '') + param[5].replace('^', '')
                            length = param[4].split('-')[1]
                        offset_head = param_seq_prev[1]
                        offset_base = param[3].split('-')[1]
                        seq_range_prev.pop()  # Remove previous sequence
                        split = 0
                        # Update ranges
                        dict['min_range'] = int(offset_head)
                        dict_prev['max_range'] = int(offset_head)
                    else:
                        length = param[4].split('-')[0]
                        name_id = param_seq_prev[0]
                        offset_head = param_seq_prev[2]
                        offset_base = param[3].split('-')[0]
                        split = int(param_seq_prev[1])
                        # Update number of partitions of the sequence
                        for x in range(i - split, i):  # Update previous sequences
                            results[x]['sequences'][-1] = results[x]['sequences'][-1].replace(
                                f' {split} ',
                                f' {split + 1} ')  # num_chunks_has_divided + 1 (i+1: total of current partitions of sequence)
                    dictio[0] = dictio[0].replace(f' {param[5]}', '')  # Remove 4rt param
                    dictio[0] = dictio[0].replace(f' {param[3]} ', f' {offset_base} ')  # [offset_base_0-offset_base_1|offset_base] -> offset_base
                    dictio[0] = dictio[0].replace(f' {param[4]}', f' {length}')  # [length_0-length_1|length] -> length
                    dictio[0] = dictio[0].replace(' <X> ', f' {str(split + 1)} ')  # X --> num_chunks_has_divided
                    dictio[0] = dictio[0].replace(' <Y> ', f' {offset_head} ')  # Y --> offset_head
                    dictio[0] = dictio[0].replace('>> ', f'{name_id} ')  # '>>' -> name_id
        return results

    def __generate_json(self, data, dt_dir):
        with open(f'{dt_dir}.json', 'w') as f:
            json.dump(data, f, indent=2)

    def __generate_fasta_index(self, key, workers):
        fexec = lithops.FunctionExecutor(max_workers=1000, runtime_memory=4096) # log_level='DEBUG

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

        self.__generate_json(results, f'./output_data/{pathlib.Path(key).stem}_index')

        print('... Done, generated index')

        # fexec.plot()
        fexec.clean()


def get_info_sequence(data, identifier):
    length = offset_head = offset = -1
    found = False
    if identifier != '':
        for i, dict in enumerate(data):
            for j, sequence in enumerate(dict['sequences']):
                if identifier in sequence:
                    found = True
                    param_seq = sequence.split(' ')
                    split = int(param_seq[1])
                    if split > 1:
                        length = 0
                        for x in range(i + 1, i + split):
                            length += int(data[x]['sequences'][0].split(' ')[4])
                        length += int(param_seq[4])
                    else:
                        length = int(param_seq[4])
                    offset_head = int(param_seq[2])
                    offset = int(param_seq[3])
                    break
            if found:
                break

    return {'length': length, 'offset_head': offset_head , 'offset': offset}

def get_sequences_of_range(data, min_range, max_range):
    sequences = []
    for i, dict in enumerate(data):
        if dict['min_range'] == min_range and dict['max_range'] == max_range:
            split = int(dict['sequences'][0].split(' ')[1])
            if split > 1:
                seq_prev = data[i - split + 1]['sequences'][-1]
                sequences.append({'identifier': seq_prev.split(' ')[0], 'offset': seq_prev.split(' ')[3]})
            else :
                sequences.append(dict['sequences'][0].split(' ')[0])
            for e in dict['sequences'][1::]:
                sequences.append(e.split(' ')[0])
            break
    return sequences



