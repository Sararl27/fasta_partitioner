# from pyfaidx import Fasta
# genes = Fasta('input_data/genes.fasta')
import json
import os
import pathlib
import re

import lithops


class PartitionFasta:

    def __init__(self, storage, bucket, obj, id):
        self.storage = storage
        self.bucket = bucket
        self.obj = obj
        self.id = id

    # Generate metadata from fasta file
    def generate_chunks(self, key, chunk_size, obj_size, partitions):
        min_range = self.id * chunk_size
        max_range = (self.id + 1) * chunk_size
        data = self.storage.get_object(bucket=self.bucket, key=key,
                                       extra_get_args={'Range': f'bytes={min_range}-{max_range - 1}'}).decode('utf-8')
        content = []
        ini_heads = list(re.finditer(r"\n>", data))  # If it were '>' it would also find the ones inside the head information
        heads = list(re.finditer(r">.+\n", data))

        if ini_heads or data[0] == '>':  # If the list is not empty or there is > in the first byte (if is empty it will return an empty list)
            first_sequence = True
            prev = ""
            for m in heads:
                start = min_range + m.start()
                end = min_range + m.end()
                if first_sequence:
                    first_sequence = False
                    if self.id >= 1 and start - 1 > min_range:  # If it is not the worker of the first part of the file and in addition it
                                                                # turns out that the partition begins in the middle of the base of a sequence.
                                                                # (start-1): avoid having a split sequence in the index that only has '\n'
                        first_line = re.search('[^\s\n)]+', data).group()  # Match until space or '\n'
                        # >> num_chunks_has_divided offset_head offset_bases_split
                        content.append(f">> X Y {str(min_range)} {first_line}")  # Split sequences

                if str(prev) != str(start):  # When if the current sequence base is not empty
                                             # name_id num_chunks_has_divided offset_head offset_bases
                    content.append(f"{m.group().split(' ')[0].replace('>', '')} 1 {str(start)} {str(end)}")
                prev = str(end)

            len_ini_head = len(ini_heads)
            if len_ini_head > len(heads):  # Check if the last head of the current one is cut
                last_seq_start = ini_heads[len_ini_head - 1].start() + 1
                text = data[last_seq_start::]
                # [<->|<_>]name_id_split offset_head
                content.append(
                    f"{'<-' if ' ' in text else '<_'}{text.split(' ')[0]} {str(last_seq_start)}")  # if '<->' there is all id

        return {'min_range': str(min_range), 'max_range': str(obj_size) if self.id == partitions - 1 else str(max_range),
                'sequences': content}


def run_worker_metadata(key, id, storage, bucket, chunk_size, obj_size, partitions):
    partitioner = PartitionFasta(storage, bucket, key, id)
    return partitioner.generate_chunks(key, chunk_size, obj_size, partitions)


def reduce_generate_chunks(results):
    if len(results) > 1:
        i = 0
        for dict in results:
            dictio = dict['sequences']
            dict_prev = results[i - 1]
            seq_range_prev = dict_prev['sequences']
            if i > 0 or dict['sequences'] and seq_range_prev and '>>' in dictio[0]:  # If i > 0 and not empty the current and previous dictionary and the first sequence is split
                text = dictio[0].split(' ')[4]
                i_seq_prev = len(seq_range_prev) - 1  # Index of the last sequence of the previous range
                seq_prev = seq_range_prev[i_seq_prev]
                if '<->' in seq_prev or '<_>' in seq_prev:
                    if '<->' in seq_range_prev[i_seq_prev]:  # If the split was after a space, then there is all id
                        name_id = seq_prev.split(' ')[0].replace('<->', '')
                    else:
                        name_id = seq_prev.split(' ')[0].replace('<_>', '') + text
                    offset_head = seq_prev.split(' ')[1]
                    seq_range_prev.pop()  # Remove previous sequence
                    split = 0
                    # Update ranges
                    dict['min_range'] = offset_head
                    dict_prev['max_range'] = offset_head
                else:
                    name_id = seq_prev.split(' ')[0]
                    offset_head = seq_prev.split(' ')[2]
                    split = int(seq_prev.split(' ')[1])
                    # Update number of partitions of the sequence
                    for x in range(i - split, i):  # update previous sequences
                        results[x]['sequences'][i_seq_prev] = results[x]['sequences'][i_seq_prev].replace(
                            f' {split} ',
                            f' {split + 1} ')  # num_chunks_has_divided + 1 (i+1: total of current partitions of sequence)
                dictio[0] = dictio[0].replace(' X ', f' {str(split + 1)} ')  # X --> num_chunks_has_divided
                dictio[0] = dictio[0].replace(' Y ', f' {offset_head} ')  # Y --> offset_head
                dictio[0] = dictio[0].replace(f' {text}', '')  # remove 4rt param
                dictio[0] = dictio[0].replace('>> ', f'{name_id} ')  # '>>' -> name_id
            i += 1
    return results  # results


def push_object_funct(local_input_path, data_bucket, input_data_prefix):
    bucket_objects = storage.list_keys(bucket=data_bucket)
    for subdir, dirs, files in os.walk(local_input_path):
        print(subdir)
        for file_name in files:
            key = os.path.join(input_data_prefix, file_name)  # Added
            if key not in bucket_objects:  # Changed: if file_name not in bucket_objects:
                with open(os.path.join(subdir, file_name), 'rb') as file:  # Changed
                    print(f'\tUploading {key}...')
                    data = file.read()
                    storage.put_object(bucket=data_bucket, key=key, body=data)
                    print('\tOk!')
            else:  # Added
                print(f'\tIt is already uploaded: {key}...')  # Added
    print("Finished!")

    return storage.list_keys(bucket=data_bucket)


def generate_json(data, dt_dir):
    with open(f'{dt_dir}.json', 'w') as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    fexec = lithops.FunctionExecutor(log_level='DEBUG', max_workers=1000, runtime_memory=4096)
    storage = lithops.Storage()

    # Params
    data_bucket = 'cloud-fasta-partitioner-1'  # Change-me
    prefix = 'fasta/'  # Change-me
    local_input_path = './input_data/'  # Change-me

    key = f'fasta/genes.fasta'  # Change-me
    obj = f'{data_bucket}/{key}'  # f'cos://{data_bucket}/{key}'  # Change-me

    workers = 100  # Change-me

    # Execution
    path_obj = push_object_funct(local_input_path, data_bucket, prefix)

    # location = obj.split('//')[1].split('/')  # obj.split('/') #
    # for_head = location[1:]
    # for_head = '/'.join(for_head)
    # data_bucket_name = location[0]
    fasta = storage.head_object(data_bucket, key)
    chunk_size = int(int(fasta['content-length']) / workers)
    # seq_name = location[-1].split('.')[0]

    print('===================================================================================')
    print('object: ' + obj)
    print('metadata chunks: ' + str(chunk_size))
    print('bucket to access data: ' + str(data_bucket))
    print('reference file name: ' + pathlib.Path(key).stem)
    print('fasta size: ' + str(fasta['content-length']) + ' bytes')
    print('===================================================================================')

    map_iterdata = [{'key': key} for i in range(workers)]
    extra_args = {'bucket': data_bucket, 'chunk_size': chunk_size, 'obj_size': fasta['content-length'],
                  'partitions': workers}
    fexec.map_reduce(map_function=run_worker_metadata, map_iterdata=map_iterdata, extra_args=extra_args,
                     reduce_function=reduce_generate_chunks)
    fexec.wait()
    results = fexec.get_result()

    generate_json(results, f'./output_data/{pathlib.Path(key).stem}_index')

    i = 0
    for el in results:
        print(f"{str(i)}: {[el['min_range'], el['max_range']]}")
        i += 1

    # fexec.plot()
    fexec.clean()
