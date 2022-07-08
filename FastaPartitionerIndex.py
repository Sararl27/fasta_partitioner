# from pyfaidx import Fasta
# genes = Fasta('input_data/genes.fasta')
import pathlib
import re
import os
import lithops
import json
import json


class PartitionFasta:

    def __init__(self, storage, bucket, obj, id):
        self.storage = storage
        self.bucket = bucket
        self.obj = obj
        self.id = id


    # Generate metadata from fasta file
    def generate_chunks(self, key, id, chunk_size, obj_size, partitions):
        min_range = id * chunk_size
        max_range = (id + 1) * chunk_size
        data = self.storage.get_object(bucket=self.bucket, key=key,extra_get_args={'Range': f'bytes={min_range}-{max_range}'}).decode('utf-8')
        FASTA_HEADER_START = r"\n>" # If it were '>' it would also find the ones inside the head information
        FASTA_HEADER = r">.+\n"
        content = []
        ini_heads = list(re.finditer(FASTA_HEADER_START, data))
        heads = list(re.finditer(FASTA_HEADER, data))

        if ini_heads or data[0] == '>':  # If the list is not empty or there is > in the first byte
            first_sequence = True
            prev = ""
            for m in heads:
                start = min_range + m.start()
                end = min_range + m.end()
                if first_sequence:
                    first_sequence = False
                    if id >= 1 and start - 1 > min_range:  # If it is not the worker of the first part of the file and in addition it
                                                           # turns out that the partition begins in the middle of the base of a sequence.
                                                           # (start-1): avoid having a split sequence in the index that only has '\n'
                        # >> num_chunks_has_divided offset_head offset_bases_split
                        content.append(f">> X Y {str(min_range)}")  # Split sequences
                # ------------------------------------------------
                # TODO saber que fa i si es útil consevar-ho
                if (id == int(partitions) - 1) and (int(obj_size) - end < int(partitions) - 1):  # Si es tracta del worker que te l'última part del fitxer
                    # original i ?????
                    content.append('This should not be')
                    #content = content + str(start) + ',' + str(obj_size) + 'x\n'  # Head of the current
                # ------------------------------------------------
                else:  # When it is not the worker who has the last part of the original file
                    if str(prev) != str(start):  # When if the current sequence base is not empty
                        # name_id num_chunks_has_divided offset_head offset_bases
                        content.append(f"{m.group().split(' ')[0].replace('>', '')} 1 {str(start)} {str(end)}")

                    prev = str(end)
            
            len_ini_head = len(ini_heads)
            if len_ini_head > len(heads):  # Check if the last head of the current one is cut
                last_seq_start = ini_heads[len_ini_head - 1].start() + 1
                text = data[last_seq_start::]
                if ' ' in text:  # The split is after a space (there all id)
                    # <->name_id_split num_chunks_has_divided offset_head offset_bases
                    content.append(f"<-{text.split(' ')[0]} 1 {str(last_seq_start)} X")
                else:  # The split not is after a space (the id is split)
                    # <_>name_id_split num_chunks_has_divided offset_head offset_bases
                    content.append(f"<_{text} 1 {str(last_seq_start)} X")

        else: # If the list is empty and there is no > in the first byte
            content.append('<Sequencies not found>')

        dict = {'min_range': str(min_range), 'max_range': str(obj_size) if id == partitions - 1 else str(max_range),
                'sequences': content}

        return [dict, [[str(len(ini_heads)), str(len(heads))], [str(ini_heads[len(ini_heads)-1].start() + 1), str(heads[len(heads)-1].start())]]]


def run_worker_metadata(key, id, storage, bucket, chunk_size, obj_size, partitions):
    partitioner = PartitionFasta(storage, bucket, key, id)
    dades = partitioner.generate_chunks(key, id, chunk_size, obj_size, partitions)
    return dades


def update_split_sequence(results, dict, i, i_seq_prev, seq_prev):
    split = int(seq_prev.split(' ')[1])
    # TODO i - split: Comprovar que fa el tinc pensat
    for x in range(i - split, split + 1):
        results[x]['sequences'][i_seq_prev] = results[x]['sequences'][i_seq_prev].replace(split, split + 1)  # num_chunks_has_divided + 1 (i+1: total of current partitions of sequence)
    dict['sequences'][0] = dict['sequences'][0].replace('X', split + 1)

def reduce_generate_chunks(results):
    # · Trobar si hi ha algú worker que ha començat amb una secuencia partida i dir-li a quina secuencia pertany (>>
    #   num_particions offset -> nam_id  num_particions offset, i actualitzar num_particions de la secuencia que si
    #   tenia la part del head (0 -> x))
    # · S'ha d'actualitzar el rang en cas que just en un worker s'ha quedat el head (o parcialment) i un altre worker té tota la secuencia (i potser part del head)
    # · Comprovar només mirant la primera secuencia, si esta dividida, comprovar l'ultima secuencia de l'anterior chunk:
    #     - si té '>>>name_id': passar el head al següent chuhnk (i eliminar els simbols '>>>') i actualitzar rangs
    #     - si té '>>: copiar el head a l'index i donar-li el offset del head
    '''if len(results) > 1:
        i = 0
        for dict in results:
            if i > 0 or not ('<Sequencies not found>' in dict['sequences'] and '<Sequencies not found>' in results[i-1]['sequences']):
                i_seq_prev = len(results[i-1]['sequences']) - 1     # Index of the last sequence of the previous range
                seq_prev = results[i-1]['sequences'][i_seq_prev]
                dict_prev = results[i-1]
                if '>>' in dict['sequences'][0]:  # If the first sequence is split
                    dict['sequences'][0] = dict['sequences'][0].replace('Y', seq_prev.split(' ')[2])   # Y --> offset_head
                    update_split_sequence(results,dict, i, i_seq_prev, seq_prev)
                    dict['sequences'][0] = dict['sequences'][0].replace('>>', seq_prev.split(' ')[0])  # '>>' -> name_id
                    # TODO ->
                    # Actualizar rangos

                elif '<->' in dict['sequences'][0]: # Comprovar que el primer index no porta >>:
                elif '<_>' in dict['sequences'][0]:
                    print('')
            i += 1'''
    return results # results


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
        '''f.write(
            '[' +
            ',\n'.join(json.dumps(i) for i in data) +
            ']\n')'''


# s3://sra-pub-src-2/SRR8774337/cleaned.fasta

if __name__ == "__main__":
    # Params
    data_bucket = 'cloud-fasta-partitioner-1'  # Change-me
    prefix = 'fasta/'  # Change-me
    local_input_path = './input_data/'  # Change-me

    key = f'fasta/genes.fasta'  # Change-me
    obj = f'cos://{data_bucket}/{key}'  # f'{my_bucket_name}/{path_obj[0]}'  # Change-me
    print(obj)

    workers = 100  # Change-me

    # Execution
    fexec = lithops.FunctionExecutor(log_level='DEBUG', max_workers=1000, runtime_memory=4096)
    storage = lithops.Storage()

    path_obj = push_object_funct(local_input_path, data_bucket, prefix)

    #location = obj.split('//')[1].split('/')  # obj.split('/') #
    #for_head = location[1:]
    #for_head = '/'.join(for_head)
    #data_bucket_name = location[0]
    fasta = storage.head_object(data_bucket, key)
    chunk_size = int(int(fasta['content-length']) / workers)
    # seq_name = location[-1].split('.')[0]

    print('===================================================================================')
    print('metadata chunks: ' + str(chunk_size))
    print('bucket to access data: ' + str(data_bucket))
    print('reference file name: ' + pathlib.Path(key).stem)
    print('fasta size: ' + str(fasta['content-length']) + ' bytes')
    print('===================================================================================')

    #map_iterdata = {'bucket': data_bucket, 'obj': obj, 'chunk_size': chunk_size,'obj_size': fasta['content-length']}
    #fexec.map_reduce(map_function=run_worker_metadata, map_iterdata=map_iterdata,reduce_function=reduce_generate_chunks)

    map_iterdata = [{'key': key} for i in range(workers)]
    extra_args = {'bucket': data_bucket,'chunk_size': chunk_size, 'obj_size': fasta['content-length'], 'partitions': workers}
    fexec.map_reduce(map_function=run_worker_metadata,map_iterdata=map_iterdata, extra_args=extra_args, reduce_function=reduce_generate_chunks)
    fexec.wait()
    results = fexec.get_result()

    res = [data[0] for data in results]
    generate_json(res, f'./output_data/{pathlib.Path(key).stem}_index')
    i = 0
    for el in results:
        # print(f"{str(i)}: {[el['min_range'], el['max_range']]}")
        '''with open(f'./output_data/a/a_{i}.txt', 'w') as f:
            f.write(el)'''
        print(f"{str(i)}: {[el[0]['min_range'], el[0]['max_range']]} | {str(el[1])}")
        i += 1

    fexec.plot()
    fexec.clean()
