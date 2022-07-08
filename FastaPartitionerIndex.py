# from pyfaidx import Fasta
# genes = Fasta('input_data/genes.fasta')
import pathlib
import re
import os
import lithops
import json
import json
import gc


class PartitionFasta:

    def __init__(self, storage, bucket, obj, id):
        self.storage = storage
        self.bucket = bucket
        self.obj = obj
        self.id = id


    # Generate metadata from fasta file
    def generate_chunks(self, key, id, chunk_size, obj_size):
        min_range = id * chunk_size
        max_range = (id + 1) * chunk_size
        data = self.storage.get_object(bucket=self.bucket, key=key,extra_get_args={'Range': f'bytes={min_range}-{max_range}'}).decode('utf-8')
        #data = obj.data_stream.sb.read().decode('utf-8')
        '''FASTA_HEADER_START = r">"
        FASTA_HEADER = r">.+\n"
        content = []
        found = False
        treated_first_sequence = False
        ini_heads = list(re.finditer(FASTA_HEADER_START, data))
        heads = list(re.finditer(FASTA_HEADER, data))
        i = 0
        start = 0
        prev = ""
        aux = []
        # TODO  : Mirar si hacer loops en base a cual puede que a través de '>' pero ir incrementando un indice paa acceder a
        # TODO  | '>.+\n' (teniendo en cuenta el tamaño de la lista len(heads))
        for m in ini_heads:
            found = True
            start = min_range + m.start()
            end = min_range + m.end()
            if start < (obj.part + 1) * obj.chunk_size: # Check that the initial byte of the head is within the range
                if not treated_first_sequence:
                    treated_first_sequence = True
                    if obj.part >= 1 and start > obj.part * obj.chunk_size:  # If it is not the worker of the first part of the file and in addition it
                                                                             # turns out that the partition begins in the middle of the data of a sequence

                        # content.append(data[0:m.start()])
                        # >> num_chunks_has_divided offset_head offset_bases_split
                       content.append('>>' + ' ' + 'X' + ' ' + 'Y' + ' ' + str(obj.part * obj.chunk_size))  # Split sequences
                # ------------------------------------------------
                # TODO saber que fa i si es útil consevar-ho
                if (obj.part == int(total_obj_size)) and (int(obj_size) - int(start + len(
                        m.group()))) < total_obj_size:  # Si es tracta del worker que te l'última part del fitxer
                    # original i ?????
                    content.append('This should not be')
                    #content = content + str(start) + ',' + str(obj_size) + 'x\n'  # Head of the current
                # ------------------------------------------------
                else:  # When it is not the worker who has the last part of the original file
                    if str(prev) != str(start) and prev != "":  # When if there is a previous sequence already read and is not
                                                                # identical to the current one and there is a previous sequence read
                        # name_id num_chunks_has_divided offset_head offset_bases
                        content.append(m.group().split(' ')[0].replace('>','') + ' ' + '1' + ' ' + str(start) + ' ' + str(end))

                    prev = str(end)

        if last_{end >= (obj.part + 1) * obj.chunk_size and found:  # Check if the last head of the current one is cut
            # >->name_id num_chunks_has_divided offset_head offset_bases
            content.append(f">->{m.group().split(' ')[0]}" + ' ' + '1' + ' ' + str(start) + ' ' + str(end))
            # content.append(str(end) + ',' + data[(m.end()::])


        if not found:
            content.append('<Sequencies not found>')


        dict = {'min_range': str(min_range), 'max_range': str(obj_size) if obj.part == int(total_obj_size) else str(max_range),
                'sequences': content}
        self.storage.put_object(self.my_bucket_name, f'cache/chunk_{self.obj.part}_index.json', str(dict))'''

        return [str(id), str(min_range), str(max_range), str(len(data)), chunk_size, obj_size]


def run_worker_metadata(key, id, storage, bucket, chunk_size, obj_size):
    partitioner = PartitionFasta(storage, bucket, key, id)
    dades = partitioner.generate_chunks(key, id, chunk_size, obj_size)
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

                elif '>->' in dict['sequences'][0]: # Comprovar que el primer index no porta >>:
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
    data_bucket = 'cloud-fasta-partitioner'  # Change-me
    prefix = 'fasta/'  # Change-me
    local_input_path = './input_data/'  # Change-me

    key = f'fasta/genes.fasta'  # Change-me
    obj = f'cos://{data_bucket}/{key}'  # f'{my_bucket_name}/{path_obj[0]}'  # Change-me
    print(obj)

    workers = 20  # Change-me

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
    extra_args = {'bucket': data_bucket,'chunk_size': chunk_size, 'obj_size': fasta['content-length']}
    fexec.map_reduce(map_function=run_worker_metadata,map_iterdata=map_iterdata, extra_args=extra_args, reduce_function=reduce_generate_chunks)
    fexec.wait()
    results = fexec.get_result()

    #generate_json(results, f'./output_data/{pathlib.Path(elem).stem}_index')
    print(results)
    for el in results:
        print(el)

    fexec.plot()
    fexec.clean()
