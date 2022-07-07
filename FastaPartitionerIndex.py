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

    def __init__(self, storage, my_bucket_name, obj):
        self.storage = storage
        self.my_bucket_name = my_bucket_name
        self.obj = obj

    # Generate metadata from fasta file
    def generate_chunks(self, obj, total_obj_size, obj_size):
        # JSON con los rangos y dentro todas las secuencias??
        # x = {'range_min': '0', 'range_max': '50', secu: '[ksndksnsdk,dnejnde,eoeifjriofjeif]'}

        ## m.group().replace('>', '')

        data = obj.data_stream.sb.read().decode('utf-8')
        FASTA_HEADER = r">.+\n"
        content = []
        found = False
        treated_first_sequence = False
        titles = list(re.finditer(FASTA_HEADER, data))
        i = 0
        start = 0
        last_el = 0
        prev = ""
        aux = []
        for m in titles:
            found = True
            start = obj.data_byte_range[0] + m.start()
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
                        content.append(m.group().split(' ')[0].replace('>','') + ' ' + '1' + ' ' + str(start) + ' ' + str(start + len(m.group())))

                    prev = str(start + len(m.group()))

                last_el = [start + len(m.group()), m.group()]

        if last_el[0] >= (obj.part + 1) * obj.chunk_size and found:  # Check if the last head of the current one is cut
            # >->name_id num_chunks_has_divided offset_head offset_bases
            content.append(f">->{last_el[1].split(' ')[0]}" + ' ' + '1' + ' ' + str(start) + ' ' + str(start + len(m.group())))
            # content.append(str(last_el[0])+','+data[(m.start()+len(m.group()))::])


        if not found:
            if obj.part == int(total_obj_size): # If it is the worker who has the last part of the original file
                content.append(str(obj.data_byte_range[0]) + '<->' + str(obj_size))
            else:  # If it is not the worker who has the last part of the original file
                content.append(str(obj.data_byte_range[0]) + '<->' + str((obj.part + 1) * obj.chunk_size))

        dict = {'min_range': str(obj.data_byte_range[0]), 'max_range': str(obj_size) if obj.part == int(total_obj_size) else str((obj.part + 1) * obj.chunk_size),
                'secuences': content}
        self.storage.put_object(self.my_bucket_name, f'cache/chunk_{obj.part}_index.json', str(dict))

        return dict


def print_chunk_info(obj):
    # obj is a CloudObject
    print("===================")
    print(obj)
    print(obj.key)
    print('part:' + str(obj.part))
    print(obj.data_byte_range)
    print(obj.chunk_size)  # in bytes
    print("===================")


def run_worker_metadata(obj, storage, my_bucket_name, total_obj_size, obj_size):
    print_chunk_info(obj)
    partitioner = PartitionFasta(storage, my_bucket_name, obj)
    dades = partitioner.generate_chunks(obj, total_obj_size, obj_size)
    return dades


def reduce_generate_chunks(results):
    # · Trobar si hi ha algú worker que ha començat amb una secuencia partida i dir-li a quina secuencia pertany (>>
    #   num_particions offset -> nam_id  num_particions offset, i actualitzar num_particions de la secuencia que si
    #   tenia la part del head (0 -> x))
    # · S'ha d'actualitzar el rang en cas que just en un worker s'ha quedat el head (o parcialment) i un altre worker té tota la secuencia (i potser part del head)
    # · Comprovar només mirant la primera secuencia, si esta dividida, comprovar l'ultima secuencia de l'anterior chunk:
    #     - si té '>>>name_id': passar el head al següent chuhnk (i eliminar els simbols '>>>') i actualitzar rangs
    #     - si té '>>: copiar el head a l'index i donar-li el offset del head
    if len(results) > 1:
        i = 0
        for dict in results:
            if i > 0:
                if '>>' in dict['secuences'][0]:  # Comprovar que el primer index no porta >>
                    leng_prev = len(results[i-1]['secuences'])
                    sec_prev = results[i-1]['secuences'][leng_prev - 1]
                    dict['secuences'][0] = dict['secuences'][0].replace('Y', sec_prev.split(' ')[2])   # Y -> offset_head
                    # Comprovar que no hagi més d'una partició en una mateixa secuencia (s'ha d'actualitzar el nombre de particions de tots els afectats)
                    if int(results[i-1]['secuences'][leng_prev - 1].split(' ')[1]) > 1: # If there are more than one split in the sequence
                        for x in range(0, i + 1):
                            dict['secuences'][0] = dict['secuences'][0].replace(i, i + 1)  # num_chunks_has_divided + 1 (i+1: total of current partitions of sequence)
                    else: # If sequences is only split into two parts
                        results[i - 1]['secuences'][leng_prev - 1].replace('1', i + 1)
                        dict['secuences'][0] = dict['secuences'][0].replace('X', i + 1)
                    dict['secuences'][0] = dict['secuences'][0].replace('>>', sec_prev.split(' ')[0])  # '>>' -> name_id

                    # TODO ->
                    # Actualizar rangos

                elif '>->' in dict['secuences'][0]: # Comprovar que el primer index no porta >>:
                    print('')
            i += 1
    else:
        final_results = results
    return final_results # results


def push_object_funct(input_dir, data_bucket, input_data_prefix):
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
    my_bucket_name = 'cloud-fasta-partitioner'  # Change-me
    local_input_path = './input_data/'  # Change-me
    prefix = 'fasta/'  # Change-me
    elem = f'{prefix}genes.fasta'  # Change-me

    obj = f'cos://{my_bucket_name}/{elem}'  # f'{my_bucket_name}/{path_obj[0]}'  # Change-me

    # Execution
    fexec = lithops.FunctionExecutor(log_level='DEBUG', max_workers=1000, runtime_memory=4096)
    storage = lithops.Storage()



    path_obj = push_object_funct(local_input_path, my_bucket_name, prefix)
    print(path_obj)
    print(obj)

    location = obj.split('//')[1].split('/')  # obj.split('/') #
    for_head = location[1:]
    for_head = '/'.join(for_head)
    data_bucket_name = location[0]
    ovlp = 300
    worker_chunk_size = 500000000  # 500000000
    fasta = storage.head_object(data_bucket_name, for_head)
    total_fasta_size = int(fasta['content-length']) / worker_chunk_size
    seq_name = location[-1].split('.')[0]

    print('===================================================================================')
    print('metadata chunks: ' + str(int(total_fasta_size)))
    print('bucket to access data: ' + str(data_bucket_name))
    print('reference genome name: ' + str(seq_name))
    print('fasta size: ' + str(fasta['content-length']) + ' bytes')
    print('===================================================================================')

    map_iterdata = {'my_bucket_name': my_bucket_name, 'obj': obj, 'total_obj_size': total_fasta_size, 'obj_size': fasta['content-length']}
    fexec.map_reduce(map_function=run_worker_metadata,map_iterdata=map_iterdata, reduce_function=reduce_generate_chunks, obj_chunk_size=worker_chunk_size)
    fexec.wait()
    results = fexec.get_result()

    #generate_json(results, f'./output_data/{pathlib.Path(elem).stem}_index')
    print(results)
    '''for el in results:
        print(el['min_range'])'''

    fexec.plot()
    fexec.clean()
