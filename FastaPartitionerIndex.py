# from pyfaidx import Fasta
# genes = Fasta('input_data/genes.fasta')
import pathlib
import re
import os
import lithops
import json
import json
import gc


class PartitionFasta():

    def __init__(self, storage, my_bucket_name, obj):
        self.storage = storage
        self.my_bucket_name = my_bucket_name
        self.obj = obj

    # Generate metadata from fasta file
    def generate_chunks(self, obj, total_obj_size, obj_size, overlap):
        # JSON con los rangos y dentro todas las secuencias??
        # x = {'range_min': '0', 'range_max': '50', secu: '[ksndksnsdk,dnejnde,eoeifjriofjeif]'}

        ## m.group().replace('>', '')

        data = obj.data_stream.sb.read().decode('utf-8')
        FASTA_HEADER = r">.+\n"
        content = []
        found = False
        titles = list(re.finditer(FASTA_HEADER, data))
        i = 0
        start = 0
        last_el = 0
        prev = ""
        for m in titles:
            found = True
            start = obj.data_byte_range[0] + m.start()
            if i == 0 and obj.part >= 1 and start > obj.part * obj.chunk_size:  # If it is not the worker of the first part of the file and in addition it
                                                                                # turns out that the partition begins in the middle of the data of a sequence
                # >> num_chunks_has_divided offset_head offset_bases_split
                content.append('>>' + ' ' + 'X' + ' ' + 'Y' + ' ' + str(obj.part * obj.chunk_size))  # Split sequence
            # ------------------------------------------------
            # TODO saber que fa i si es útil consevar-ho
            if (obj.part == int(total_obj_size)) and (int(obj_size) - int(start + len(
                    m.group()))) < total_obj_size:  # Si es tracta del worker que te l'última part del fitxer
                # original i ?????
                content = content + str(start) + ',' + str(obj_size) + 'x\n'  # Head of the current
            # ------------------------------------------------
            else:  # When it is not the worker who has the last part of the original file
                if str(prev) != str(start) and prev != "":  # When if there is a previous sequence already read and is not
                                                            # identical to the current one and there is a previous sequence read
                    # name_id num_chunks_has_divided offset_head offset_bases
                    content.append(m.split(' ')[0] + ' ' + '1' + ' ' + str(start) + ' ' + str(start + len(m.group())))

                #prev = str(start + len(m.group()))

            last_el = start + len(m.group())
            i = +1
        # TODO Ja no: Falta obtener la longitud de los datos del ultimo elemento
        if last_el < (obj.part + 1) * obj.chunk_size and found:  # Check if the last head of the current one is not cut
            if obj.part == int(total_obj_size):  # If it is the worker who has the last part of the original file
                content = content + str(last_el) + ',' + str(int(obj_size)) + '\n'  # Sequence of the current
            else:  # If it is not the worker who has the last part of the original file
                content = content + str(last_el) + ',' + str(
                    str((obj.part + 1) * obj.chunk_size + overlap)) + '\n'  # Sequence of the current
        else:   # The header is cut
            # >>>name_id num_chunks_has_divided offset_head offset_bases
            content.append(f"{m}" + ' ' + '1' + ' ' + str(start) + ' ' + str(start + len(m.group())))



        if not found:
            if obj.part == int(total_obj_size):
                content = content + str(obj.data_byte_range[0]) + ',' + str(obj_size) + '\n'
            else:
                content = content + str(obj.data_byte_range[0]) + ',' + str(
                    ((obj.part + 1) * obj.chunk_size) + overlap) + '\n'

        self.storage.put_object(self.my_bucket_name, f'cache/obj{obj.part}.data', content)

        return {'min_range': str(obj.data_byte_range[0]), 'max_range': str(obj.data_byte_range[1]),
                'secuences': content}


def reduce_generate_chunks():
    # Trobar si hi ha algú worker que ha començat amb una secuencia partida i dir-li a quina secuencia pertany (>>
    # num_particions offset -> nam_id  num_particions offset, i actualitzar num_particions de la secuencia que si
    # tenia la part del head (0 -> x))
    # S'ha d'actualitzar el rang en cas que just en un worker s'ha quedat el head (o parcialment) i un altre worker té tota la secuencia (i potser part del head)
    # Comprovar només mirant la primera secuencia, si esta dividida, comprovar l'ultima secuencia de l'anterior chunk:
    #     si té '>>>name_id': passar el head al següent chuhnk (i eliminar els simbols '>>>') i actualitzar rangs
    #     si té '>>: copiar el head a l'index i donar-li el offset del head
    return None


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


def get_objext_funct(data_bucket, elem, dt_dir):
    data = json.load(storage.get_object(bucket=data_bucket, key=elem))
    with open(f'{dt_dir}.json', 'w') as f:
        f.write(
            '[' +
            ',\n'.join(json.dumps(i) for i in dict) +
            ']\n')


# s3://sra-pub-src-2/SRR8774337/cleaned.fasta

if __name__ == "__main__":
    # Param
    my_bucket_name = 'cloud-fasta-partitioner'  # Change-me
    local_input_path = './input_data/'  # Change-me
    prefix = 'fasta/'
    elem = f'{prefix}genes.fasta'

    # Execution
    fexec = lithops.FunctionExecutor(log_level='DEBUG', max_workers=1000)
    storage = lithops.Storage()

    # obj = 's3://ayman-lithops-meta-cloudbutton-hutton/fasta/hg19.fa'  # Change-me

    path_obj = push_object_funct(local_input_path, my_bucket_name, prefix)
    print(path_obj)
    obj = f'cos://{my_bucket_name}/{elem}'  # f'{my_bucket_name}/{path_obj[0]}'
    print(obj)

    location = obj.split('//')[1].split('/')  # obj.split('/') #
    for_head = location[1:]
    for_head = '/'.join(for_head)
    data_bucket_name = location[0]
    ovlp = 300
    worker_chunk_size = 50000000  # 500000000
    fasta = storage.head_object(data_bucket_name, for_head)
    total_fasta_size = int(fasta['content-length']) / worker_chunk_size
    seq_name = location[-1].split('.')[0]

    print('===================================================================================')
    print('metadata chunks: ' + str(int(total_fasta_size)))
    print('bucket to access data: ' + str(data_bucket_name))
    print('reference genome name: ' + str(seq_name))
    print('fasta size: ' + str(fasta['content-length']) + ' bytes')
    print('===================================================================================')

    # map_reduce
    fexec.map(run_worker_metadata, {'my_bucket_name': my_bucket_name, 'obj': obj, 'total_obj_size': total_fasta_size,
                                    'obj_size': fasta['content-length']}, obj_chunk_size=worker_chunk_size)
    fexec.wait()
    results = fexec.get_result()

    get_objext_funct(my_bucket_name, elem, f'./output_data/{pathlib.Path(elem).stem}_index')

    print(results)

    for el in results:
        print(el)
        '''for k in el.split('\n'):
            print(k)
            if 'x' in k:
                print(k)'''

    fexec.plot()
    fexec.clean()
