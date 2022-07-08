import os
import pathlib
import re
import lithops
import json
import gc
import time


class PartitionFasta():

    def __init__(self, storage, my_bucket_name, obj):
        self.storage = storage
        self.my_bucket_name = my_bucket_name
        self.obj = obj

    # Generate metadata from fasta file
    def generate_chunks(self, obj, total_obj_size, obj_size, overlap):
        data = obj.data_stream.sb.read().decode('utf-8')
        '''FASTA_HEADER = r">.+\n"
        content = ""
        found = False
        titles = list(re.finditer(FASTA_HEADER, data))
        i = 0
        start = 0
        last_el = 0
        prev = ""
        auxi = ""
        for m in titles:
            found = True
            start = obj.data_byte_range[0] + m.start()
            if i == 0:
                aux = f'm: {m}\nm.start(): {m.start()}\nm.gorup() len : {len(m.group())}\nm.gorup() : {m.group()}\n\n'
                if obj.part >= 1:
                    if start > (
                            obj.part) * obj.chunk_size:  # Si no es tracta del worker de la primera part del fixer i a més resulta que la particio comença a la meitat
                        # de la data d'una secuencia , entra
                        content = content + str((obj.part) * obj.chunk_size) + ',' + str(start + overlap) + '\n'
            if (obj.part == int(total_obj_size)) and (int(obj_size) - int(start + len(
                    m.group()))) < total_obj_size:  # Si es tracta del worker que te l'última part del fitxer original i ?????
                auxi = str(int(obj_size) - int(start + len(m.group()))) + ' - '
                content = content + str(start) + ',' + str(obj_size) + 'x\n'  # Head de l'actual
            else:  # Quan no es tracta de l'última
                if str(prev) != str(
                        start) and prev != "":  # Quan si hi a secuencia previa i no és identica a la actual, Entra
                    # auxi = 'else if: ' + content + prev +',' + str(start)+ ' 1\n'
                    content = content + 'secuencie data: ' + prev + ',' + str(
                        start) + '\n'  # Secuencia de l'anterior (prev)
                    content = content + 'head: ' + str(start) + ',' + str(
                        start + len(m.group())) + 'x\n'  # Head de l'actual
                else:  # Quan comença i encara no hi ha anterior secuencia, Entra
                    # auxi = 'else else: ' + content + str(start) + ',' + str(start + len(m.group())) + ' x2\n'
                    content = content + 'head: ' + str(start) + ',' + str(
                        start + len(m.group())) + 'x\n'  # Head de l'actual
                prev = str(start + len(m.group()))

            last_el = start + len(m.group())
            i = +1

        if last_el < (obj.part + 1) * obj.chunk_size and found:  # Check if the last head of the current one is not cut
            if obj.part == int(total_obj_size):  # If it is the worker who has the last part of the original file
                content = content + str(last_el) + ',' + str(int(obj_size)) + '\n'  # Sequence of the current
            else:  # If it is not the worker who has the last part of the original file
                content = content + str(last_el) + ',' + str(
                    str((obj.part + 1) * obj.chunk_size + overlap)) + '\n'  # Sequence of the current

        if not found:
            if obj.part == int(total_obj_size):
                content = content + str(obj.data_byte_range[0]) + ',' + str(obj_size) + '\n'
            else:
                content = content + str(obj.data_byte_range[0]) + ',' + str(
                    ((obj.part + 1) * obj.chunk_size) + overlap) + '\n'

        self.storage.put_object(self.my_bucket_name, f'cache/obj{obj.part}.data', content)'''

        return [[str(obj.data_byte_range[0]), str((obj.part) * obj.chunk_size)], [str(obj.data_byte_range[1]), str((obj.part + 1) * obj.chunk_size), str(len(data) + obj.data_byte_range[0])], obj_size, str(obj.chunk_size - (obj.data_byte_range[1] - obj.data_byte_range[0]))]  # data.partition('\n')[0] #auxi # content #f'cache/obj{obj.part}.data'


def print_chunk_info(obj):
    # obj is a CloudObject
    aux = f"===================\n" + 'obj: ' + str(obj) + "\n" + 'key: ' + str(obj.key) + "\n" + 'part: ' + str(
        obj.part) + "\n" + 'range: ' + str(obj.data_byte_range) + "\n" + 'chunk_size: ' + str(
        obj.chunk_size) + "\n===================\n"
    return aux
    '''print("===================")
    print(obj)
    print(obj.key)
    print('part:'+str(obj.part))
    print(obj.data_byte_range)  
    print(obj.chunk_size)  # in bytes
    print("===================")'''


def run_worker_metadata(obj, storage, my_bucket_name, total_obj_size, obj_size):
    aux = print_chunk_info(obj)
    # return aux
    partitioner = PartitionFasta(storage, my_bucket_name, obj)
    dades = partitioner.generate_chunks(obj, total_obj_size, obj_size, 0)
    return dades


def push_object_funct(input_dir, data_bucket, input_data_prefix):
    bucket_objects = storage.list_keys(bucket=data_bucket)
    # storage.delete_objects(data_bucket,bucket_objects)
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

            '''if key not in bucket_objects:  # Changed: if file_name not in bucket_objects:
                with open(os.path.join(subdir, file_name), 'rb') as file:  # Changed
                    print(f'\tUploading {key}...')
                    data = file.read()
                    storage.put_object(bucket=data_bucket, key=key + '_1', body=data)
                    print('\tOk!')
            else:  # Added
                print(f'\tIt is already uploaded: {key}...')  # Added'''
            print("Finished!")

    return storage.list_keys(bucket=data_bucket)


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

    print(pathlib.Path('fasta/genes.fasta').name)
    print(pathlib.Path('fasta/genes.fasta').stem)

    location = obj.split('//')[1].split('/')  # obj.split('/') #
    for_head = location[1:]
    for_head = '/'.join(for_head)
    data_bucket_name = location[0]
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

    fexec.map(run_worker_metadata, {'my_bucket_name': my_bucket_name, 'obj': obj, 'total_obj_size': total_fasta_size,
                                    'obj_size': fasta['content-length']}, obj_chunk_size=worker_chunk_size)
    fexec.wait()
    results = fexec.get_result()
    print(results)

    for el in results:
        print(el)
        '''for k in el.split('\n'):
            print(k)
            if 'x' in k:
                print(k)'''

    # fexec.plot()
    # fexec.clean()
