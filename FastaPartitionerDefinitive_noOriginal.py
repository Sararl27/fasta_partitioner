import os
import re
import lithops
import json
import gc
import time
import json

dict = [{
  "id": "04",
  "name": "sunil",
  "department": "HR"
} , {
  "id": "03",
  "name": "sunil33",
  "department": "H3R"
} ]

with open(f'./output_data/test.json', 'w') as f:
    #json.dump(dict, f)
    f.write(
      '[' +
      ',\n'.join(json.dumps(i) for i in dict) +
      ']\n')

with open('./output_data/test.json') as f:
  data = json.load(f)

print(data)
print(data[0]['id'])

print()

characters = ['Tokyo', 'Lisbon', 'Moscow', 'Berlin']
new_characters = ['Nairobi', 'Denver', 'Rio']
characters.append(new_characters)
print('"append" list:', characters)
characters = ['Tokyo', 'Lisbon', 'Moscow', 'Berlin']
characters.extend(new_characters)
print('"extend" list:', characters)

print()

head_sec = '>>>tr|A0A7X8NEN1|A0A7X8NEN1_9ACTN>>> Glutamate racemase OS=Olsenella sp. KGMB02461 OX=2726201 GN=murI PE=3 SV=1'
print(f'Primera part de {head_sec}:\n\t'
      + head_sec.split(' ')[0])
print(f'Tret > de {head_sec}:\n\t'
      + head_sec.split(' ')[0].replace('>',''))

print()

print(f'There is >>> in >>{head_sec}?')
if '>>>' in f'>>{head_sec}':
   print('success')
else:
    print('no')

results = [{"min_range": "0",
    "max_range": "119101845",
    "secuences": ["tr|A0A024R230|A0A024R230_HUMAN 1 2430 2538","tr|A0A024R244|A0A024R244_HUMAN 1 3390 3487","tr|A0A024R328|A0A024R328_HUMAN 1 4133 4236",
      "tr|A0A024R6R4|A0A024R6R4_HUMAN 1 4924 5025",
      "tr|A0A024R7T2|A0A024R7T2_HUMAN 1 5696 5821",
      "tr|A0A024R809|A0A024R809_HUMAN 1 6954 7053",
      "tr|A0A024R813|A0A024R813_HUMAN 1 7939 8033",
      "sp|A0A024RBG1|NUD4B_HUMAN 1 8635 8762",
      "tr|A0A024RC06|A0A024RC06_HUMAN 1 8947 9034",
      "tr|A0A061K6E7|A0A061K6E7_ECOLX 1 9648 9768",
      "tr|A0A061QB49|A0A061QB49_9VIBR 1 10510 10638"]}, {"min_range": "0",
    "max_range": "119101845",
    "secuences": ["tr|A0A024R230|A0A024R230_HUMAN 1 2430 2538","tr|A0A024R244|A0A024R244_HUMAN 1 3390 3487","tr|A0A024R328|A0A024R328_HUMAN 1 4133 4236",
      "tr|A0A024R6R4|A0A024R6R4_HUMAN 1 4924 5025",
      "tr|A0A024R7T2|A0A024R7T2_HUMAN 1 5696 5821",
      "tr|A0A024R809|A0A024R809_HUMAN 1 6954 7053",
      "tr|A0A024R813|A0A024R813_HUMAN 1 7939 8033",
      "sp|A0A024RBG1|NUD4B_HUMAN 1 8635 8762",
      "tr|A0A024RC06|A0A024RC06_HUMAN 1 8947 9034",
      "tr|A0A061K6E7|A0A061K6E7_ECOLX 1 9648 9768",
      "tr|A0A061QB49|A0A061QB49_9VIBR 1 10510 10638"]}]

x = 1
for dict in results:
    dict['secuences'][x] = 'a'   # Y -> offset_head
    print(dict['secuences'][x])
    x += 1

print(results)



'''class PartitionFasta():

        


    def __init__(self, storage, my_bucket_name,obj):
        self.storage = storage
        self.my_bucket_name = my_bucket_name
        self.obj = obj
        

    #Generate metadata from fasta file
    def generate_chunks(self, obj, obje,key,total_obj_size,obj_size,overlap):
        data = storage.get_object(obje, key).decode('utf-8')
        # data = obj.data_stream.sb.read().decode('utf-8')
        # with open(obj, "r", encoding='utf-8') as f:
        #    data = f.read()
        FASTA_HEADER = r">.+\n"
        content = ""
        found = False
        titles = list(re.finditer(FASTA_HEADER,data))
        i = 0
        start = 0
        last_el = 0
        prev = ""
        for m in titles:
            found = True
            start = obj.data_byte_range[0] + m.start()
            if i == 0:
                if obj.part >= 1:
                    if start > (obj.part)*obj.chunk_size: 
                        content =  content +  str((obj.part)*obj.chunk_size) +',' + str(start+overlap)+ '\n'
            if (obj.part == int(total_obj_size)) and  (int(obj_size) - int(start + len(m.group()))) < total_obj_size:
                content =  content + str(start) +',' + str(obj_size)+ 'x\n'
            else:
                if str(prev) != str(start) and prev != "":
                    content =  content + prev +',' + str(start)+ '\n'
                    content =  content + str(start) +',' + str(start + len(m.group()))+ 'x\n'
                else:
                    content =  content + str(start) +',' + str(start + len(m.group()))+ 'x\n'
                prev = str(start + len(m.group()))

                
                
            last_el = start + len(m.group())
            i=+1
        
        if last_el < (obj.part+1)*obj.chunk_size and found:
            if obj.part == int(total_obj_size):
                content = content + str(last_el) +','+ str(int(obj_size))+ '\n'
            else:
                content = content + str(last_el) +','+ str(str((obj.part+1)*obj.chunk_size + overlap))+ '\n'

        if not found :
            if obj.part == int(total_obj_size):
                content = content + str(obj.data_byte_range[0]) +',' + str(obj_size)+ '\n'
            else:
                content = content + str(obj.data_byte_range[0]) +','+ str(((obj.part+1)*obj.chunk_size) + overlap)+ '\n'
        

        self.storage.put_object(self.my_bucket_name, f'cache/obj{obj.part}.data',content)
        
        return f'cache/obj{obj.part}.data'
       


def print_chunk_info(obj):
    # obj is a CloudObject
    print("===================")
    print(obj)
    print(obj.key)
    print('part:'+str(obj.part))
    print(obj.data_byte_range)  
    print(obj.chunk_size)  # in bytes
    print("===================")

def run_worker_metadata(my_bucket_name, key, total_obj_size,obj_size):
    #print_chunk_info(obj)
    partitioner = PartitionFasta(storage,my_bucket_name,obj)
    # dades = partitioner.generate_chunks(obj, total_obj_size, obj_size, '')
    dades = partitioner.generate_chunks(obj, my_bucket_name, key, total_obj_size, obj_size, '')
    return dades


def push_object_funct(input_dir, data_bucket, input_data_prefix):
    bucket_objects = storage.list_keys(bucket=data_bucket)
    for subdir, dirs, files in os.walk(local_input_path):
        print(subdir)
        for file_name in files:
            key = os.path.join(input_data_prefix, file_name)  # Added
            if key not in bucket_objects:   # Changed: if file_name not in bucket_objects:
                with open(os.path.join(subdir, file_name), 'rb') as file: #Changed
                    print(f'\tUploading {key}...')
                    data = file.read()
                    storage.put_object(bucket=data_bucket, key=key, body=data)
                    print('\tOk!')
            else:   # Added
                print(f'\tIt is already uploaded: {key}...')   # Added
    print("Finished!")

    return storage.list_keys(bucket=data_bucket)


    
#s3://sra-pub-src-2/SRR8774337/cleaned.fasta

if __name__ == "__main__":
    storage = lithops.Storage()
    # obj = 's3://ayman-lithops-meta-cloudbutton-hutton/fasta/hg19.fa'  # Change-me
    my_bucket_name = 'cloudbutton-hutton' # Change-me
    local_input_path = './input_data/' # Change-me
    path_obj = push_object_funct(local_input_path, my_bucket_name, 'fasta')
    print(path_obj[0])
    obj = f'cos://{my_bucket_name}/{path_obj[0]}'

    fexec = lithops.LocalhostExecutor()

    location = obj.split('//')[1].split('/')
    for_head = location[1:]
    for_head = '/'.join(for_head)
    data_bucket_name = location[0]
    ovlp = 300
    storage = lithops.Storage()
    worker_chunk_size = 500000000
    fasta = storage.head_object(data_bucket_name, for_head)
    total_fasta_size = int(fasta['content-length']) / worker_chunk_size
    seq_name = location[-1].split('.')[0]
    print('===================================================================================')
    print('# metadata chunks: ' + str(int(total_fasta_size)))
    print('bucket to access data: ' + str(data_bucket_name))
    print('reference genome name: ' + str(seq_name))
    print('fasta size: ' + str(fasta['content-length']) + ' bytes')
    print('===================================================================================')

    # local = './input_data/genes.fasta'
    run_worker_metadata(my_bucket_name, path_obj[0], total_fasta_size, fasta['content-length'])
    fexec.map(run_worker_metadata, {'my_bucket_name':my_bucket_name,'obj':obj,'total_obj_size':total_fasta_size,'obj_size':fasta['content-length']},
            obj_chunk_size=worker_chunk_size)
    
    fexec.wait()
    results = fexec.get_result()

    for el in results:
        for k in el.splitlines():
            if 'x' in k:
                print(k)

    fexec.plot()
    fexec.clean()'''
