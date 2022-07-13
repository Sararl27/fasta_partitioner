import re
import unittest

from pyfaidx import Fasta

import lithops
import os
import fastaPartitionerIndex as fp
import testsPartitionerFasta


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


def generate_fasta_index_pyfaidx():
    genes = Fasta('input_data/genes.fasta')
    with open('./input_data/genes.fasta.fai', "w") as f:
        f.write(genes)


if __name__ == "__main__":
    storage = lithops.Storage()

    # Params
    data_bucket = 'cloud-fasta-partitioner'  # Change me
    prefix = 'fasta/'  # Change me
    local_input_path = './input_data/genes.fasta'  # Change me

    key = f'fasta/genes.fasta'  # Change me
    obj = f'{data_bucket}/{key}'  # f'cos://{data_bucket}/{key}'  # Change me

    workers = 100  # Change me

    # Execution
    #push_object_funct(local_input_path, data_bucket, prefix)
    #generate_fasta_index_pyfaidx()
    #fp.FastaPartitioner(storage, data_bucket, key, workers)

    # run all tests with verbosity
    unittest.TextTestRunner().run(unittest.TestLoader().loadTestsFromModule(testsPartitionerFasta))
    # TODO descubrir pq da diferente offset de base entre ['tr|A0A068S3P6|A0A068S3P6_9FUNG', '1054', '1191018'] y ['tr|A0A068S3P6|A0A068S3P6_9FUNG', '1054', '1191156']