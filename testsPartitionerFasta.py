import unittest
import fastaPartitionerIndex as fp
import random

result = []

#TODO: actualizar test con el nuevo formato
class TestPartitionOptions(unittest.TestCase):
    def setUp(self):
        self.path_data = './output_data/genes_index.json'
        self.path_data_assert = './input_data/genes.fasta.fai'

    def test_get_info_sequence(self):
        with open(self.path_data_assert, 'r') as d:
            seq_d_a = d.readline()
        # Not updated
        i = 0
        list = []
        for el in self.data_assert:
            if random.randint(0,4) == 0:
                list.append(el.split('\t'))
                if i > 500:
                    break
                i += 1
        list.extend([['tr|IAmFalse'], ['']])

        for i, el in enumerate(list):
            info = fp.get_info_sequence(self.data, el[0])
            if info['length'] != -1 and info['offset'] != -1 and info['offset_head'] != -1:
                self.assertEqual([el[0], info['length'], info['offset']], [el[0], int(el[1]), int(el[2])])

    def test_get_range_sequence(self):
        with open(self.path_data_assert, 'r') as d:
            seq_d = d.readline()
        # Not updated
        #list = [{'min': dict['min_range'], 'max': dict['max_range']} for dict in self.data]
        #range = list[random.randint(0, len(list))]
        result.append(f"Test 'test_get_range_sequence'")
        for el in self.data:
            sequences = fp.get_sequences_of_range(self.data, el['min_range'], el['max_range'])
            self.assertEqual(len(sequences), len(el['sequences']))
            result.append(f"{el['min_range']}-{el['max_range']}: {sequences[0]}")

    def test_index_generated(self):
        j = 0
        last_seq = ''
        for i, dict in enumerate(self.data):
            for sequence in dict['sequences']:
                info = sequence.split(' ')
                if last_seq != info[0]:
                    split = int(info[1])
                    if split > 1:
                        length = 0
                        for x in range(i + 1, i + split):
                            length += int(self.data[x]['sequences'][0].split(' ')[4])
                        length += int(info[4])
                        info[4] = str(length)
                    info_assert = self.data_assert[j].split('\t')
                    self.assertEqual([info[0], info[4], info[3]], [info_assert[0],info_assert[1],info_assert[2]])

                    last_seq = info[0]
                    j += 1

if __name__ == '__main__':
    unittest.main()
