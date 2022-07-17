import unittest
import random
import json
import re
import fastaPartitionerIndex as fp
import random

results = []
class TestPartitionOptions(unittest.TestCase):
    def setUp(self):
        with open('./output_data/genes_index.json', "r") as f:
            self.data = json.load(f)
        with open('./input_data/genes.fasta.fai', "r") as f:
            self.data_assert = f.read().splitlines()

    def test_get_info_sequence(self):
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
        results.append(f"Test 'test_get_range_sequence'")
        for el in self.data:
            sequences = fp.get_sequences_of_range(self.data, el['min_range'], el['max_range'])
            self.assertEqual(len(sequences), len(el['sequences']))
            #results.append(f"{el['min_range']}-{el['max_range']}: {sequences[0]}")
        list = [[9643, 36109], [60411, 70827], [113494, 120950], [473060, 717220], [949179, 957690]]
        for i in list:
            sequences = fp.get_sequences_of_range(self.data, i[0], i[1])
            results.append(f"{i[0]}-{i[1]}: {sequences[0]}\t\t{sequences[-1]}")
            results.append(str(sequences) + '\n')
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
