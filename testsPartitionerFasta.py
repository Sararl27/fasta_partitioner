import unittest
import random
import json
import re
import fastaPartitionerIndex as fp

class TestPartitionOptions(unittest.TestCase):
    def setUp(self):
        with open('./output_data/genes_index.json', "r") as f:
            self.data = json.load(f)
        with open('./input_data/genes.fasta.fai', "r") as f:
            self.data_assert = f.read().splitlines()

    def test_get_info_sequence(self):
        list = ['tr|IAmFalse', 'tr|A0A1W6BZY2|A0A1W6BZY2_ECOLX', 'tr|A0A7H0NQ17|A0A7H0NQ17_9ACTN','tr|A0A8H7Z4U6|A0A8H7Z4U6_AJECA',
                'tr|A0A8H7Z4X4|A0A8H7Z4X4_AJECA', 'tr|A0A8H7Z523|A0A8H7Z523_AJECA','' , 'tr|Z4YJ09|Z4YJ09_DANRE',
                'tr|A0A8D5HB68|A0A8D5HB68_RHOHA', 'tr|G4V624|G4V624_SCHMA', 'tr|A0A417FG44|A0A417FG44_9FIRM',
                'tr|A0A0F7CYC8|A0A0F7CYC8_HHV1']
        for el in list:
            info = fp.get_info_sequence(self.data, el)
            if info['length'] != -1 and info['offset'] != -1 and info['offset_head'] != -1:
                for x in self.data_assert:
                    if el in x:
                        sequence = x
                        break
                info_seq = sequence.split('\t')
                self.assertEqual([el, info['length'], info['offset']], [info_seq[0], int(info_seq[1]), int(info_seq[2])])

    def test_get_range_sequence(self):
        list = [{'min': dict['min_range'], 'max': dict['max_range']} for dict in self.data]
        range = list[random.randint(0, len(list))]
        print(f"{range['min']}-{range['max']}: {fp.get_sequences_of_range(self.data, range['min'], range['max'])}")

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
                        info[4] = length
                    info_assert = self.data_assert[j].split('\t')
                    self.assertEqual([info[0], info[4], info[3]], [info_assert[0],info_assert[1],info_assert[2]])

                    last_seq = info[0]
                    j += 1


if __name__ == '__main__':
    unittest.main()
