import unittest
import numpy as np
from SpliceJunctionDiscovery import num_of_matches_before_first_intronic_section, length_of_first_intronic_section


class TestDeconstructSigs(unittest.TestCase):
    def test_num_of_matches_before_first_intronic_section_s(self):
        cigar_string = '8S13M221400N34M2658N29M'
        self.assertEqual(num_of_matches_before_first_intronic_section(cigar_string), 13)

    def test_num_of_matches_before_first_intronic_section_double(self):
        cigar_string = '13M221400N34M2658N29M'
        self.assertEqual(num_of_matches_before_first_intronic_section(cigar_string), 13)

    def test_num_of_matches_before_first_intronic_section_single(self):
        cigar_string = '33M20N8M'
        self.assertEqual(num_of_matches_before_first_intronic_section(cigar_string), 33)

    def test_num_of_matches_before_first_intronic_section_deletion(self):
        cigar_string = '3M1D40M20N'
        self.assertEqual(num_of_matches_before_first_intronic_section(cigar_string), 3+1+40)

    def test_num_of_matches_before_first_intronic_section_insertion(self):
        cigar_string = '22M1I19M1893N6M'
        self.assertEqual(num_of_matches_before_first_intronic_section(cigar_string), 22+19)

    def test_length_of_first_intronic_section_double(self):
        cigar_string = '13M221400N34M2658N29M'
        self.assertEqual(length_of_first_intronic_section(cigar_string), 221400)

    def test_length_of_first_intronic_section_single(self):
        cigar_string = '33M20N8M'
        self.assertEqual(length_of_first_intronic_section(cigar_string), 20)

    def test_length_of_first_intronic_section_deletion(self):
            cigar_string = '3M1D40M20N'
            self.assertEqual(length_of_first_intronic_section(cigar_string), 20)

    def test_length_of_first_intronic_section_insertion(self):
        cigar_string = '22M1I19M1893N6M'
        self.assertEqual(length_of_first_intronic_section(cigar_string), 1893)


if __name__ == '__main__':
    unittest.main()
