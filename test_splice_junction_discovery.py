import unittest
import numpy as np
from SpliceJunctionDiscovery import num_of_matches_before_first_intronic_section, length_of_first_intronic_section,\
    get_first_splice_junction, get_end_position_from_junction, get_second_splice_junction


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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def test_get_first_splice_junction_s(self):
        # Read positions are inclusive at start, not inclusive at end
        cigar_string = '8S13M221400N34M2658N29M'
        pos = 5
        t_chrom = 2
        sj = get_first_splice_junction(cigar_string, 2, pos)
        self.assertEqual(sj, '{},{},{}'.format(t_chrom, pos+13, pos+13+221400))

    def test_get_first_splice_junction_double(self):
        # Read positions are inclusive at start, not inclusive at end
        cigar_string = '13M221400N34M2658N29M'
        pos = 5
        t_chrom = 2
        sj = get_first_splice_junction(cigar_string, 2, pos)
        self.assertEqual(sj, '{},{},{}'.format(t_chrom, pos+13, pos+13+221400))

    def test_get_first_splice_junction_deletion(self):
        # Read positions are inclusive at start, not inclusive at end
        cigar_string = '3M1D40M20N'
        pos = 5
        t_chrom = 2
        sj = get_first_splice_junction(cigar_string, 2, pos)
        self.assertEqual(sj, '{},{},{}'.format(t_chrom, pos+3+40+1, pos+3+40+1+20))

    def test_get_first_splice_junction_insertion(self):
        # Read positions are inclusive at start, not inclusive at end
        cigar_string = '22M1I19M1893N6M'
        pos = 5
        t_chrom = 2
        sj = get_first_splice_junction(cigar_string, 2, pos)
        self.assertEqual(sj, '{},{},{}'.format(t_chrom, pos+22+19, pos+22+19+1893))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def test_get_end_position_from_junction(self):
        self.assertEqual(get_end_position_from_junction('1,3,5'), 5)
        self.assertEqual(get_end_position_from_junction('1,3,35'), 35)
        self.assertEqual(get_end_position_from_junction('1,3,20005'), 20005)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def test_get_second_splice_junction_double(self):
        # Read positions are inclusive at start, not inclusive at end
        cigar_string = '13M221400N34M2658N29M'
        pos = 5
        t_chrom = 2

        sj = get_first_splice_junction(cigar_string, 2, pos)
        intron_end_position = get_end_position_from_junction(sj)  # not inclusive
        self.assertEqual(intron_end_position, pos+13+221400)

        sj = get_second_splice_junction(cigar_string, 2, pos, intron_end_position)
        self.assertEqual(sj, '{},{},{}'.format(t_chrom, intron_end_position+34, intron_end_position+34+2658))


if __name__ == '__main__':
    unittest.main()
