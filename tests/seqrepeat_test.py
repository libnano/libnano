# -*- coding: utf-8 -*-
import random
import unittest
import time

import _setup
from libnano.metric import seqrepeat
from libnano.fileio import fasta


class TestSeqRepeat(unittest.TestCase):

    def test_gapCheck(self):
        """ test the gapCheck
        70 base window no repeats more than 8 base pairs
        must be a 70 bp gap between first base of repeat and
        the start of the next instance
        3'  and 5' GC content high
        no repeat at 5' end
        """
        key = 'CATTACTCTAAGTAAGCTGGTCTGTCGCCA'
        seq = ( key +   # repeat 1 30bp
                'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
                key + # repeat 2 30 bp
                "GGCCGAAGGGGGTGACCGACAAACGGCGGCGGGATTTCGCGCCAGGCGTAGTGTAACATGACCGCTAGGA"
            )
        a = seqrepeat.RepeatCheck(seq, 8)
        check = a.gapCheck(70)
        self.assertEqual(len(check), 1)
        i, j, size = check[0]
        self.assertEqual(size, 8)
        self.assertEqual(seq[i:i + size], seq[j:j + size])
        check = a.gapCheck(72)
        self.assertEqual(len(check), 2)
        i, j, size = check[0]
        self.assertEqual(size, 30)
        self.assertEqual(seq[i:i + size], seq[j:j + size])

    def test_allRepeats(self):
        key = 'CATTACTCTAAGTAAGCTGGTCTGTCGCCA'
        seq = ( key +   # repeat 1 30bp
                'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
                key + # repeat 2 30 bp
                "GGCCGAAGGGGGTGACCGACAAACGGCGGCGGGATTTCGCGCCAGGCGTAGTGTAACATGACCGCTAGGA"
            )
        a = seqrepeat.RepeatCheck(seq, 8)
        check = a.allRepeats()
        self.assertIn(key, check)
        self.assertEqual(check[key], [0, 71])
        key2 = 'CGCCAGGC'
        self.assertIn(key2, check)
        self.assertEqual(check[key2], [96, 140])
        self.assertEqual(len(check), 2)

    def test_indicesOf(self):
        """ test basic construction and indicesOf
        """
        a = seqrepeat.RepeatCheck('GGGGAAAATTTTTTATTTTTGGGGAAAAT', 8)
        self.assertEqual(a.indicesOf('AAAATTTT'), [4])
        self.assertEqual(a.indicesOf('GGGGAAAAT'), [0, 20])
        self.assertEqual(a.indicesOf('ZZZZZZZZZ'), [])

        a = seqrepeat.RepeatCheck('GGGGAAAATTTTTTATTTTTGGGGAAAAT', 6)
        self.assertEqual(a.indicesOf('AAAATTTT'), [4])
        self.assertEqual(a.indicesOf('GGGGAAA'), [0, 20])
        self.assertEqual(a.indicesOf('ZZZZZZZ'), [])

        longest = 'GGGGAAAATTTTTTATTTTTGGGGAAAATGGA'    # 32 bases
        a = seqrepeat.RepeatCheck(longest, 32)
        self.assertEqual(a.indicesOf(longest), [0])

        with self.assertRaises(ValueError):
            seqrepeat.RepeatCheck("AAAAAAAA", 32)

    def test_screenEnds(self):
        key = 'CATTACTCTAAGTAAGCTGGTCTGTCGCCA'
        # (True, True)
        seq = ( key +   # repeat 1 30bp
                'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
                key # repeat 2 30 bp
            )
        a = seqrepeat.RepeatCheck(seq, 8)
        check = a.screenEnds(8)
        self.assertEqual(check, (True, True))

        # (False, True)
        seq = ( "AAACCCCGC" + key +   # repeat 1 30bp
                'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
                key # repeat 2 30 bp
            )
        a = seqrepeat.RepeatCheck(seq, 8)
        check = a.screenEnds(8)
        self.assertEqual(check, (False, True))

        # (False, False)
        seq = ( "AAACCCCGC" + key +   # repeat 1 30bp
                'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
                key + "GGGAACCTTTA" # repeat 2 30 bp
            )
        a = seqrepeat.RepeatCheck(seq, 8)
        check = a.screenEnds(8)
        self.assertEqual(check, (False, False))

        # (True, False)
        seq = ( key +   # repeat 1 30bp
                'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
                key + "GGGAACCTTTA" # repeat 2 30 bp
            )
        a = seqrepeat.RepeatCheck(seq, 8)
        check = a.screenEnds(8)
        self.assertEqual(check, (True, False))

    def test_windowCheck(self):
        a = seqrepeat.RepeatCheck('GGGGAAAATGGGGAAAA', 8)
        violations, counts = a.windowCheck(17, 0)
        self.assertEqual(violations, [0])
        self.assertEqual(counts, [1])
        out = a.getRepeatWindow(violations[0], 17)
        self.assertIn('GGGGAAAA', out)

if __name__ == '__main__':
    unittest.main(verbosity=2)