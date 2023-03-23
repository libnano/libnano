# Copyright (C) 2014-2018. Nick Conway & Ben Pruitt; Wyss Institute
# Copyright (C) 2023 Nick Conway & Ben Pruitt;
# See LICENSE.TXT for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
tests.seqrepeat_test
~~~~~~~~~~~~~~~~~~~~

'''
# import random
# import time
import unittest

# from libnano.fileio import fasta
from libnano.metric import seqrepeat  # type: ignore

# import conftest


class TestSeqRepeat(unittest.TestCase):

    def test_gap_check(self):
        '''Test the gapCheck
        70 base window no repeats more than 8 base pairs
        must be a 70 bp gap between first base of repeat and
        the start of the next instance
        3'  and 5' GC content high
        no repeat at 5' end
        '''
        key = 'CATTACTCTAAGTAAGCTGGTCTGTCGCCA'  # Repeat 30bp
        # Space 41 bp
        gap1_str = 'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT'
        gap2_str = (
            'GGCCGAAGGGGGTGACCGACAAACGGCGGCGGGATTTCGCGCCAGGCGTAGTGTAACA'
            'TGACCGCTAGGA'
        )
        seq = f'{key}{gap1_str}{key}{gap2_str}'
        a = seqrepeat.RepeatCheck(
            seq,
            8,
        )
        check = a.gapCheck(70)
        self.assertEqual(
            len(check),
            1,
        )
        i, j, size = check[0]
        self.assertEqual(
            size,
            8,
        )
        self.assertEqual(
            seq[i:i + size],
            seq[j:j + size],
        )
        check = a.gapCheck(72)
        self.assertEqual(
            len(check),
            2,
        )
        i, j, size = check[0]
        self.assertEqual(
            size,
            30,
        )
        self.assertEqual(
            seq[i:i + size],
            seq[j:j + size],
        )

    def test_all_repeats(self):
        # Repeat  30bp
        key = 'CATTACTCTAAGTAAGCTGGTCTGTCGCCA'
        gap1_str = 'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT'  # space 41 bp
        gap2_str = (
            'GGCCGAAGGGGGTGACCGACAAACGGCGGCGGGATTTCGCGCCAGGCGTAGTGTAACA'
            'TGACCGCTAGGA'
        )
        seq = f'{key}{gap1_str}{key}{gap2_str}'

        a = seqrepeat.RepeatCheck(
            seq,
            8,
        )
        check = a.allRepeats()
        self.assertIn(
            key,
            check,
        )
        self.assertEqual(
            check[key],
            [0, 71],
        )
        key2 = 'CGCCAGGC'
        self.assertIn(
            key2,
            check,
        )
        self.assertEqual(
            check[key2],
            [96, 140],
        )
        self.assertEqual(
            len(check),
            2,
        )

    def test_indices_of(self):
        '''Test basic construction and indicesOf
        '''
        a = seqrepeat.RepeatCheck(
            'GGGGAAAATTTTTTATTTTTGGGGAAAAT',
            8,
        )
        self.assertEqual(
            a.indicesOf('AAAATTTT'),
            [4],
        )
        self.assertEqual(
            a.indicesOf('GGGGAAAAT'),
            [0, 20],
        )
        self.assertEqual(
            a.indicesOf('ZZZZZZZZZ'),
            [],
        )

        a = seqrepeat.RepeatCheck(
            'GGGGAAAATTTTTTATTTTTGGGGAAAAT',
            6,
        )
        self.assertEqual(
            a.indicesOf('AAAATTTT'),
            [4],
        )
        self.assertEqual(
            a.indicesOf('GGGGAAA'),
            [0, 20],
        )
        self.assertEqual(
            a.indicesOf('ZZZZZZZ'),
            [],
        )

        longest = 'GGGGAAAATTTTTTATTTTTGGGGAAAATGGA'    # 32 bases
        a = seqrepeat.RepeatCheck(
            longest,
            32,
        )
        self.assertEqual(
            a.indicesOf(longest),
            [0],
        )

        with self.assertRaises(ValueError):
            seqrepeat.RepeatCheck('AAAAAAAA', 32)

    def test_screen_ends(self):
        # Repeat 30 bp
        key = 'CATTACTCTAAGTAAGCTGGTCTGTCGCCA'

        gap1_str = 'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT'  # Space 41 bp

        seq = f'{key}{gap1_str}{key}'

        a = seqrepeat.RepeatCheck(seq, 8)
        check = a.screenEnds(8)
        self.assertEqual(check, (True, True))

        # (False, True)
        seq = (
            f'AAACCCCGC{key}'
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT'  # space 41 bp
            f'{key}'
        )
        a = seqrepeat.RepeatCheck(
            seq,
            8,
        )
        check = a.screenEnds(8)
        self.assertEqual(
            check,
            (False, True),
        )

        # (False, False)
        seq = (
            f'AAACCCCGC{key}'  # repeat 1 30bp
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT'  # space 41 bp
            f'{key}GGGAACCTTTA'  # repeat 2 30 bp
        )
        a = seqrepeat.RepeatCheck(
            seq,
            8,
        )
        check = a.screenEnds(8)
        self.assertEqual(
            check,
            (False, False),
        )

        # (True, False)
        seq = (
            f'{key}CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT'  # space 41 bp
            f'{key}GGGAACCTTTA'  # repeat 2 30 bp
        )
        a = seqrepeat.RepeatCheck(
            seq,
            8,
        )
        check = a.screenEnds(8)
        self.assertEqual(
            check,
            (True, False),
        )

    def test_window_check(self):
        a = seqrepeat.RepeatCheck(
            'GGGGAAAATGGGGAAAA',
            8,
        )
        violations, counts = a.windowCheck(
            17,
            0,
        )
        self.assertEqual(
            violations,
            [0],
        )
        self.assertEqual(
            counts,
            [1],
        )
        out = a.getRepeatWindow(
            violations[0],
            17,
        )
        self.assertIn(
            'GGGGAAAA',
            out,
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
