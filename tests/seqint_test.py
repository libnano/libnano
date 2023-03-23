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
tests.seqint_test
~~~~~~~~~~~~~~~~~

'''
import random
import unittest

import conftest

from libnano import seqint  # type: ignore

# ~~~~~ Python implementations of sequence manipulations for comparison ~~~~~ #

_DNACOMP = str.maketrans(
    'ACGTacgt',
    'TGCATGCA',
)


def reverseComplement(seq):
    return seq.translate(_DNACOMP)[::-1]


def complement(seq):
    return seq.translate(_DNACOMP)


def add_to_window(seq, base):
    return seq[1:] + base


def addBase(seq, base):
    return seq + base

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class TestSeqInt(unittest.TestCase):

    def _rand_seq(self):
        ''' Generate random seq of length between 2 and 30 '''
        return ''.join([
            random.choice('ATGCatgc') for x in
            range(random.randint(2, 30))
        ])

    def test_basic_conversions(self):
        ''' Test conversion to seqint and back for various random seqs '''
        for x in range(1000):
            seq = self._rand_seq()
            seq_int = seqint.seq2Int(seq)
            out_seq = seqint.int2Seq(
                seq_int,
                len(seq),
            )
            self.assertEqual(
                out_seq,
                seq.upper(),
            )

    def test_window_addition(self):
        ''' Test use of seqint as circular buffer via base addition '''
        for x in range(1000):
            seq = self._rand_seq()
            base = random.choice('ATGC')
            seq_int = seqint.seq2Int(seq)
            seq_int_add = seqint.addToWindow(
                seq_int,
                base,
                len(seq),
            )
            seq_add = seqint.int2Seq(
                seq_int_add,
                len(seq),
            )
            py_seq_add = add_to_window(
                seq,
                base,
            )
            self.assertEqual(
                seq_add,
                py_seq_add.upper(),
            )

    def test_reverse_complement(self):
        ''' Test reverse complement of seqint '''
        for x in range(1000):
            seq = self._rand_seq()
            seq_int = seqint.seq2Int(seq)
            seq_int_rc = seqint.reverseComplement(
                seq_int,
                len(seq),
            )
            seq_rc = seqint.int2Seq(
                seq_int_rc,
                len(seq),
            )
            py_seq_rc = reverseComplement(seq)
            self.assertEqual(
                seq_rc,
                py_seq_rc.upper(),
                msg=seq,
            )

    def test_substring(self):
        ''' Test seqint substring method '''
        for x in range(1000):
            seq = self._rand_seq()
            seq_int = seqint.seq2Int(seq)
            sidx = random.randint(
                0,
                len(seq) - 1,
            )
            eidx = random.randint(
                sidx + 1,
                len(seq),
            )
            sub_seq_int = seqint.getSubstring(
                seq_int,
                sidx,
                eidx,
                len(seq),
            )
            self.assertEqual(
                seqint.int2Seq(
                    sub_seq_int,
                    eidx - sidx,
                ),
                seq[sidx:eidx].upper(),
            )

    def test_edge_cases(self):
        # Should not allow input sequence of over 30 bases
        self.assertRaises(
            ValueError,
            seqint.seq2Int,
            'A' * 31,
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
