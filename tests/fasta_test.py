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
tests.fasta_test
~~~~~~~~~~~~~~~~~~~~

'''
import os.path as op
# import sys
import unittest

from libnano.fileio import fasta

# import conftest


LOCAL_DIR: str = op.abspath(op.dirname(__file__))


class TestFastaParsing(unittest.TestCase):

    def setUp(self):
        self.ex1 = op.join(LOCAL_DIR, 'test_data/fasta_in1.fa')
        self.ex2 = op.join(LOCAL_DIR, 'test_data/fasta_in2.fa')

    def _checkParsedRecord(self, rec_id, rec_seq):
        self.assertTrue(
            '\n' not in rec_seq and ' ' not in rec_seq,
        )
        self.assertTrue(
            '>' not in rec_seq,
        )
        self.assertTrue(
            len(rec_id) > 0 and len(rec_seq) > 0,
        )

    def test_parseFasta(self):
        parsed_input = fasta.parseFasta(
            self.ex1,
            alphabet='AMINO_ACID',
        )
        for rec_id, rec_seq in parsed_input:
            self._checkParsedRecord(
                rec_id,
                rec_seq,
            )
        self.assertEqual(
            len(parsed_input),
            3,
        )
        parsed_input = fasta.parseFasta(
            self.ex2,
            alphabet='AMINO_ACID',
        )
        self.assertEqual(
            len(parsed_input),
            3,
        )

    def test_parseFastaGen(self):
        num_recs = 0
        for rec_id, rec_seq in fasta.parseFastaGen(
            self.ex1,
            alphabet='AMINO_ACID',
        ):
            num_recs += 1
            self._checkParsedRecord(
                rec_id,
                rec_seq,
            )
        self.assertEqual(num_recs, 3)
        num_recs = 0
        for rec_id, rec_seq in fasta.parseFastaGen(
            self.ex1,
            alphabet='AMINO_ACID',
        ):
            num_recs += 1
            self._checkParsedRecord(
                rec_id,
                rec_seq,
            )
        self.assertEqual(
            num_recs,
            3,
        )

    def test_sanitizeRecSeq(self):
        for seq_id, req_seq in fasta.parseFastaGen(
            self.ex1,
            alphabet='AMINO_ACID',
        ):
            pass
        with self.assertRaises(ValueError):
            for seq_id, req_seq in fasta.parseFastaGen(
                self.ex1,
                alphabet='DNA',
            ):
                pass
        with self.assertRaises(ValueError):
            for seq_id, req_seq in fasta.parseFastaGen(
                self.ex1,
                not_allowed='N',
            ):
                pass


if __name__ == '__main__':
    unittest.main(verbosity=2)
