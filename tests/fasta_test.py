# -*- coding: utf-8 -*-
import unittest
import sys
from os.path import (
    join,
    abspath,
    dirname
)

import conftest
from libnano.fileio import fasta

LOCAL_DIR: str = abspath(dirname(__file__))

class TestFastaParsing(unittest.TestCase):

    def setUp(self):
        self.ex1 = join(LOCAL_DIR, 'test_data/fasta_in1.fa')
        self.ex2 = join(LOCAL_DIR, 'test_data/fasta_in2.fa')

    def _checkParsedRecord(self, rec_id, rec_seq):
        self.assertTrue('\n' not in rec_seq and ' ' not in rec_seq)
        self.assertTrue('>' not in rec_seq)
        self.assertTrue(len(rec_id)>0 and len(rec_seq)>0)

    def test_parseFasta(self):
        parsed_input = fasta.parseFasta(self.ex1, alphabet='AMINO_ACID')
        for rec_id, rec_seq in parsed_input:
            self._checkParsedRecord(rec_id, rec_seq)
        self.assertEqual(len(parsed_input), 3)
        parsed_input = fasta.parseFasta(self.ex2, alphabet='AMINO_ACID')
        self.assertEqual(len(parsed_input), 3)

    def test_parseFastaGen(self):
        num_recs = 0
        for rec_id, rec_seq in fasta.parseFastaGen(self.ex1,
                                                   alphabet='AMINO_ACID'):
            num_recs += 1
            self._checkParsedRecord(rec_id, rec_seq)
        self.assertEqual(num_recs, 3)
        num_recs = 0
        for rec_id, rec_seq in fasta.parseFastaGen(self.ex1,
                                                   alphabet='AMINO_ACID'):
            num_recs += 1
            self._checkParsedRecord(rec_id, rec_seq)
        self.assertEqual(num_recs, 3)

    def test_sanitizeRecSeq(self):
        for seq_id, req_seq in fasta.parseFastaGen(self.ex1,
                                                       alphabet='AMINO_ACID'):
            pass
        with self.assertRaises(ValueError):
            for seq_id, req_seq in fasta.parseFastaGen(self.ex1,
                                                       alphabet='DNA'):
                pass
        with self.assertRaises(ValueError):
            for seq_id, req_seq in fasta.parseFastaGen(self.ex1,
                                                       not_allowed='N'):
                pass


if __name__ == '__main__':
    unittest.main(verbosity=2)