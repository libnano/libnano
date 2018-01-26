import sys

import random
import unittest
from os.path import join, abspath, dirname

# For package imports
sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import _setup

from libnano import seqint

# ~~~~~ Python implementations of sequence manipulations for comparison ~~~~~ #

_PY3 = sys.version_info[0] == 3

if _PY3:
    maketrans = str.maketrans
else:
    from string import maketrans

_DNAcomp = maketrans('ACGTacgt','TGCATGCA')

def reverseComplement(seq):
    return seq.translate(_DNAcomp)[::-1]

def complement(seq):
    return seq.translate(_DNAcomp)

def addToWindow(seq, base):
    return seq[1:] + base

def addBase(seq, base):
    return seq + base

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class TestSeqInt(unittest.TestCase):

    def _randSeq(self):
        ''' Generate random seq of length between 2 and 30 '''
        return ''.join([random.choice('ATGCatgc') for x in
                        range(random.randint(2, 30))])

    def test_basicConversions(self):
        ''' Test conversion to seqint and back for various random seqs '''
        for x in range(1000):
            seq = self._randSeq()
            seq_int = seqint.seq2Int(seq)
            out_seq = seqint.int2Seq(seq_int, len(seq))
            self.assertEqual(out_seq, seq.upper())

    def test_windowAddition(self):
        ''' Test use of seqint as circular buffer via base addition '''
        for x in range(1000):
            seq = self._randSeq()
            base = random.choice('ATGC')
            seq_int = seqint.seq2Int(seq)
            seq_int_add = seqint.addToWindow(seq_int, base, len(seq))
            seq_add = seqint.int2Seq(seq_int_add, len(seq))
            py_seq_add = addToWindow(seq, base)
            self.assertEqual(seq_add, py_seq_add.upper())

    def test_reverseComplement(self):
        ''' Test reverse complement of seqint '''
        for x in range(1000):
            seq = self._randSeq()
            seq_int = seqint.seq2Int(seq)
            seq_int_rc = seqint.reverseComplement(seq_int, len(seq))
            seq_rc = seqint.int2Seq(seq_int_rc, len(seq))
            py_seq_rc = reverseComplement(seq)
            self.assertEqual(seq_rc, py_seq_rc.upper(), msg=seq)

    def test_substring(self):
        ''' Test seqint substring method '''
        for x in range(1000):
            seq = self._randSeq()
            seq_int = seqint.seq2Int(seq)
            sidx = random.randint(0, len(seq)-1)
            eidx = random.randint(sidx+1, len(seq))
            sub_seq_int = seqint.getSubstring(seq_int, sidx, eidx, len(seq))
            self.assertEqual(seqint.int2Seq(sub_seq_int, eidx-sidx), seq[sidx:eidx].upper())

    def test_edgeCases(self):
        # Should not allow input sequence of over 30 bases
        self.assertRaises(ValueError, seqint.seq2Int, 'A' * 31)


if __name__ == '__main__':
    unittest.main(verbosity=2)
