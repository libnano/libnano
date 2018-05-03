# -*- coding: utf-8 -*-
import timeit
import random
import unittest

import _setup
from libnano import seqstr

# ~~~~~ Python implementations of sequence manipulations for comparison ~~~~~ #

_compLUT = str.maketrans(
    'ACGTUMRWSYKVHDBN\nacgtumrwsykvhdbn',
    'TGCAAKYWSRMBDHVN\ntgcaakywsrmbdhvn'
)

def reverse(seq):
    return seq[::-1]

def reverseComplement(seq):
    return seq.translate(_compLUT)[::-1]

def complement(seq):
    return seq.translate(_compLUT)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

class TestSeqStr(unittest.TestCase):

    def _randSeq(self, min_len=2, max_len=100):
        ''' Generate random seq of length between 2 and 1000 '''
        return ''.join([random.choice('ACGTUMRWSYKVHDBN\nacgtumrwsykvhdbn')
                        for x in range(random.randint(min_len, max_len))])

    def test_rev(self):
        for x in range(1000):
            seq = self._randSeq()
            py_rev = reverse(seq)
            rev = seqstr.reverse(seq)
            self.assertEqual(py_rev, rev)

    def test_comp(self):
        for x in range(1000):
            seq = self._randSeq()
            py_comp = complement(seq)
            comp = seqstr.complement(seq)
            self.assertEqual(py_comp, comp)

    def test_revComp(self):
        for x in range(1000):
            seq = self._randSeq()
            py_revcomp = reverseComplement(seq)
            revcomp = seqstr.reverseComplement(seq)
            self.assertEqual(py_revcomp, revcomp)

    def test_revCompProfile(self):
        py_time = timeit.timeit(lambda: reverseComplement("ACGTUMRWSYKVHDBNACGTUMRWSYKVHDBN"),
                                number=10000)
        capi_time = timeit.timeit(lambda: seqstr.reverseComplement(
                                  "ACGTUMRWSYKVHDBNACGTUMRWSYKVHDBN"), number=10000)
        print('\nNative python rev_comp (time for 10000X):', py_time)
        print('C API python rev_comp (time for 10000X): ', capi_time)

    def test_hamming(self):
        for x in range(1000):
            s1 = self._randSeq(20, 100)
            num_mismatches = random.randint(0, len(s1))
            shuffled_bases = list(range(0, len(s1)))
            random.shuffle(shuffled_bases)
            s2 = list(s1)
            for x in range(num_mismatches):
                random_idx = shuffled_bases[x]
                s2[random_idx] = random.choice('ATGC'.replace(s1[random_idx], ''))
            s2 = ''.join(s2)
            hd = seqstr.hammingDistance(s1, s2)
            self.assertEqual(num_mismatches, hd, msg='\ns1: {}\ns2: {}'.format(s1, s2))

    # def test_can3pMisPrime(self):
    #         num_bases = 15
    #         thresholds = [0 for x in range(num_bases)]
    #         thresholds[0] = -1
    #         thresholds[1] = -1
    #         thresholds[2] = 2
    #         thresholds[3] = 2
    #         thresholds[4] = 2
    #         thresholds[5] = 1
    #         thresholds[6] = 1

    #         genome =        "ACGTTTACGTAATTTAGTTTGC"

    #         primer_pass =   "TCTTTAAATGGTGGC"
    #         self.assertTrue(seqstr.can3pMisprime(primer_pass,
    #                         genome,
    #                         thresholds,
    #                         False))

    #         primer_fail =   "CCGTAATTTAGTTTG"
    #         self.assertFalse(seqstr.can3pMisprime(primer_fail,
    #             genome,
    #             thresholds,
    #             False))


    #         primer_check_must_3p_mismatch_false =     "TCTTTAAATGTTTGC"
    #         self.assertFalse(seqstr.can3pMisprime(
    #             primer_check_must_3p_mismatch_false,
    #             genome,
    #             thresholds,
    #             True))

    #         primer_check_must_3p_mismatch_true =     "TCCGATGCTGGTCCT"
    #         self.assertTrue(seqstr.can3pMisprime(
    #             primer_check_must_3p_mismatch_true,
    #             genome,
    #             thresholds,
    #             True))

    # # end def


if __name__ == '__main__':
    unittest.main(verbosity=2)

