from __future__ import print_function
import sys

import random
import unittest
import re
import time
from os.path import join, abspath, dirname

# For package imports
sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import _setup

from libnano.search import submerpool, seedfinder
from libnano.fileio import fasta

LOCAL_DIR = abspath(dirname(__file__))

""" TODO Fix this file for the submer module

This also tests
"""

@unittest.skip("Submer testing needs to be updated / refactored")
class TestSubmer(unittest.TestCase):

    def setUp(self):
        recs = fasta.parseFasta(join(LOCAL_DIR, 'test_data/MG1655.fa'))
        self.mg1655 = recs[0][1]

    def findMatchesNaive(self, target, search_seq, k):
        """
        Naive match finding by generating list of regular expressions
        returns list of
        """
        wc = ['[ATGC]']
        target = list(target)
        regex_list = []
        target_len = len(target)
        def _inner(cur_regex, idx_prev, k_remain):
            if k_remain == 0:
                regex_list.append(''.join(cur_regex))
            for idx in range(idx_prev+1, target_len):
                _inner(cur_regex[:idx] + wc + cur_regex[idx+1:], idx, k_remain-1)
        _inner(target, -1, k)
        match_idxs = []
        for r in regex_list:
            match_idxs += [m.start() for m in re.finditer(r, search_seq)]
        # return match_idxs
        if match_idxs:
            match_idxs = list(set(match_idxs))
            match_idxs.sort()
        return match_idxs

    def test_createSubmerSearch(self):
        # test mismatch creation
        for x in range(10):
            m = random.randint(5, 20)
            k = random.randint(1,3)
            sd_list = seedfinder.findSeed(m, k)
            tbl = submerpool.SubmerPoolSearch([self.mg1655], m, mismatches=k)
            self.assertNotEqual(tbl, None)
        # test no mismatchs
        tbl = submerpool.SubmerPoolSearch([self.mg1655], m, mismatches=0)
        self.assertNotEqual(tbl, None)
        # test no mismatchs
        tbl = submerpool.SubmerPoolSearch([self.mg1655], m, force_hamming=True)
        self.assertEqual(tbl.getMatchers(), None)
    # end def

    def test_findMatches(self):
        print('')
        m = 10  # initialize to ensure a seed check
        k = 1
        for i in range(2):
            print("* Round %d: %d, %d" % (i, m, k))
            target = ''.join([random.choice('ATGC') for j in range(m)])
            tbl = submerpool.SubmerPoolSearch([self.mg1655], m, mismatches=k)
            mt_start = time.time()
            mt = tbl.find(target)
            if len(mt) == 0:
                continue
            print("    reg:", len(mt))
            mt_end = time.time()
            print('    mt_time: %0.2f' % (mt_end - mt_start))
            mt_naive = self.findMatchesNaive(target, self.mg1655, k)
            mt_naive_end = time.time()
            print("    naive", len(mt_naive))
            print('    mt_naive_time: %0.2f' % (mt_naive_end-mt_end))
            mt_offsets = []
            i = 0
            for ref_idx, offset in mt:
                # print(self.mg1655[item:item+len(target)])
                # self.assertIn(offset, mt_naive)
                mt_offsets.append(offset)
                if offset not in mt_naive:
                    print("    REG", i)
                    i += 1
                    print(target)
                    print(self.mg1655[offset:offset+len(target)])
            for offset in mt_naive:
                # print(self.mg1655[item:item+len(target)])
                # self.assertIn(offset, mt_naive)
                if offset not in mt_offsets:
                    print("    NAIVE", i)
                    i += 1
                    print(target)
                    print(self.mg1655[offset:offset+len(target)])

            print('    target: ' + target + ', k:', k)
            # apparently the edge case is missing something
            # so >= is in order, look for:
            # AAAACTGGGA
            # AAAACTGGAA
            # in mg1655, that is naive doesn't find it
            self.assertGreaterEqual(len(mt), len(mt_naive))
            # if len(mt) > 0:
            #     break
            m = random.randint(10, 14)
            k = random.randint(1,2)
    # end def

if __name__ == '__main__':
    unittest.main(verbosity=2)

