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
tests.test_submerpool
~~~~~~~~~~~~~~~~~~~~~

'''
from __future__ import print_function

import os.path as op
import random
import re
# import sys
import time
import unittest

import pytest

from libnano.fileio import fasta
from libnano.search import submerpool  # type: ignore

LOCAL_DIR = op.dirname(op.abspath(__file__))

'''TODO: Fix this file for the submer module

'''


@pytest.mark.skip(reason='Submer testing needs to be updated / refactored')
class TestSubmer(unittest.TestCase):

    def setUp(self):
        recs = fasta.parse_fasta(
            op.join(
                LOCAL_DIR,
                'test_data/MG1655.fa',
            ),
        )
        self.mg1655 = recs[0][1]

    def find_matches_naive(
            self,
            target,
            search_seq,
            k,
    ):
        '''
        Naive match finding by generating list of regular expressions
        returns list

        '''
        wc = ['[ATGC]']
        target = list(target)
        regex_list = []
        target_len = len(target)

        def _inner(cur_regex, idx_prev, k_remain):
            if k_remain == 0:
                regex_list.append(''.join(cur_regex))
            for idx in range(idx_prev + 1, target_len):
                _inner(
                    cur_regex[:idx] + wc +
                    cur_regex[idx + 1:], idx, k_remain - 1,
                )
        _inner(target, -1, k)
        match_idxs = []
        for r in regex_list:
            match_idxs += [
                m.start()
                for m in re.finditer(r, search_seq)
            ]
        # return match_idxs
        if match_idxs:
            match_idxs = list(set(match_idxs))
            match_idxs.sort()
        return match_idxs

    def test_create_submer_search(self):
        # Test mismatch creation
        for x in range(10):
            m = random.randint(5, 20)
            k = random.randint(1, 3)
            # sd_list = seedfinder.findSeed(m, k)
            tbl = submerpool.SubmerPoolSearch(
                [self.mg1655],
                m,
                mismatches=k,
            )
            self.assertNotEqual(tbl, None)
        # Test no mismatchs
        tbl = submerpool.SubmerPoolSearch(
            [self.mg1655],
            m,
            mismatches=0,
        )
        self.assertNotEqual(tbl, None)
        # Test no mismatchs
        tbl = submerpool.SubmerPoolSearch(
            [self.mg1655],
            m,
            force_hamming=True,
        )
        self.assertEqual(tbl.get_matchers(), None)

    def test_find_matches(self):
        m = 10  # Initialize to ensure a seed check
        k = 1
        for i in range(2):
            print('* Round %d: %d, %d' % (i, m, k))
            target = ''.join(
                [
                    random.choice('ATGC')
                    for j in range(m)
                ],
            )
            tbl = submerpool.SubmerPoolSearch(
                [self.mg1655],
                m,
                mismatches=k,
            )
            mt_start = time.time()
            mt = tbl.find(target)
            if len(mt) == 0:
                continue
            print('    reg:', len(mt))
            mt_end = time.time()
            print('    mt_time: %0.2f' % (mt_end - mt_start))
            mt_naive = self.find_matches_naive(target, self.mg1655, k)
            mt_naive_end = time.time()
            print('    naive', len(mt_naive))
            print('    mt_naive_time: %0.2f' % (mt_naive_end - mt_end))
            mt_offsets = []
            i = 0
            for ref_idx, offset in mt:
                # print(self.mg1655[item:item+len(target)])
                # self.assertIn(offset, mt_naive)
                mt_offsets.append(offset)
                if offset not in mt_naive:
                    print('    REG', i)
                    i += 1
                    print(target)
                    print(self.mg1655[offset:offset + len(target)])
            for offset in mt_naive:
                # print(self.mg1655[item:item+len(target)])
                # self.assertIn(offset, mt_naive)
                if offset not in mt_offsets:
                    print('    NAIVE', i)
                    i += 1
                    print(target)
                    print(self.mg1655[offset:offset + len(target)])

            print('    target: ' + target + ', k:', k)
            # Apparently the edge case is missing something
            # so >= is in order, look for:
            # AAAACTGGGA
            # AAAACTGGAA
            # in mg1655, that is naive doesn't find it
            self.assertGreaterEqual(
                len(mt),
                len(mt_naive),
            )
            # if len(mt) > 0:
            #     break
            m = random.randint(10, 14)
            k = random.randint(1, 2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
