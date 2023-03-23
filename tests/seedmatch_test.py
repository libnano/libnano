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
tests.seedmatch_test
~~~~~~~~~~~~~~~~~~~~

'''
from __future__ import print_function

# import random
import unittest

from libnano.search.seedmatch import SeedMatcher  # type: ignore

# import conftest


DEBUG: bool = False


class TestSeedMatch(unittest.TestCase):

    def setUp(self):
        self.test_seed = '#--#--#'  # solves the 10 2 problem
        # test_seed_pair = ['####---', '-#--###'] # solves the 10, 2 as a pair

        self.target = 'GCACGCTGAT'

        library = ['AAAAAAAGCACGCTGATAAAAAAAAAA']
        library_cut_left = ['ACGCTGATAAAAAAAAAA']
        library_cut_right = ['AAAAAAAGCACGCTG']
        library_cut_right_too_much = ['AAAAAAAGCACGCT']
        library_mid = ['AAAAAAAGCAAGATGATAAAAAAAAAA']

        self.library = (
            library_cut_left +
            library +
            library_cut_right +
            library_cut_right_too_much +
            library_mid
        )

    def testMatch(self):
        '''Test on python basestring input
        '''
        self.simpleRunner(self.target, self.library)

    def testBytesMatch(self):
        '''Test on byte strings input
        '''
        bytes_lib = [x.encode('utf-8') for x in self.library]
        self.simpleRunner(self.target.encode('utf-8'), bytes_lib)

    def simpleRunner(self, target, library):
        sm = SeedMatcher(self.test_seed, library)
        outlist = sm.match(target, 2)
        if DEBUG:
            print('target begin', target)
            for i, match in enumerate(outlist):
                lib_idx, offset = match
                print(
                    'match %d. : lib_idx: %d, offset: %d' %
                    (i, lib_idx, offset),
                )
                if offset < 0:
                    print(target)
                    print(b' ' * abs(offset) + library[lib_idx])
                else:
                    print(b' ' * offset + target)
                    print(library[lib_idx])
            print(outlist)
        # Brute force match
        brute_matches = [(-1, 0)] * len(library)
        for i, seq in enumerate(library):
            if target in seq:
                brute_matches[i] = (i, seq.index(target))
        # Only compare the brute force exact match
        self.assertEqual(brute_matches[1], outlist[1])


if __name__ == '__main__':
    unittest.main(verbosity=2)
