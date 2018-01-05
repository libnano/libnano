from __future__ import print_function

# For package imports
import sys
from os.path import join, abspath, dirname
sys.path.insert(0, abspath(join(dirname(__file__), '..')))

from libnano.core.seqsearch.seedmatch import SeedMatcher
from libnano.core.seqsearch import seedfinder

import random
import unittest
import _setup

DEBUG = False

class TestSeedMatch(unittest.TestCase):

    def setUp(self):
        self.test_seed = '#--#--#' # solves the 10 2 problem
        # test_seed_pair = ['####---', '-#--###'] # solves the 10, 2 as a pair

        self.target = "GCACGCTGAT"

        library = ["AAAAAAAGCACGCTGATAAAAAAAAAA"]
        library_cut_left = ["ACGCTGATAAAAAAAAAA"]
        library_cut_right = ["AAAAAAAGCACGCTG"]
        library_cut_right_too_much = ["AAAAAAAGCACGCT"]
        library_mid = ["AAAAAAAGCAAGATGATAAAAAAAAAA"]

        self.library = library_cut_left + library + library_cut_right + \
                    library_cut_right_too_much + library_mid

    # end def

    def testMatch(self):
        """ test on python basestring input
        """
        self.simpleRunner(self.target, self.library)
    # end def

    def testBytesMatch(self):
        """ test on byte strings input
        """
        bytes_lib = [x.encode('utf-8') for x in self.library]
        self.simpleRunner(self.target.encode('utf-8'), bytes_lib)
    # end def

    def simpleRunner(self, target, library):
        sm  = SeedMatcher(self.test_seed, library)
        outlist = sm.match(target, 2)
        if DEBUG:
            print("target begin", target)
            for i, match in enumerate(outlist):
                lib_idx, offset = match
                print("match %d. : lib_idx: %d, offset: %d" % (i, lib_idx, offset))
                if offset < 0:
                    print(target)
                    print(b" "*abs(offset) + library[lib_idx])
                else:
                    print(b" "*offset+target)
                    print(library[lib_idx])
            print(outlist)
        # brute force match
        brute_matches = [(-1, 0)]*len(library)
        for i, seq in enumerate(library):
            if target in seq:
                brute_matches[i] = (i, seq.index(target))
        # Only compare the brute force exact match
        self.assertEqual(brute_matches[1], outlist[1])
    # end def


if __name__ == '__main__':
    unittest.main(verbosity=2)

