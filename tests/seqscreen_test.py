import sys

import random
import unittest
from os.path import join, abspath, dirname

# For package imports
sys.path.insert(0, abspath(join(dirname(__file__), '..')))
import _setup

from libnano.metric import seqscreen

class TestSeqFilter(unittest.TestCase):

    def _randSeq(self):
        ''' Generate random seq of length between 2 and 30 '''
        return ''.join([random.choice('ATGCatgc') for x in
                        range(random.randint(10, 300))])

    def test_containsRun(self):
        self.assertTrue(seqscreen.containsRun('AATGC', 2,2,2,2,3,2))
        self.assertFalse(seqscreen.containsRun('AATGC', 1,2,2,2,3,2))
        self.assertFalse(seqscreen.containsRun('AATGC', 2,2,2,2,2,2))

    def test_gcWindow(self):
        for _ in range(1000):
            seq = self._randSeq()
            gc_min = random.randint(0, 99)
            gc_max = random.randint(gc_min, 100)
            window_size = random.randint(4, len(seq))
            window_pass, window_start, gc_count = seqscreen.gcWindow(
                seq, gc_min, gc_max, window_size)
            if not window_pass:
                v_gc_count = seq[window_start:window_start+window_size].count('G') + \
                             seq[window_start:window_start+window_size].count('C')
                self.assertEqual(gc_count, v_gc_count)


if __name__ == '__main__':
    unittest.main(verbosity=2)