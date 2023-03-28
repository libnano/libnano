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
tests.test_seqrepeat
~~~~~~~~~~~~~~~~~~~~

'''
import random
import unittest

from libnano.metric import seqscreen  # type: ignore

# import conftest


class TestSeqFilter(unittest.TestCase):

    def _rand_seq(self):
        '''Generate random seq of length between 2 and 30

        '''
        return ''.join([
            random.choice('ATGCatgc') for x in
            range(random.randint(10, 300))
        ])

    def test_contains_run(self):
        self.assertTrue(
            seqscreen.containsRun(
                'AATGC', 2, 2, 2, 2, 3, 2,
            ),
        )
        self.assertFalse(
            seqscreen.containsRun(
                'AATGC', 1, 2, 2, 2, 3, 2,
            ),
        )
        self.assertFalse(
            seqscreen.containsRun(
                'AATGC', 2, 2, 2, 2, 2, 2,
            ),
        )

    def test_gc_window(self):
        for _ in range(1000):
            seq = self._rand_seq()
            gc_min = random.randint(
                0, 99,
            )
            gc_max = random.randint(
                gc_min, 100,
            )
            window_size = random.randint(
                4, len(seq),
            )
            window_pass, window_start, gc_count = seqscreen.gcWindow(
                seq,
                gc_min,
                gc_max,
                window_size,
            )
            if not window_pass:
                v_gc_count = (
                    seq[window_start:window_start + window_size].count('G') +
                    seq[window_start:window_start + window_size].count('C')
                )
                self.assertEqual(
                    gc_count,
                    v_gc_count,
                )


if __name__ == '__main__':
    unittest.main(verbosity=2)
