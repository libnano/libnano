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
tests.test_seqstr
~~~~~~~~~~~~~~~~~

'''
import random
import timeit
import unittest

import conftest

from libnano import seqstr

# ~~~~~ Python implementations of sequence manipulations for comparison ~~~~~ #

_COMPLUT = str.maketrans(
    'ACGTUMRWSYKVHDBN\nacgtumrwsykvhdbn',
    'TGCAAKYWSRMBDHVN\ntgcaakywsrmbdhvn',
)


def reverse(seq):
    return seq[::-1]


def reverseComplement(seq):
    return seq.translate(_COMPLUT)[::-1]


def complement(seq):
    return seq.translate(_COMPLUT)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


class TestSeqStr(unittest.TestCase):

    def _rand_seq(
            self,
            min_len=2,
            max_len=100,
    ):
        ''' Generate random seq of length between 2 and 1000 '''
        return ''.join([
            random.choice('ACGTUMRWSYKVHDBN\nacgtumrwsykvhdbn')
            for x in range(random.randint(min_len, max_len))
        ])

    def test_rev(self):
        for x in range(1000):
            seq = self._rand_seq()
            py_rev = reverse(seq)
            rev = seqstr.reverse(seq)
            self.assertEqual(py_rev, rev)

    def test_comp(self):
        for x in range(1000):
            seq = self._rand_seq()
            py_comp = complement(seq)
            comp = seqstr.complement(seq)
            self.assertEqual(
                py_comp,
                comp,
            )

    def test_rev_comp(self):
        for x in range(1000):
            seq = self._rand_seq()
            py_revcomp = reverseComplement(seq)
            revcomp = seqstr.reverseComplement(seq)
            self.assertEqual(
                py_revcomp,
                revcomp,
            )

    def test_rev_compProfile(self):
        py_time = timeit.timeit(
            lambda: reverseComplement('ACGTUMRWSYKVHDBNACGTUMRWSYKVHDBN'),
            number=10000,
        )
        capi_time = timeit.timeit(
            lambda: seqstr.reverseComplement(
                'ACGTUMRWSYKVHDBNACGTUMRWSYKVHDBN',
            ),
            number=10000,
        )
        print('\nNative python rev_comp (time for 10000X):', py_time)
        print('C API python rev_comp (time for 10000X): ', capi_time)

    def test_hamming(self):
        for x in range(1000):
            s1 = self._rand_seq(20, 100)
            num_mismatches = random.randint(0, len(s1))
            shuffled_bases = list(range(0, len(s1)))
            random.shuffle(shuffled_bases)
            s2 = list(s1)
            for x in range(num_mismatches):
                random_idx = shuffled_bases[x]
                s2[random_idx] = random.choice(
                    'ATGC'.replace(s1[random_idx], ''),
                )
            s2 = ''.join(s2)
            hd = seqstr.hammingDistance(
                s1,
                s2,
            )
            self.assertEqual(
                num_mismatches, hd,
                msg=f'\ns1: {s1}\ns2: {s2}',
            )


if __name__ == '__main__':
    unittest.main(verbosity=2)
