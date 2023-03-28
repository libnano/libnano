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
tests.test_restriction
~~~~~~~~~~~~~~~~~~~~~~

'''
import unittest

from libnano.datasets import (
    dataset_container,
    qc_enzyme_dataset,
)
from libnano.search import restriction  # type: ignore


class TestRestriction(unittest.TestCase):

    def test_enzyme_dataset(self):
        ''' Test the integrity of the current dataset '''
        qc_enzyme_dataset(
            dataset_container.ENZYME_DATASET,
        )

    def test_bytes_support(self):
        rs = restriction.RestrictionSearcher(
            'EcoRI',
            'BsaI',
            'BaeI',
        )
        self.assertRaises(
            TypeError,
            rs.find_sites,
            (b'ATTTGGACCAG'),
        )
        rs = restriction.RestrictionSearcher(
            b'EcoRI',
            b'BsaI',
            b'BaeI',
        )
        self.assertRaises(
            TypeError,
            rs.find_sites,
            ('ATTTGGACCAG'),
        )

    def test_sites_present(self):
        ''' Test RestrictionSearcher.sites_present '''
        rs = restriction.RestrictionSearcher(
            'EcoRI',
            'BsaI',
            'BaeI',
        )
        # Contains the full EcoRI site
        ecori_seq = 'GAATTC'
        self.assertTrue(
            rs.sites_present(ecori_seq),
        )
        self.assertTrue(
            rs.sites_present(
                ecori_seq,
                full_sites=False,
            ),
        )
        # Contains only the core BsaI site
        bsai_core_only_seq = 'GGTCTC'
        self.assertFalse(
            rs.sites_present(bsai_core_only_seq),
        )
        self.assertTrue(
            rs.sites_present(
                bsai_core_only_seq,
                full_sites=False,
            ),
        )
        # Contains both EcoRI and BsaI full sites, overlapping
        ecori_bsai_seq = 'GGTCTCGAATTC'
        self.assertTrue(
            rs.sites_present(ecori_bsai_seq),
        )
        self.assertTrue(
            rs.sites_present(
                ecori_bsai_seq,
                full_sites=False,
            ),
        )
        # Compare instance method to function
        self.assertEqual(
            rs.sites_present(ecori_bsai_seq),
            restriction.sites_present(
                ecori_bsai_seq,
                ['EcoRI', 'BsaI', 'BaeI'],
            ),
        )

    def test_count_sites(self):
        ''' Test RestrictionSearcher.count_sites '''
        rs = restriction.RestrictionSearcher(
            'EcoRI',
            'BsaI',
            'BaeI',
        )
        # Contains one copy of the EcoRI site (palindromic, count should be 2)
        ecori_seq = 'GAATTC'
        counts = rs.count_sites(ecori_seq)
        self.assertTrue(
            counts[0] == 2 and
            counts[1] == 0 and
            counts[2] == 0,
        )
        # Contains only the core BsaI site (not palindromic, count should be 1)
        bsai_core_only_seq = 'GGTCTC'
        # Looking for full sites so all counts should be 0
        counts = rs.count_sites(bsai_core_only_seq)
        self.assertTrue(
            counts[0] == 0 and
            counts[1] == 0 and
            counts[2] == 0,
        )
        # Looking for core sites so BsaI count should be 1
        counts = rs.count_sites(
            bsai_core_only_seq,
            full_sites=False,
        )
        self.assertTrue(
            counts[0] == 0 and
            counts[1] == 1 and
            counts[2] == 0,
        )
        # Contains both EcoRI and BsaI full sites, overlapping
        ecori_bsai_seq = 'GGTCTCGAATTC'
        counts = rs.count_sites(
            ecori_bsai_seq,
        )
        self.assertTrue(
            counts[0] == 2 and
            counts[1] == 1 and
            counts[2] == 0,
        )
        counts = rs.count_sites(
            ecori_bsai_seq,
            full_sites=False,
        )
        self.assertTrue(
            counts[0] == 2 and
            counts[1] == 1 and
            counts[2] == 0,
        )
        # Compare instance method to function
        self.assertEqual(
            rs.count_sites(ecori_bsai_seq),
            restriction.count_sites(
                ecori_bsai_seq,
                ['EcoRI', 'BsaI', 'BaeI'],
            ),
        )

    def test_find_sites(self):
        '''Test RestrictionSearcher.find_sites

        '''
        rs = restriction.RestrictionSearcher(
            'EcoRI',
            'BsaI',
            'BaeI',
        )
        # Contains both EcoRI and BsaI full sites, overlapping
        ecori_bsai_seq = 'GGTCTCGAATTC'
        sites = rs.find_sites(ecori_bsai_seq)
        self.assertTrue(
            len(sites[0]) == 2 and
            len(sites[1]) == 1 and
            len(sites[2]) == 0,
        )
        # EcoRI on fwd strand
        self.assertTrue(
            sites[0][0] == (
                1, 6, 12,
                'EcoRI', 0,
                ((1,), (5,)),
                'GAATTC',
            ),
        )
        # EcoRI on rev strand (indexed from 5' of fwd strand)
        self.assertTrue(
            sites[0][1] == (
                -1, 6, 12,
                'EcoRI', 0,
                ((1,), (5,)),
                'GAATTC',
            ),
        )
        # BsaI on fwd strand
        self.assertTrue(
            sites[1][0] == (
                1, 0, 11,
                'BsaI', 1,
                ((7,), (11,)),
                '[ATGC]{5}GAGACC',
            ),
        )
        # Compare instance method to function
        self.assertEqual(
            rs.find_sites(ecori_bsai_seq),
            restriction.find_sites(
                ecori_bsai_seq,
                ['EcoRI', 'BsaI', 'BaeI'],
            ),
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
