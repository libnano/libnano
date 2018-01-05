
import sys
import unittest

from os.path import join, abspath, dirname

# For package imports
sys.path.insert(0, abspath(join(dirname(__file__), '..')))

from libnano.core.seqsearch import restriction
from libnano.datasets import dataset_container, qcEnzymeDataset
import _setup


class TestRestriction(unittest.TestCase):

    def test_enzymeDataset(self):
        ''' Test the integrity of the current dataset '''
        qcEnzymeDataset(dataset_container.ENZYME_DATASET)

    def test_bytesSupport(self):
        if sys.version_info[0] > 2:
            rs = restriction.RestrictionSearcher('EcoRI', 'BsaI', 'BaeI')
            self.assertRaises(TypeError, rs.findSites, (b'ATTTGGACCAG'))
            rs = restriction.RestrictionSearcher(b'EcoRI', b'BsaI', b'BaeI')
            self.assertRaises(TypeError, rs.findSites, ('ATTTGGACCAG'))

    def test_sitesPresent(self):
        ''' Test RestrictionSearcher.sitesPresent '''
        rs = restriction.RestrictionSearcher('EcoRI', 'BsaI', 'BaeI')
        # Contains the full EcoRI site
        ecori_seq = 'GAATTC'
        self.assertTrue(rs.sitesPresent(ecori_seq))
        self.assertTrue(rs.sitesPresent(ecori_seq, full_sites=False))
        # Contains only the core BsaI site
        bsai_core_only_seq = 'GGTCTC'
        self.assertFalse(rs.sitesPresent(bsai_core_only_seq))
        self.assertTrue(rs.sitesPresent(bsai_core_only_seq, full_sites=False))
        # Contains both EcoRI and BsaI full sites, overlapping
        ecori_bsai_seq = 'GGTCTCGAATTC'
        self.assertTrue(rs.sitesPresent(ecori_bsai_seq))
        self.assertTrue(rs.sitesPresent(ecori_bsai_seq, full_sites=False))
        # Compare instance method to function
        self.assertEqual(
            rs.sitesPresent(ecori_bsai_seq), 
            restriction.sitesPresent(
                ecori_bsai_seq,
                ['EcoRI', 'BsaI', 'BaeI'],
            )
        )

    def test_countSites(self):
        ''' Test RestrictionSearcher.countSites '''
        rs = restriction.RestrictionSearcher('EcoRI', 'BsaI', 'BaeI')
        # Contains one copy of the EcoRI site (palindromic, count should be 2)
        ecori_seq = 'GAATTC'        
        counts = rs.countSites(ecori_seq)
        self.assertTrue(counts[0] == 2 and counts[1] == 0 and counts[2] == 0)
        # Contains only the core BsaI site (not palindromic, count should be 1)
        bsai_core_only_seq = 'GGTCTC'
        # Looking for full sites so all counts should be 0
        counts = rs.countSites(bsai_core_only_seq)
        self.assertTrue(counts[0] == 0 and counts[1] == 0 and counts[2] == 0)
        # Looking for core sites so BsaI count should be 1
        counts = rs.countSites(bsai_core_only_seq, full_sites=False)
        self.assertTrue(counts[0] == 0 and counts[1] == 1 and counts[2] == 0)
        # Contains both EcoRI and BsaI full sites, overlapping
        ecori_bsai_seq = 'GGTCTCGAATTC'
        counts = rs.countSites(ecori_bsai_seq)
        self.assertTrue(counts[0] == 2 and counts[1] == 1 and counts[2] == 0)
        counts = rs.countSites(ecori_bsai_seq, full_sites=False)
        self.assertTrue(counts[0] == 2 and counts[1] == 1 and counts[2] == 0)
        # Compare instance method to function
        self.assertEqual(
            rs.countSites(ecori_bsai_seq), 
            restriction.countSites(
                ecori_bsai_seq,
                ['EcoRI', 'BsaI', 'BaeI'],
            )
        )

    def test_findSites(self):
        ''' Test RestrictionSearcher.findSites '''
        rs = restriction.RestrictionSearcher('EcoRI', 'BsaI', 'BaeI')
        # Contains both EcoRI and BsaI full sites, overlapping
        ecori_bsai_seq = 'GGTCTCGAATTC' 
        sites = rs.findSites(ecori_bsai_seq)
        self.assertTrue(len(sites[0]) == 2 and len(sites[1]) == 1 and 
                        len(sites[2]) == 0)
        # EcoRI on fwd strand
        self.assertTrue(sites[0][0] == (1, 6, 12))
        # EcoRI on rev strand (indexed from 5' of fwd strand)
        self.assertTrue(sites[0][1] == (-1, 6, 12))
        # BsaI on fwd strand
        self.assertTrue(sites[1][0] == (1, 0, 11))
        # Compare instance method to function
        self.assertEqual(
            rs.findSites(ecori_bsai_seq), 
            restriction.findSites(
                ecori_bsai_seq,
                ['EcoRI', 'BsaI', 'BaeI'],
            )
        )

if __name__ == '__main__':
    unittest.main(verbosity=2)

