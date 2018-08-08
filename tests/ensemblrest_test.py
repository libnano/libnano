# -*- coding: utf-8 -*-
import unittest
import os
from typing import (
    Set
)

import conftest

import libnano.ensemblrest as er

if os.environ.get('IS_TRAVIS') is None:
    class TestEnsembleRest(unittest.TestCase):
        def test_getProbes(self):
            GENE: str = 'SLC17A8'
            TRANSCRIPT: str = 'ENST00000323346.9'
            PROBE_NAME: str = 'ILMN_1767842'
            GENE_EID: str = 'ENSG00000179520'

            out_list: list = er.getArrayProbes(TRANSCRIPT)
            self.assertIsInstance(out_list, list)
            self.assertTrue(len(out_list) > 0)
            self.assertIsInstance(out_list[0], dict)
            out_unique: Set[str] = set()
            for x in out_list:
                p_hash: str = er._probe_dict_2_str(x)
                if p_hash in out_unique:
                    raise ValueError("non-unique_hash: %s" % (p_hash))
                out_unique.add(p_hash)
            grouped_out_list: list = er.probeListGroupByProbeName(out_list)
            self.assertTrue(len(grouped_out_list) > 0)

            probe: dict = er.getProbeFromList(PROBE_NAME, grouped_out_list)
            self.assertIsInstance(probe, dict)

        def test_getRegionSequence(self):
            REF_SEQ: str = 'ATCCATGCAAGCCCCATAAAACAGTTCCTAGCATGCAGAAAATGCCCACG'
            SPECIES: str = 'human'
            GENE: str = 'SLC17A8'
            GENE_EID: str = 'ENSG00000179520'
            TRANSCRIPT: str = 'ENST00000323346.9'
            THREE_P_EXON: str = 'ENSE00001244923'

            PROBE: dict = {
                'probe_length': 50,
                'probe_set': '',
                'feature_type': 'array_probe',
                'end': 100421903,
                'seq_region_name': '12',
                'strand': 1,
                'probe_name': 'ILMN_1767842',
                'start': 100421854,
                'microarrays': ['HumanWG_6_V2', 'HumanWG_6_V3', 'HumanHT-12_V3', 'HumanHT-12_V4', 'HumanRef-8_V3']
            }
            # 1. Get the probe sequence
            seq: str = er.getRegionSequence(SPECIES,
                                chromosome=PROBE['seq_region_name'],
                                start_idx=PROBE['start'],
                                end_idx=PROBE['end'],
                                strand=PROBE['strand']
                                )
            self.assertEqual(REF_SEQ, seq)

            # 2. confirm that getRegionSequence matches slicing the whole gene
            # for a forward strand (strand == 1) (no reverse complement check)
            lookup: dict = er.lookUpID(GENE_EID)
            gene_start_idx: int = lookup['start']
            full_gene_seq: str = er.getSequence(GENE_EID, seq_type='genomic')
            sl: slice = slice(          PROBE['start'] - gene_start_idx,
                                        PROBE['end'] - gene_start_idx + 1)
            slice_seq: str = full_gene_seq[sl]

            self.assertEqual(slice_seq, REF_SEQ)

            # 3. Check that the sequence exists in the 3' exon
            three_p_seq: str = er.getSequence(THREE_P_EXON)
            self.assertTrue(len(three_p_seq) > 0)
            self.assertIn(REF_SEQ, three_p_seq)
        # end def
    # end class

if __name__ == '__main__':
    tr = TestEnsembleRest()
    tr.test_getRegionSequence()


