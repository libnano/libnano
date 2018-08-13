# -*- coding: utf-8 -*-
import unittest
import os
from typing import (
    Set,
    List,
    NamedTuple
)
import time

import pandas as pd

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

        def test_getProbesForID(self):
            SPECIES: str = 'human'
            GeneInfo = NamedTuple('GeneInfo', [
                ('symbol', str),
                ('barcode', str),
                ('gene_id', str),
                ('transcript_id', str), # canonical id
                ('utr_id', str)  #
                ]
            )
            GENES_AND_BARCODES: List[GeneInfo] = [
                GeneInfo('SLC17A8', 'ACAGC', 'ENSG00000179520',
                    'ENST00000323346', 'ENSE00001244923'),
                GeneInfo('GFAP', 'TACAT', 'ENSG00000131095',
                    'ENST00000638281', 'ENSE00003806990')
            ]
            NUMBER_PROBES_PER_GENE: int = 5
            SLEEP_TIME: float = 0.1

            PROBE_CUTOFF_LENGTH: int = 36

            out_columns: List[str] = [
                'symbol',
                'gene_id',
                'exon_id',
                'probe_name',
                'probe_seq',
                'probe_start',
                'probe_end',
                'probe_strand',
                'probe_length',
                'barcode',
                'array_freq'
            ]
            df_out: pd.DataFrame = pd.DataFrame(columns=out_columns)

            todfdict = lambda x: {a: b for a, b in zip(out_columns, x)}

            for item in GENES_AND_BARCODES:
                three_p_exon_id: str = item.utr_id
                transcript_id: str = item.transcript_id
                filtered_probes: pd.DataFrame = er.getProbesForID(three_p_exon_id,
                    keep_n=NUMBER_PROBES_PER_GENE)
                self.assertTrue(len(filtered_probes) > 0)
                barcode: str = item.barcode
                exon_seq: str = er.getSequence(three_p_exon_id)
                transcript: dict = er.lookUpID(transcript_id)

                for exon in transcript['Exon']:
                    if exon['id'] == three_p_exon_id:
                        exon_start_idx: int  = exon['start']
                        exon_end_idx: int  = exon['end']


                for i in range(len(filtered_probes)):
                    probe = filtered_probes.iloc[i]
                    p_start: int =  probe['start']
                    p_end: int =    probe['end']
                    p_strand: int = probe['strand']
                    p_length: int = probe['probe_length']

                    # sometimes p_length doesn't match start and end indices so let's filter those out
                    if (p_end - p_start + 1) > p_length:
                        print(probe['probe_name'])
                        continue
                    elif p_length < PROBE_CUTOFF_LENGTH:
                        # Try extending from the 5' end of the sequences
                        delta: int = PROBE_CUTOFF_LENGTH - p_length
                        if p_strand == 1:
                            if (p_start > exon_start_idx and
                                (p_start - exon_start_idx) > delta):
                                p_start -= delta
                            else:
                                continue
                        else:
                            if (p_start < exon_start_idx and
                                (exon_start_idx - p_start) > delta):
                                p_start += delta # increase the index
                            else:
                                continue

                    try:
                        seq: str = er.getRegionSequence(
                            SPECIES,
                            chromosome=probe['seq_region_name'],
                            start_idx=p_start,
                            end_idx=p_end,
                            strand=p_strand
                        )
                    except:
                        print(item.symbol, transcript_id, probe['probe_name'])
                        print(p_start, p_end)
                        raise
                    try:
                        was_rc: bool
                        seq, was_rc = er.filterRegionSequence(
                                seq,
                                p_strand,
                                transcript_id,
                                transcript,
                                exon_seq
                        )
                    except ValueError as ex:
                        continue # skip this probe
                    row: list = [
                        item.symbol,
                        item.gene_id,
                        item.utr_id,
                        probe['probe_name'],
                        seq,
                        p_start,
                        p_end,
                        p_strand,
                        probe['probe_length'],
                        barcode,
                        probe['array_freq']
                    ]
                    df_out = df_out.append(todfdict(row), ignore_index=True)
                    time.sleep(SLEEP_TIME)  # sleep such that requests don't time out
                # end for
            # end for
            self.assertTrue(len(df_out) > 0)
        # end def
    # end class

if __name__ == '__main__':
    tr = TestEnsembleRest()
    tr.test_getRegionSequence()


