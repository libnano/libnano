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
tests.ensemblrest_test
~~~~~~~~~~~~~~~~~~~~~~

'''
import os
import time
import unittest
from typing import (
    List,
    NamedTuple,
    Set,
)

# import conftest
import pandas as pd
import pytest

import libnano.ensemblrest as er


@pytest.mark.skipif(os.environ.get('IS_TRAVIS') is None)
class TestEnsembleRest(unittest.TestCase):
    def test_get_probes(self):
        # gene = 'SLC17A8'
        transcript = 'ENST00000323346.9'
        probe_name = 'ILMN_1767842'
        # gene_eid = 'ENSG00000179520'

        out_list = er.getArrayProbes(transcript)
        self.assertIsInstance(
            out_list,
            list,
        )
        self.assertTrue(
            len(out_list) > 0,
        )
        self.assertIsInstance(
            out_list[0],
            dict,
        )
        out_unique: Set[str] = set()

        for x in out_list:
            p_hash = er._probe_dict_2_str(x)
            if p_hash in out_unique:
                raise ValueError(
                    f'non-unique_hash: {p_hash}',
                )
            out_unique.add(p_hash)
        grouped_out_list = er.probeListGroupByProbeName(
            out_list,
        )
        self.assertTrue(
            len(grouped_out_list) > 0,
        )

        probe_dict = er.getProbeFromList(
            probe_name,
            grouped_out_list,
        )
        self.assertIsInstance(
            probe_dict,
            dict,
        )

    def test_getRegionSequence(self):
        ref_seq = 'ATCCATGCAAGCCCCATAAAACAGTTCCTAGCATGCAGAAAATGCCCACG'
        species = 'human'
        # gene = 'SLC17A8'
        gene_eid = 'ENSG00000179520'
        # transcript = 'ENST00000323346.9'
        three_p_exon = 'ENSE00001244923'

        probe_dict = {
            'probe_length': 50,
            'probe_set': '',
            'feature_type': 'array_probe',
            'end': 100421903,
            'seq_region_name': '12',
            'strand': 1,
            'probe_name': 'ILMN_1767842',
            'start': 100421854,
            'microarrays': [
                'HumanWG_6_V2',
                'HumanWG_6_V3',
                'HumanHT-12_V3',
                'HumanHT-12_V4',
                'HumanRef-8_V3',
            ],
        }
        # 1. Get the probe sequence
        seq = er.getRegionSequence(
            species,
            chromosome=probe_dict['seq_region_name'],
            start_idx=probe_dict['start'],
            end_idx=probe_dict['end'],
            strand=probe_dict['strand'],
        )
        self.assertEqual(ref_seq, seq)

        # 2. confirm that getRegionSequence matches slicing the whole gene
        # for a forward strand (strand == 1) (no reverse complement check)
        lookup: dict = er.lookUpID(gene_eid)
        gene_start_idx: int = lookup['start']
        full_gene_seq = er.getSequence(gene_eid, seq_type='genomic')
        slc = slice(
            probe_dict['start'] - gene_start_idx,
            probe_dict['end'] - gene_start_idx + 1,
        )
        slice_seq = full_gene_seq[slc]

        self.assertEqual(
            slice_seq,
            ref_seq,
        )

        # 3. Check that the sequence exists in the 3' exon
        three_p_seq = er.getSequence(
            three_p_exon,
        )
        self.assertTrue(
            len(three_p_seq) > 0,
        )
        self.assertIn(
            ref_seq,
            three_p_seq,
        )

    def test_get_probes_for_id(self):
        species = 'human'

        class GeneInfo(NamedTuple):
            symbol: str
            barcode: str
            gene_id: str
            transcript_id: str
            utr_id: str

        genes_and_barcodes_list = [
            GeneInfo(
                symbol='SLC17A8',
                barcode='ACAGC',
                gene_id='ENSG00000179520',
                transcript_id='ENST00000323346',
                utr_id='ENSE00001244923',
            ),
            GeneInfo(
                symbol='GFAP',
                barcode='TACAT',
                gene_id='ENSG00000131095',
                transcript_id='ENST00000638281',
                utr_id='ENSE00003806990',
            ),
        ]
        number_probes_per_gene = 5
        sleep_time = 0.1

        probe_cutoff_length = 36

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
            'array_freq',
        ]
        df_out: pd.DataFrame = pd.DataFrame(
            columns=out_columns,
        )

        def todfdict(x):
            return {a: b for a, b in zip(out_columns, x)}

        for item in genes_and_barcodes_list:
            three_p_exon_id = item.utr_id
            transcript_id = item.transcript_id
            filtered_probes: pd.DataFrame = er.getProbesForID(
                three_p_exon_id,
                keep_n=number_probes_per_gene,
            )
            self.assertTrue(len(filtered_probes) > 0)
            barcode = item.barcode
            exon_seq = er.getSequence(three_p_exon_id)
            transcript: dict = er.lookUpID(transcript_id)

            for exon in transcript['Exon']:
                if exon['id'] == three_p_exon_id:
                    exon_start_idx: int = exon['start']
                    # exon_end_idx: int = exon['end']

            for i in range(len(filtered_probes)):
                probe = filtered_probes.iloc[i]
                p_start: int = probe['start']
                p_end: int = probe['end']
                p_strand: int = probe['strand']
                p_length: int = probe['probe_length']

                # Sometimes p_length doesn't match start and end indices
                # so let us filter those out
                if (p_end - p_start + 1) > p_length:
                    continue
                elif p_length < probe_cutoff_length:
                    # Try extending from the 5' end of the sequences
                    delta: int = probe_cutoff_length - p_length
                    if p_strand == 1:
                        if (
                            p_start > exon_start_idx and
                            (p_start - exon_start_idx) > delta
                        ):
                            p_start -= delta
                        else:
                            continue
                    else:
                        if (
                            p_start < exon_start_idx and
                            (exon_start_idx - p_start) > delta
                        ):
                            p_start += delta  # Increase the index
                        else:
                            continue

                try:
                    seq = er.getRegionSequence(
                        species,
                        chromosome=probe['seq_region_name'],
                        start_idx=p_start,
                        end_idx=p_end,
                        strand=p_strand,
                    )
                except ValueError:
                    raise
                try:
                    seq, _ = er.filterRegionSequence(
                        seq,
                        p_strand,
                        transcript_id,
                        transcript,
                        exon_seq,
                    )
                except ValueError:
                    continue  # Skip this probe
                row = [
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
                    probe['array_freq'],
                ]
                df_out = df_out.append(
                    todfdict(row),
                    ignore_index=True,
                )
                # Sleep such that requests don't time out
                time.sleep(sleep_time)

        self.assertTrue(len(df_out) > 0)


if __name__ == '__main__':
    tr = TestEnsembleRest()
    tr.test_getRegionSequence()
