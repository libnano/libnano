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
tests.probe_checkdev
~~~~~~~~~~~~~~~~~~~~

'''
# import pprint
from typing import (  # Tuple,
    List,
    NamedTuple,
)

# import conftest
import pandas as pd

import libnano.ensemblrest as er

GENE: str = 'SLC17A8'
GENE_EID: str = 'ENSG00000179520'
TRANSCRIPT: str = 'ENST00000323346.9'
THREE_P_EXON: str = 'ENSE00001244923'
PROBE_NAME: str = 'ILMN_1767842'
SPECIES: str = 'human'


class GeneInfo(NamedTuple):
    symbol: str
    barcode: str
    gene_id: str
    transcript_id: str
    utr_id: str


GENES_AND_BARCODES_LIST = [
    GeneInfo(
        'SLC17A8', 'ACAGC', 'ENSG00000179520',
        'ENST00000323346', 'ENSE00001244923',
    ),
    GeneInfo(
        'GFAP', 'TACAT', 'ENSG00000131095',
        'ENST00000638281', 'ENSE00003806990',
    ),
    GeneInfo(
        'FOXO1', 'TTTGC', 'ENSG00000150907',
        'ENST00000379561', 'ENSE00001481591',
    ),
    GeneInfo(
        'PSEN2', 'CATTA', 'ENSG00000143801',
        'ENST00000366783', 'ENSE00001380688',
    ),
    GeneInfo(
        'DAXX', 'AACCG', 'ENSG00000204209',
        'ENST00000374542', 'ENSE00001815438',
    ),
    GeneInfo(
        'CDK5R1', 'CGAGA', 'ENSG00000176749',
        'ENST00000313401', 'ENSE00001271015',
    ),
]

GENE_LIST = [x.symbol for x in GENES_AND_BARCODES_LIST]
print(er.getThreePrimeUTRs(SPECIES, GENE_LIST))

out_columns = [
    'symbol',
    'gene_id',
    'exon_id'
    'probe_name',
    'probe_seq',
    'probe_start',
    'probe_end',
    'probe_strand',
]
df_out: pd.DataFrame = pd.DataFrame(
    columns=out_columns,
)


def todfdict(x):
    return {a: b for a, b in zip(out_columns, x)}


for item in GENES_AND_BARCODES_LIST:
    three_p_exon_id: str = item.utr_id
    filtered_probes: pd.DataFrame = er.getProbesForID(
        three_p_exon_id, keep_n=5,
    )
    for i in range(len(filtered_probes)):
        probe = filtered_probes.iloc[i]
        print(probe)
        p_start: int = probe['start']  # type: ignore
        p_end: int = probe['end']  # type: ignore
        p_strand: int = probe['strand']  # type: ignore
        seq: str = er.getRegionSequence(
            SPECIES,
            chromosome=probe['seq_region_name'],   # type: ignore
            start_idx=p_start,
            end_idx=p_end,
            strand=probe['strand'],  # type: ignore
        )
        print(seq)
        row: list = [
            item.symbol,
            item.gene_id,
            item.utr_id,
            probe['probe_name'],
            seq,
            p_start,
            p_end,
            p_strand,
        ]
        df_out = df_out.append(   # type: ignore
            todfdict(row),
            ignore_index=True,
        )
