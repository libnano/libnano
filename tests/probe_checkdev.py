from typing import (
    List,
    Tuple
)
import pprint

import pandas as pd

import conftest
import libnano.ensemblrest as er

GENE: str = 'SLC17A8'
GENE_EID: str = 'ENSG00000179520'
TRANSCRIPT: str = 'ENST00000323346.9'
THREE_P_EXON: str = 'ENSE00001244923'
PROBE_NAME: str = 'ILMN_1767842'
SPECIES: str = 'human'



filtered_probes: pd.DataFrame = er.getProbesForID(THREE_P_EXON, keep_n=5)
for i in range(len(filtered_probes)):
    probe = filtered_probes.iloc[i]
    print(probe)
    seq: str = er.getRegionSequence(SPECIES,
                                    chromosome=probe['seq_region_name'],
                                    start_idx=probe['start'],
                                    end_idx=probe['end'],
                                    strand=probe['strand']
                                    )
    print(seq)
# end for
GENES_AND_BARCODES: List[Tuple[str, str]] = [
    ('SLC17A8', 'ACAGC', 'ENST00000323346', 'ENSE00001244923'),
    ('GFAP',    'TACAT', 'ENST00000638281', 'ENSE00003806990'),
    ('FOXO1',   'TTTGC', 'ENST00000379561', 'ENSE00001481591'),
    ('PSEN2',   'CATTA', 'ENST00000366783', 'ENSE00001380688'),
    ('DAXX',    'AACCG', 'ENST00000374542', 'ENSE00001815438'),
    ('CDK5R1',  'CGAGA', 'ENST00000313401', 'ENSE00001271015')
]

GENE_LIST = [x[0] for x in GENES_AND_BARCODES]
print(er.getThreePrimeUTRs(SPECIES, GENE_LIST))

