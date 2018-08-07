from typing import (
    List
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

out_list: list = er.getProbes(THREE_P_EXON)
# pprint.pprint(out_list)
# grouped_out_list: list = er.probeListGroupByProbeName(out_list)
# pprint.pprint(grouped_out_list)
# out_max_idx: int = -1
# max_len: int = 0
# for i, item in enumerate(grouped_out_list):
#     len_item: int = len(item['microarrays'])
#     if len_item > max_len:
#         max_len = len_item
#         out_max_idx = i
# print(max_len)
# print(grouped_out_list[out_max_idx])

columns = [
'probe_name',
'start',
'end',
'probe_length',
'feature_type',
'probe_set',
'seq_region_name',
'strand',
'microarray'
]
df: pd.DataFrame = pd.DataFrame(out_list, columns=columns)
DataFrameGroupBy_T = pd.core.groupby.groupby.DataFrameGroupBy
grouped: DataFrameGroupBy_T = df.groupby('probe_name')
count_of_probe_use: pd.Series = grouped.size()
print(count_of_probe_use.max())
largest = count_of_probe_use.nlargest(5)
# probe: dict = er.getProbeFromList(PROBE_NAME, grouped_out_list)


filtered_probes: pd.DataFrame = df[df['probe_name'].isin(largest.index.values)]
columns_to_keep: List[str] = [
    'probe_name',
    'start',
    'end',
    'strand',
    'probe_length',
    'seq_region_name'
]
filtered_probes = filtered_probes.loc[:, columns_to_keep].drop_duplicates()
print(filtered_probes)
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

# print(grouped.get_group(PROBE_NAME))
# er.PROBE_KEYS

# pk = (
# 'end',
# 'feature_type',
# 'microarray',
# 'probe_length',
# 'probe_name',
# 'probe_set',
# 'seq_region_name',
# 'start',
# 'strand'
# )
# # df['probe_name'].tolist()
# indx = pd.MultiIndex.from_tuples(df['probe_name'].tolist(), names=['probe_name'])
# cols = pk
# df = pd.DataFrame(out_list, df['probe_name'].tolist(), cols)
# print(df.head())
# piv = df.pivot(index='probe_name', columns=['start', 'end', 'strand'])
# print(piv)
pndf = df['probe_name']
# print(df.sort_values(by='probe_name')['probe_name'])
# ATCCATGCAAGCCCCATAAAACAGTTCCTAGCATGCAGAAAATGCCCACG
