from __future__ import print_function

import sys
import os
import re
LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(os.path.dirname(LOCAL_DIR),'tests', 'test_data')
fn = os.path.join(data_path, 'NEB_restriction_enzymes.txt')

"""
There should be

{ 'name': { 'base_seq': [ fwd, revcomp],
            'expanded_seq': [ fwd, revcomp],
            'base_regex': [ fwd, revcomp],
            'expanded_regex': [ fwd, revcomp],
            'cutsites_idxs' [fwd, rev],
            'synonyms': []
            'alias': [names],
            'shorthand': <NEB format>,
        },
}
"""


enzyme_re_raw = ('(?:\((?P<s1>[\d|-]+)/(?P<s2>[\d|-]+)\))?'
                 '(?P<seq>[^\(]+)'
                 '(?:\((?P<e1>[\d|-]+)/(?P<e2>[\d|-]+)\))?')
enzyme_re = re.compile(enzyme_re_raw)

def convertRestrictionRe(raw_data):
    re_res = re.search(enzyme_re, raw_data)
    raw_seq = re_res.group('seq')
    seq = raw_seq.replace('/', '')
    if re_res.group('s1') is not None:
        if re_res.group('e1') is not None:
            idxs_5p = [int(re_res.group('s1')), int(re_res.group('s2'))]
            idxs_3p = [int(re_res.group('e1')), int(re_res.group('e2'))]
            return seq, (idxs_5p, idxs_3p)
        else:
            return seq, [int(re_res.group('s1')), int(re_res.group('s2'))]
    elif re_res.group('e1') is not None:
            return seq, [int(re_res.group('e1')), int(re_res.group('e2'))]
    else:
        # no end stuff vanilla
        try:
            idx = raw_seq.index('/')
        except ValueError:
            # no / in the sequence for some reason like Nb.BbvCI
            idx = None
        return (seq, idx)


def convertRestriction(seq):
    start, mid, end = seq.partition(')')
    if mid == '':
        # no end stuff vanilla
        try:
            idx = start.index('/')
        except ValueError:
            # no / in the sequence for some reason like Nb.BbvCI
            return (seq, None)
        if idx != len(start)-1:
            seq = start[:idx]+seq[idx+1:]
        else:
            seq = seq[:-1]
        return (seq, idx)
    elif end == '':
        # single ended like ACNNNNGTAYC(12/7)
        seq, _, idxs_3p = start.partition('(') # get the second end
        idxs_3p = idxs_3p.rstrip(')')
        idxs_3p = idxs_3p.split('/')
        idxs_3p = [int(x) for x in idxs_3p]
        return (seq, idxs_3p)
    else: # double ended like (10/15)ACNNNNGTAYC(12/7)
        idxs_5p = start.lstrip('(')
        idxs_5p = idxs_5p.split('/')
        idxs_5p = [int(x) for x in idxs_5p]

        seq, sep, idxs_3p = end.partition('(') # get the second end
        if sep == '':
            idxs_3p = None # 5p only
        else:
            idxs_3p = idxs_3p.rstrip(')')
            idxs_3p = idxs_3p.split('/')
            idxs_3p = [int(x) for x in idxs_3p]
        return (seq, (idxs_5p, idxs_3p))
# end def

if __name__ == '__main__':
    # fn = sys.argv[1]
    d = {}
    with open(fn, 'r') as fd:
        lines = fd.readlines()
    for line in lines:
        # print(line)
        line = line.rstrip('\n')
        items = line.split(' ')
        seq = items[0]
        names = items[1:]
        for name in names:
            d[name] = convertRestrictionRe(seq)
    print(d)
    from libnano.matcher import visualizeCutSite
    for re_name, (re_seq, idx_list) in d.items():
        print(re_name, re_seq, idx_list)
        visualizeCutSite(re_seq, idx_list)
        print()
