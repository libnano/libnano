# -*- coding: utf-8 -*-
'''libnano.padlock

Generation + filtering of padlock probes / MIPs from a target region sequence


Padlock structure reminder, left and right are in terms of the hybridized sequence

LINEAR VERSION:

5'    Right Arm       Scaffold Seq (aka Loop)      Left Arm     3'
+------------------>+-----------~-----------++------------------>

HYBRIDIZED VERSION

                 Scaffold Seq (aka Loop)
        -------------------~--------------------
        |                                      |
        <     Left Arm    3' 5'   Right Arm    +
3'      +------------------>+------------------>     5'
<----------------------------------------------------+
              copied RT'd cDNA reverse strand

'''
import logging
import random
import io
from typing import (
    List,
    Tuple,
    Dict,
    NamedTuple
)
from pprint import pprint
import os.path

import yaml
from primer3 import (
    calcTm,
    calcHeterodimerTm
)

from libnano.seqstr import reverseComplement as rc
from libnano.metric.seqmetric import gcContent
from libnano import DATASET_DIR

pjoin = os.path.join

with io.open(pjoin(DATASET_DIR, 'padlock.yaml'), 'r') as fd:
    PADLOCK_SEQS: dict = yaml.load(fd)

T2S_SEQ: str = PADLOCK_SEQS['T2S_SEQ']
SCAFFOLD_SEQ_SOLID: str = ( PADLOCK_SEQS['SCAFFOLD_SEQ_SOLID'][0] +
                            T2S_SEQ +
                            PADLOCK_SEQS['SCAFFOLD_SEQ_SOLID'][1] )
SCAFFOLD_SEQ_ILUMINA: str = (   PADLOCK_SEQS['SCAFFOLD_SEQ_ILUMINA'][0] +
                                T2S_SEQ +
                                PADLOCK_SEQS['SCAFFOLD_SEQ_ILUMINA'][1] )
ILLLUMINA_SEQ: str = PADLOCK_SEQS['ILLLUMINA_SEQ']
SCAFFOLD_SEQ_HYBRID: str = (PADLOCK_SEQS['SCAFFOLD_SEQ_HYBRID'][0] +
                            ILLLUMINA_SEQ +
                            PADLOCK_SEQS['SCAFFOLD_SEQ_HYBRID'][1] +
                            T2S_SEQ +
                            PADLOCK_SEQS['SCAFFOLD_SEQ_HYBRID'][2] )


PadHit = NamedTuple("PadHit", [ ('name0', str),         # identifier string
                                ('name1', str),         # identifier string
                                ('strand_dir', str),    # fwd or rev
                                ('genome_idx', int),    # the starting genomic index of this base
                                ('idx', int),           # the index into this sequence string
                                ('padlock_seq', str),
                                ('barcode', str),
                                ('seq_r', str),
                                ('scaffold', str),
                                ('seq_l', str)])

DEFAULT_PADLOCK_CONFIG = lambda: {
    'species': 'human',                 # target species
    'padding': 25,                      # NOTE UNUSED: Padding from 3' end of gene (nt)
    'spacing': 20,                      # Spacing between probe starts (nt)
    'arm_length': 20,                   # Padlock / MIP arm length (nt)
    'gap_size': 0,                      # MIP gap size
    'arm_gc_min': 0.4,                  # Minimum arm GC content
    'arm_gc_max': 0.6,                  # Maximum arm GC content
    'scaffold': SCAFFOLD_SEQ_HYBRID,    # Adapter sequence (for struct checks)
    'exclude_seqs': ['GGGGG'],          # Subsequences to avoid
    'structure_tm_max': 30,             # Max melting temp of secondary struct
    'keep_random_n': None,              # If not None, keep only `n` probes
    'thermo_params': {                  # Primer3 thermo params
        'mv_conc': 50,
        'dv_conc': 0,
        'dntp_conc': 0.8,
        'dna_conc': 50
    },
    'arm_tm_min': 50
}
P_PARAMS: dict = DEFAULT_PADLOCK_CONFIG()

def createScaffold(barcode: str, scaf_type: str ='solid') -> str:
    if scaf_type == 'solid':
        scaffold = SCAFFOLD_SEQ_SOLID
    elif scaf_type == 'illumina':
        scaffold = SCAFFOLD_SEQ_ILUMINA
    elif scaf_type == 'hybrid':
        scaffold = SCAFFOLD_SEQ_HYBRID
    else:
        raise ValueError("Unknown scaf_type, %s" % (scaf_type))
    return  scaffold.format(barcode=barcode,
                            armr='',
                            t2s5p='',
                            t2s3p='',
                            il5p='',
                            il3p='',
                            arml=''
                            )
# end def

def screenPadlockArms(  p_l_seq: str,
                        p_r_seq: str,
                        loop_seq: str,
                        p_params: dict,
                        do_print: bool = False) -> Tuple[bool, dict]:
    is_good = True
    tp = p_params['thermo_params']
    report = {
        'arm_gc_min_l': 0,
        'arm_gc_max_l': 0,
        'arm_gc_min_r': 0,
        'arm_gc_max_r': 0,
        'l_clamp': True,
        'tm_arm_min_l': 0,
        'tm_arm_min_r': 0,
        'ex_seq': [],
        'tm_hairpin_l': 0,
        'tm_hairpin_r': 0,
        'tm_hetero_0': 0,
        'tm_hetero_1': 0,
        'tm_hetero_2': 0
    }

    "1. GC content checks"
    p_l_gc_content = gcContent(p_l_seq)
    p_r_gc_content = gcContent(p_r_seq)
    if p_l_gc_content < p_params['arm_gc_min']:
        if do_print:
            print("\tgc content L min fail %0.3f" % p_l_gc_content)
        is_good = False
    report['arm_gc_min_l'] = p_l_gc_content
    if p_r_gc_content < p_params['arm_gc_min']:
        if do_print:
            print("\tgc content R min fail %0.3f" % p_r_gc_content)
        is_good = False
    report['arm_gc_min_r'] = p_r_gc_content
    if p_l_gc_content > p_params['arm_gc_max']:
        if do_print:
            print("\tgc content L max fail %0.3f" % p_l_gc_content)
        is_good = False
    report['arm_gc_max_l'] = p_l_gc_content
    if p_r_gc_content > p_params['arm_gc_max']:
        if do_print:
            print("\tgc content R max fail %0.3f" % p_r_gc_content)
        is_good = False
    report['arm_gc_max_r'] = p_r_gc_content


    "2. GC clamp checks"
    l_3p_check = padlockLeftArmGCClamp(p_l_seq)
    if l_3p_check > 3:
        if do_print:
            print("\tl clamp fail")
        is_good = False
    report['l_clamp'] = False

    "3. Arm Tm check"
    p_arm_tm_l = calcTm(p_l_seq, **tp)
    p_arm_tm_r = calcTm(p_r_seq, **tp)
    if p_arm_tm_l < p_params['arm_tm_min']:
        if do_print:
            print("\tArm L fail %2.3f" % p_arm_tm_l)
        is_good = False
    report['tm_arm_min_l'] = p_arm_tm_l
    if p_arm_tm_r < p_params['arm_tm_min']:
        if do_print:
            print("\tArm R fail %2.3f" % p_arm_tm_r)
        is_good = False
    report['tm_arm_min_r'] = p_arm_tm_r

    p_seq = (
        p_r_seq + loop_seq + p_l_seq
    )
    "4. Check for excluded seqs"
    ex_fail = False
    for ex_seq in p_params['exclude_seqs']:
        if ex_seq in p_seq:
            ex_fail = True
            report['ex_seq'].append(ex_seq)
            break
    if ex_fail:
        is_good = False

    "5. Secondary structure / primer dimer checks"
    p_het_tm_0 = calcHeterodimerTm(p_l_seq, p_r_seq, **tp)
    p_het_tm_1 = calcHeterodimerTm(p_l_seq, loop_seq, **tp)
    p_het_tm_2 = calcHeterodimerTm(p_r_seq, loop_seq, **tp)
    if p_het_tm_0 > p_params['structure_tm_max']:
        if do_print:
            print("\thetero 0 fail")
        is_good = False
    report['tm_hetero_0'] = p_het_tm_0
    if p_het_tm_1 > p_params['structure_tm_max']:
        if do_print:
            print("\thetero 1 fail")
        is_good = False
    report['tm_hetero_1'] = p_het_tm_1
    if p_het_tm_2 > p_params['structure_tm_max']:
        if do_print:
            print("\thetero 2 fail")
        is_good = False
    report['tm_hetero_2'] = p_het_tm_2
    return is_good, report
# end def

def splitHitList(   items: List[Tuple[int, dict]],
                    arm_length,
                    spacing) -> List[List[Tuple[int, dict]]]:
    '''Split hits into groups by a spacing and an arm_length
    '''
    # split into groups by spacing
    if len(items) > 0:
        delta: int = items[0][0] + 2*arm_length + spacing
        group: List[Tuple[int, dict]] = []
        hit_lists: List[List[Tuple[int, dict]]] = [group]
        for i, report in items:
            if i > delta:
                group = []
                hit_lists.append(group)
                # increment delta for next group
                delta = i + 2*arm_length + spacing
            group.append((i, report))
        return hit_lists
    else:
        return []
# end def

def sortHitList(items: List[Tuple[int, dict]]) -> List[Tuple[int, dict]]:
    '''Sort max sum of arm Tms
    '''
    max_tm_f = lambda x: x[1]['tm_arm_min_l'] + 0.9*x[1]['tm_arm_min_r']
    return sorted(items, key=max_tm_f, reverse=True)
# end def

def writePadlocksToCSV(padlock_results: Dict[str, List[PadHit]], filename: str):
    '''Write padlocks to to a CSV file
    '''
    tp = P_PARAMS['thermo_params']
    with io.open(filename, 'w') as fd:
        fd.write(   'gene_name, name0, name1, strand_dir, genome_idx, index, '
                    'gap_size, sequence, barcode, right_arm, scaffold, '
                    'left_arm, right_tm, left_tm\n')
        temp = '%s, %s, %s, %s, %d, %d, %d, %s, %s, %s, %s, %s, %2.3f, %2.3f\n'
        for gene, seq_list in padlock_results.items():
            for seq_tuple in seq_list:
                seq_r, seq_l = seq_tuple.seq_r, seq_tuple.seq_l
                tm_tuple = (calcTm(seq_r, **tp), calcTm(seq_l, **tp))
                fd.write(temp % ((gene,) + seq_tuple + tm_tuple) )
    print('Wrote padlocks to %s' % filename)
# end def

def generatePadlocks(seq: str,
                    name0: str,
                    name1: str,
                    strand_dir: str,
                    barcodes: List[str],
                    genome_idx: int = -1,
                    arm_length: int = 20,
                    params: dict = None,
                    do_print: bool = False
                    ) -> List[PadHit]:
    '''Screen poly G's in the padlock to avoid warnings from IDT synthesis

    Args:
        seq: the length of a sequence
        name0: name to identify the padlocks by.  For example, you could use
            the Ensembl Transcript ID
        name1: name to identify the padlocks by.  For example, you could use
            the Ensembl Exon ID
        barcodes: list of one or more barcodes to try
        genome_idx: start index of a sequence in the genome
        arm_length: the length of a padlock arm
        params: default is P_PARAMS. parameters for padlock screening. Add in
            things like ``gap_size`` here
        do_print: debug the design printing
    '''
    if params is None:
        params = P_PARAMS
    else:
        default_p = DEFAULT_PADLOCK_CONFIG()
        default_p.update(params)
        params = default_p

    arm_length2: int =  2*arm_length
    gap_size: int =     params['gap_size']
    spacing: int =      params['spacing']

    scaffold: str = None
    if len(barcodes) == 0 or not isinstance(barcodes, (tuple, list)):
        raise ValueError("barcodes length must be non-zero")
    for barcode in barcodes:
        candidate_scaffold: str = createScaffold(barcode, scaf_type='hybrid')
        if 'GGGG' not in candidate_scaffold:
            scaffold = candidate_scaffold
    if scaffold is None:
        raise ValueError('polyG in scaffold for all barcodes')

    items = []
    search_range = range(len(seq) - arm_length2)
    for i in search_range:
        if 'GGGG' not in seq[i:i + arm_length2 + gap_size]: # including the gap_size
            l_primer = seq[i:i + arm_length]
            r_primer = seq[i + arm_length + gap_size:i + arm_length2 + gap_size]
            is_good, report = screenPadlockArms(    l_primer,
                                                    r_primer,
                                                    scaffold,
                                                    params)
            if is_good:
                '''add the start index of the padlock and the report dictionary
                to the items list'''
                items.append((i, report))
            # elif do_print:
            #     print("FAILURE")
            #     pprint(report)
    # end for
    # print(items)
    hit_lists = splitHitList(items, arm_length=arm_length, spacing=spacing)
    hit_lists = [sortHitList(x) for x in hit_lists]

    if do_print:
        print('The number of hits', len(hit_lists))

    # pick the first element in each group
    sampled_list: List = [x[0] for x in hit_lists]
    if do_print:
        print('$HIT_COUNT:', [len(x) for x in hit_lists])
    sequences_list: List[PadHit] = []
    for i, b in sampled_list:

        seq_l = seq[i:i + arm_length]
        seq_r = seq[i + arm_length + gap_size:i + arm_length2 + gap_size]

        if do_print:
            print("%d,\t %2.3f, %2.3f" % (i, b['tm_arm_min_l'], b['tm_arm_min_r']))
            print("%s, %s" % (seq_l, seq_r))
            print("%s" % (seq[i:i + arm_length2]))

        sequences_list.append(PadHit(
                                    name0,
                                    name1,
                                    strand_dir,
                                    genome_idx,
                                    i,
                                    gap_size,
                                    seq_r + scaffold + seq_l,
                                    barcode,
                                    seq_r,
                                    scaffold,
                                    seq_l))
    # end for
    return sequences_list
# end def

def padlockRightArmGCClamp(p: str) -> int:
    r_3p = p[-5:]
    r_3p_check = r_3p.count('G') + r_3p.count('C')
    return r_3p_check
# end def

def padlockLeftArmGCClamp(p: str) -> int:
    l_3p = p[0:5]
    l_3p_check = l_3p.count('G') + l_3p.count('C')
    return l_3p_check
# end def

