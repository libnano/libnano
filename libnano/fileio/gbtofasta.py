#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Rips the sequence out of a genbank file and outputs a fasta file:

<filename.gb> -> <filename.fa>

'''

import os

def genbankToFasta(genbank_fp: str, seq_record: dict = None):
    if seq_record:
        primary_rec = seq_record
    else:
        rec = gb_reader.parse(seq_fp)
        primary_rec = str(rec['seq'])

    genbank_fp_no_ext = os.path.splitext(genbank_fp)[0]
    output_fp = genbank_fp_no_ext + '.fa'
    with open(output_fp, 'w', encoding='utf-8') as output_fd:
        output_fd.write('>{}\n'.format(genbank_fp_no_ext.split('/\\')[-1]))
        output_fd.write(str(primary_rec['seq']))
        output_fd.write('\n')
    return

if __name__ == '__main__':
    import sys
    genbank_fp = sys.argv[1]
    genbankToFasta(genbank_fp)
