# -*- coding: utf-8 -*-
import os.path as _ospath

try:
    import libnano
except:
    import sys
    from os.path import dirname, abspath
    MOD_PATH = dirname(dirname(dirname(abspath(__file__))))
    sys.path = [MOD_PATH] + sys.path

from libnano.fileio import (
    gbtofasta,
    xmfa,
    fasta,
    gb_reader,
    gb_writer
)

def getSeqFromFile(seq_fp):
    _, ext = _ospath.splitext(seq_fp)
    if ext in ('.gb', '.gbk', '.genbank'):
        rec = gb_reader.parse(seq_fp)
        return str(rec['seq'])
    elif ext in ('.fa', '.fasta'):
        # Just returns the first sequence in the file for now --
        # this may be undesirable as a general behavior
        return fasta.parseFasta(seq_fp)[0][1]

__all__ = ['fasta', 'gbtofasta', 'xmfa', 'getSeqFromFile',
            'gb_reader', 'gb_writer']