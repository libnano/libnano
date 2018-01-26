import sys

import random
import unittest
import time
from os.path import join, abspath, dirname

# For package imports
sys.path.insert(0, abspath(join(dirname(__file__), '..')))

import _setup
from libnano.core.seqmetric import seqrepeat
from libnano.fileio import fasta

LOCAL_DIR = abspath(dirname(__file__))

def foo():
    a = (12,12323,3214)
    return a

def test_windowing():
    a = seqrepeat.RepeatCheck('GGGGAAAATTTTTTATTTTTGGGGAAAAT', 8)
    print("1:", a.indicesOf('AAAATTTT'))
    print("2:", a.indicesOf('GGGGAAAAT'))

def testGC():
    seq = ('CATTACTCTAAGTAAGCTGGTCTGTCGCCA' +   # repeat 1 30bp
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
            'CATTACTCTAAGTAAGCTGGTCTGTCGCCA' + # repeat 2 30 bp
            "GGCCGAAGGGGGTGACCGACAAACGGCGGCGGGATTTCGCGCCAGGCGTAGTGTAACATGACCGCTAGGA"
        )
    a = seqrepeat.RepeatCheck(seq, 8)
    check = a.gapCheck(70)
    i,j,size = check[0]
    print(check, seq[i:i + size], seq[j:j + size])
    check = a.gapCheck(72)
    print(check)

def testRepeats():
    seq = ('CATTACTCTAAGTAAGCTGGTCTGTCGCCA' +   # repeat 1 30bp
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
            'CATTACTCTAAGTAAGCTGGTCTGTCGCCA' + # repeat 2 30 bp
            "GGCCGAAGGGGGTGACCGACAAACGGCGGCGGGATTTCGCGCCAGGCGTAGTGTAACATGACCGCTAGGA"
        )
    a = seqrepeat.RepeatCheck(seq, 8)
    check = a.allRepeats()
    print(check)

def testScreenEnds():
    key = 'CATTACTCTAAGTAAGCTGGTCTGTCGCCA'
    # (True, True)
    seq = ( key +   # repeat 1 30bp
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
            key # repeat 2 30 bp
        )
    a = seqrepeat.RepeatCheck(seq, 8)
    check = a.screenEnds(8)
    print("screen ends", check)
    # (False, True)
    seq = ( "AAACCCCGC" + key +   # repeat 1 30bp
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
            key # repeat 2 30 bp
        )
    a = seqrepeat.RepeatCheck(seq, 8)
    check = a.screenEnds(8)
    print("screen ends", check)
    # (False, False)
    seq = ( "AAACCCCGC" + key +   # repeat 1 30bp
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
            key + "GGGAACCTTTA" # repeat 2 30 bp
        )
    a = seqrepeat.RepeatCheck(seq, 8)
    check = a.screenEnds(8)
    print("screen ends", check)
    # (True, False)
    seq = ( key +   # repeat 1 30bp
            'CCTTATAGTTTACGAACCAGTTATGTGATATTGCCTGATCT' + # space 41 bp
            key + "GGGAACCTTTA" # repeat 2 30 bp
        )
    a = seqrepeat.RepeatCheck(seq, 8)
    check = a.screenEnds(8)
    print("screen ends", check)


if __name__ == '__main__':
    # testGC()
    # testRepeats()
    testScreenEnds()