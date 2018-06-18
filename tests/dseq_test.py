# -*- coding: utf-8 -*-
import unittest
import _setup
from libnano.datastructures.dseq import DSeq

class TestDSeq(unittest.TestCase):

    def test_DSeqSingle(self):
        dseq = DSeq("GGATCCAAA")
        print(dseq)

    def test_DSeqAdd(self):
        dseq = DSeq("GGATCCAAA")
        print(dseq + dseq)

if __name__ == '__main__':
    x = TestDSeq()
    x.test_DSeqSingle()
    x.test_DSeqAdd()