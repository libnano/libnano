# -*- coding: utf-8 -*-
import unittest
import _setup
from libnano.datastructures.dseq import DSeq
from libnano.seqstr import reverse

class TestDSeq(unittest.TestCase):

    def test_single(self):
        dseq: DSeq = DSeq("GGATCCAAA")

    def test_add(self):
        dseq: DSeq = DSeq("GGATCCAAA")
        check1 = dseq + dseq
        dseq: DSeq = DSeq("GGATCCAAAG", 'TTTGGATCCC', overhang=1)
        check2 = dseq + dseq

    def test_cut(self):
        # Test single cutter BsaI
        fwd: str = 'GGTCTCGAATTCAAA'
        rev: str = 'GAATTCGAGACCAAA'

        bsai_cut_fwd_5prime: str = 'GGTCTCG'
        bsai_cut_fwd_3prime: str = 'AATTCAAA'

        bsai_cut_rev_3prime: str = 'AATTCGAGACCAAA'
        bsai_cut_rev_5prime: str = 'G'
        bsai_cuts_check: Tuple[DSeq, DSeq] = (
            DSeq(bsai_cut_fwd_5prime,
                bsai_cut_rev_3prime,
                overhang=3
            ),
            DSeq(bsai_cut_fwd_3prime,
                bsai_cut_rev_5prime,
                overhang=-4)
        )
        ds_bsai: DSeq = DSeq(fwd, rev)
        bsai_cuts: List[DSeq] = ds_bsai.cut('BsaI')
        for ds, ds_check in zip(bsai_cuts, bsai_cuts_check):
            self.assertEqual(ds, ds_check)

        # Test double cutter BaeI
        seq_buf: str = 'CCCCCC'
        fwd = 'A'*10 + 'AC' + 'A'*4 + 'GTAYC' + 'A'*12
        rev_rev = 'T'*15 + 'TG' + 'T'*4 + 'CATRG' + 'T'*7

        fwd: str = fwd.replace('Y', 'C')
        rev_rev: str = rev_rev.replace('R', 'G')
        rev: str = reverse(rev_rev)
        ds_baei: DSeq = DSeq(seq_buf + fwd + seq_buf,
                            seq_buf + rev + seq_buf, overhang=5)
        baei_cuts: List[DSeq] = ds_baei.cut('BaeI')

        end_cut: DSeq = DSeq(seq_buf, seq_buf, overhang=5)
        center_cut: DSeq = DSeq(fwd, rev, overhang=5)

        self.assertEqual(len(baei_cuts), 3)
        self.assertEqual(baei_cuts[0], end_cut)
        self.assertEqual(baei_cuts[1], center_cut)
        self.assertEqual(baei_cuts[2], end_cut)

        # BsaI mutliple cuts
        fwd: str = 'GGTCTCGAATTCAAATTT'
        rev: str = 'GAATTCGAGACCAAATTT'
        bsai_ds2 = DSeq(fwd, rev)
        bsai_ds2_double = bsai_ds2 + bsai_ds2
        bsai_double_cuts = bsai_ds2_double.cut('BsaI')
        self.assertEqual(len(bsai_double_cuts), 3)
        bsai_ds2_checks = ( DSeq(              'GGTCTCG',
                                 reverse('TTTAAACCAGAGCTTAA'),
                                    overhang=6),
                            DSeq(   'AATTCAAATTTGGTCTCG',
                                reverse('GTTTAAACCAGAGCTTAA'),
                                overhang=-4),
                            DSeq('AATTCAAATTT',
                                     'G',
                                overhang=-4)
        )
        for x, y in zip(bsai_double_cuts, bsai_ds2_checks):
            self.assertEqual(x, y)
    # end def


if __name__ == '__main__':
    x = TestDSeq()
    x.test_single()
    x.test_add()
    x.test_cut()