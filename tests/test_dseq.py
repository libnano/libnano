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
tests.dseq_test
~~~~~~~~~~~~~~~

'''
import unittest
from typing import (
    List,
    Tuple,
)

from libnano.dseq import DSeq
from libnano.seqstr import reverse  # type: ignore


class TestDSeq(unittest.TestCase):

    def test_single(self):
        _ = DSeq('GGATCCAAA')

    def test_add(self):
        dseq = DSeq('GGATCCAAA')
        _ = dseq + dseq
        dseq = DSeq(
            'GGATCCAAAG',
            'TTTGGATCCC',
            overhang=1,
        )
        _ = dseq + dseq

    def test_circularize(self):
        dseq = DSeq('GGATCCAAA')
        self.assertTrue(
            dseq.isCircularizable(),
        )
        dseq_circle = dseq.circularize()
        self.assertFalse(
            dseq_circle.isCircularizable(),
        )

    def test_cut(self):
        # 1. Test single cutter BsaI
        fwd = 'GGTCTCGAATTCAAA'
        rev = 'GAATTCGAGACCAAA'

        bsai_cut_fwd_5prime = 'GGTCTCG'
        bsai_cut_fwd_3prime = 'AATTCAAA'

        bsai_cut_rev_3prime = 'AATTCGAGACCAAA'
        bsai_cut_rev_5prime = 'G'
        bsai_cuts_check: Tuple[DSeq, DSeq] = (
            DSeq(
                bsai_cut_fwd_5prime,
                bsai_cut_rev_3prime,
                overhang=3,
            ),
            DSeq(
                bsai_cut_fwd_3prime,
                bsai_cut_rev_5prime,
                overhang=-4,
            ),
        )
        ds_bsai_dseq = DSeq(fwd, rev)
        bsai_cuts_dseq_list = ds_bsai_dseq.cut('BsaI')
        for ds, ds_check in zip(bsai_cuts_dseq_list, bsai_cuts_check):
            self.assertEqual(
                ds,
                ds_check,
            )

        # 2. Test double cutter BaeI
        seq_buf = 'CCCCCC'
        fwd = 'A' * 10 + 'AC' + 'A' * 4 + 'GTAYC' + 'A' * 12
        rev_rev = 'T' * 15 + 'TG' + 'T' * 4 + 'CATRG' + 'T' * 7

        fwd = fwd.replace('Y', 'C')
        rev_rev = rev_rev.replace('R', 'G')
        rev = reverse(rev_rev)
        ds_baei_dseq = DSeq(
            seq_buf + fwd + seq_buf,
            seq_buf + rev + seq_buf, overhang=5,
        )
        baei_cuts: List[DSeq] = ds_baei_dseq.cut('BaeI')

        end_cut_dseq = DSeq(
            seq_buf,
            seq_buf,
            overhang=5,
        )
        center_cut_dseq = DSeq(
            fwd,
            rev,
            overhang=5,
        )

        self.assertEqual(
            len(baei_cuts),
            3,
        )
        self.assertEqual(
            baei_cuts[0],
            end_cut_dseq,
        )
        self.assertEqual(
            baei_cuts[1],
            center_cut_dseq,
        )
        self.assertEqual(
            baei_cuts[2],
            end_cut_dseq,
        )

        # 3. BsaI mutliple cuts
        fwd = 'GGTCTCGAATTCAAATTT'
        rev = 'GAATTCGAGACCAAATTT'
        bsai_ds2 = DSeq(fwd, rev)
        bsai_ds2_double = bsai_ds2 + bsai_ds2
        bsai_double_cuts = bsai_ds2_double.cut('BsaI')
        self.assertEqual(
            len(bsai_double_cuts),
            3,
        )
        bsai_ds2_checks = (
            DSeq(
                'GGTCTCG',
                reverse('TTTAAACCAGAGCTTAA'),
                overhang=6,
            ),
            DSeq(
                'AATTCAAATTTGGTCTCG',
                reverse('GTTTAAACCAGAGCTTAA'),
                overhang=-4,
            ),
            DSeq(
                'AATTCAAATTT',
                'G',
                overhang=-4,
            ),
        )
        for x, y in zip(bsai_double_cuts, bsai_ds2_checks):
            self.assertEqual(x, y)

        # 4. Test circular cutter BsaI
        fwd = 'TTTGGTCTCGAATTCAAA'
        rev = 'TTTGAATTCGAGACCAAA'

        bsai_cut_fwd = 'AATTCAAATTTGGTCTCG'
        bsai_cut_rev = reverse('GTTTAAACCAGAGCTTAA')

        bsai_circ_cut_checks_dseq_list = [
            DSeq(
                bsai_cut_fwd,
                bsai_cut_rev,
                overhang=-4,
            ),
        ]
        ds_bsai_dseq = DSeq(fwd, rev)
        ds_bsai_circular_dseq = ds_bsai_dseq.circularize()
        bsai_circ_cuts: List[DSeq] = ds_bsai_circular_dseq.cut(
            'BsaI',
        )
        for ds, ds_check in zip(
            bsai_circ_cuts,
            bsai_circ_cut_checks_dseq_list,
        ):
            self.assertEqual(ds, ds_check)
    # end def


if __name__ == '__main__':
    x = TestDSeq()
    x.test_single()
    x.test_add()
    x.test_cut()
