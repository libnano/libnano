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
tests.test_prostr
~~~~~~~~~~~~~~~~~

'''
import unittest

from libnano.datasets import dataset_container
from libnano.prostr import (
    ReverseTranslator,
    WeightedReverseTranslator,
    dnaToAA,
    rnaToAA,
)


class TestProStr(unittest.TestCase):

    def setUp(self):
        pass

    def test_ReverseTranslator(self):
        aa_seq = ''.join(dataset_container.AA_TO_DNA.keys())
        # DNA random
        a = ReverseTranslator(aa_seq)
        self.assertIsNotNone(a)
        for _ in range(1000):
            rt = a.revTranslate()
            self.assertEqual(
                dnaToAA(rt),
                aa_seq,
            )

        # RNA random
        a = ReverseTranslator(aa_seq, na_type='RNA')
        self.assertIsNotNone(a)
        for _ in range(1000):
            self.assertEqual(
                rnaToAA(
                    a.revTranslate(),
                ),
                aa_seq,
            )

        # Weighted by organism codon freqs
        for organism in dataset_container.CODON_FREQ_DATASET.keys():
            # DNA, weighted
            a = WeightedReverseTranslator(aa_seq, organism)
            self.assertIsNotNone(a)
            for _ in range(1000):
                self.assertEqual(
                    dnaToAA(a.revTranslate()),
                    aa_seq,
                )
            # RNA, weighted
            a = WeightedReverseTranslator(
                aa_seq,
                organism,
                na_type='RNA',
            )
            self.assertIsNotNone(a)
            for _ in range(1000):
                self.assertEqual(
                    rnaToAA(a.revTranslate()),
                    aa_seq,
                )

    def test_dna_to_aa(self):
        self.assertEqual(
            dnaToAA('TTTTTCAGT'),
            'FFS',
        )
        with self.assertRaises(KeyError):
            dnaToAA('TTTTTCAGTA')  # seq % 3 != 0
        with self.assertRaises(KeyError):
            dnaToAA('TTTTTCAGX')  # Invalid character
        # Exhaustive (test translation of every codon)
        dna_seq = ''
        aa_seq = ''
        for aa, codons in dataset_container.AA_TO_DNA.items():
            dna_seq += ''.join(codons)
            aa_seq += aa * len(codons)
        self.assertEqual(
            dnaToAA(dna_seq),
            aa_seq,
        )

    def test_rna_to_aa(self):
        self.assertEqual(
            rnaToAA('UUUUUCAGU'),
            'FFS',
        )
        with self.assertRaises(KeyError):
            rnaToAA('UUUUUCAGUU')  # seq % 3 != 0
        with self.assertRaises(KeyError):
            rnaToAA('UUUUUCAGX')  # Invalid character
        # Exhaustive (test translation of every codon)
        rna_seq = ''
        aa_seq = ''
        for aa, codons in dataset_container.AA_TO_RNA.items():
            rna_seq += ''.join(codons)
            aa_seq += aa * len(codons)
        self.assertEqual(
            rnaToAA(rna_seq),
            aa_seq,
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
