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
    dna_to_aa,
    rna_to_aa,
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
            rt = a.rev_translate()
            self.assertEqual(
                dna_to_aa(rt),
                aa_seq,
            )

        # RNA random
        a = ReverseTranslator(aa_seq, na_type='RNA')
        self.assertIsNotNone(a)
        for _ in range(1000):
            self.assertEqual(
                rna_to_aa(
                    a.rev_translate(),
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
                    dna_to_aa(a.rev_translate()),
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
                    rna_to_aa(a.rev_translate()),
                    aa_seq,
                )

    def test_dna_to_aa(self):
        self.assertEqual(
            dna_to_aa('TTTTTCAGT'),
            'FFS',
        )
        with self.assertRaises(KeyError):
            dna_to_aa('TTTTTCAGTA')  # seq % 3 != 0
        with self.assertRaises(KeyError):
            dna_to_aa('TTTTTCAGX')  # Invalid character
        # Exhaustive (test translation of every codon)
        dna_seq = ''
        aa_seq = ''
        for aa, codons in dataset_container.AA_TO_DNA.items():
            dna_seq += ''.join(codons)
            aa_seq += aa * len(codons)
        self.assertEqual(
            dna_to_aa(dna_seq),
            aa_seq,
        )

    def test_rna_to_aa(self):
        self.assertEqual(
            rna_to_aa('UUUUUCAGU'),
            'FFS',
        )
        with self.assertRaises(KeyError):
            rna_to_aa('UUUUUCAGUU')  # seq % 3 != 0
        with self.assertRaises(KeyError):
            rna_to_aa('UUUUUCAGX')  # Invalid character
        # Exhaustive (test translation of every codon)
        rna_seq = ''
        aa_seq = ''
        for aa, codons in dataset_container.AA_TO_RNA.items():
            rna_seq += ''.join(codons)
            aa_seq += aa * len(codons)
        self.assertEqual(
            rna_to_aa(rna_seq),
            aa_seq,
        )


if __name__ == '__main__':
    unittest.main(verbosity=2)
