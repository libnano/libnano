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
libnano.dseq
~~~~~~~~~~~~
'''
from __future__ import annotations

import copy
import os.path as op
import re
import sys
from typing import (
    List,
    Optional,
    Tuple,
    Union,
)

try:
    import libnano
except (ImportError, ModuleNotFoundError):
    LIBNANO_PATH = op.dirname(
        op.dirname(__file__),
    )
    sys.path.append(LIBNANO_PATH)


from libnano.fileio.naformat import (
    PRIME_ENUM_MAP,
    Alignment,
    PrimeEnum,
    align_complement,
    five_prime_type,
    string_align_complement,
    three_prime_type,
)
from libnano.search.restriction import (  # type: ignore
    RestrictionMatch,
    RestrictionSearcher,
)
from libnano.seqrecord.seqrecordbase import (  # type: ignore
    ALPHABETS,
    AlphaEnum,
)
from libnano.seqstr import (  # type: ignore
    reverse,
    reverseComplement,
)

NEWLINE_STR: str = '\r\n' if sys.platform == 'win32' else '\n'
DSEQ_STR: str = (
    '%s(%d) is_circular: %s' + NEWLINE_STR +
    '%s' + NEWLINE_STR +
    '%s'
)


class DSeq:
    '''
    '''

    def __init__(
        self,
        fwd: str,
        rev: str = '',
        overhang: Optional[int] = None,
        is_circular: bool = False,
        alphabet: int = AlphaEnum.DNA,
    ):
        '''
        Args:
            fwd: ``fwd`` maps to reference in ``ssw-py``
            rev: ``rev`` maps to read in ``ssw-py``
            overhang: Use to force an alignment.  + overhang means the 5' end
                of the fwd strand is shifted in the 3' direction of the rev
                strand - overhang means the 5' end of the fwd strand is shifted
                in the 3' direction of the reverse strand
            alphabet: whether this is DNA or RNA
        '''

        assert (alphabet in ALPHABETS)
        self.alphabet: int = alphabet
        self.fwd_str = fwd
        self.overhang: int = 0

        if not rev:
            if overhang is not None:
                raise ValueError(
                    "overhang can't be defined for without a reverse strand",
                )
            else:
                self.rev_str = reverseComplement(fwd)
        else:
            self.rev_str = rev

        max_idx_fwd = len(fwd) - 1
        max_idx_rev = len(self.rev_str) - 1
        the_length = max(
            max_idx_fwd,
            max_idx_rev,
        ) + 1  # default

        self.alignment: Alignment
        self.rc_rev_str = ''

        if overhang is None:
            alignment, self.rc_rev_str = align_complement(
                fwd,
                self.rev_str,
            )

            if max_idx_fwd > alignment.reference_end:  # positive overhang
                self.overhang = (
                    max_idx_fwd - alignment.reference_end
                )
                the_length = (
                    max_idx_fwd + alignment.read_start + 1
                )
            elif max_idx_rev > alignment.read_end:  # negative overhang
                self.overhang = (
                    alignment.read_end - max_idx_fwd
                )
                the_length = (
                    max_idx_rev + alignment.reference_start + 1
                )
            self.alignment = alignment

        else:
            self.overhang = overhang
            if overhang < 0:
                reference_start_idx = -overhang
                read_start_idx = 0
                delta = min(
                    max_idx_fwd + overhang,
                    max_idx_rev,
                )
                reference_end_idx = reference_start_idx + delta
                read_end_idx = delta
                the_length = max_idx_rev + reference_start_idx + 1
            elif overhang > 0:
                reference_start_idx = 0
                read_start_idx = overhang
                delta = min(
                    max_idx_fwd,
                    max_idx_rev - read_start_idx,
                )
                reference_end_idx = delta
                read_end_idx = read_start_idx + delta
                the_length = max_idx_fwd + read_start_idx + 1
            else:
                reference_start_idx = 0
                read_start_idx = 0
                delta = min(
                    max_idx_fwd,
                    max_idx_rev,
                )
                reference_end_idx = delta
                read_end_idx = delta
                the_length = max(
                    max_idx_fwd,
                    max_idx_rev,
                ) + 1
            self.alignment = Alignment(
                '',  # null
                0,  # null
                0,  # null
                reference_start_idx,
                reference_end_idx,
                read_start_idx,
                read_end_idx,
            )
            self.rc_rev_str = reverseComplement(rev)
        self.the_length = the_length
        if is_circular:
            type3, _ = self.three_prime_end()
            type5, _ = self.five_prime_end()
            if ((type3 != PrimeEnum.BLUNT) and
                    (type5 != PrimeEnum.BLUNT)):
                raise ValueError('DNA is_circular but ends cannot mate')
        self.is_circular = is_circular

    def __str__(self) -> str:
        return self.highlight(do_highlight=False)

    def highlight(self, do_highlight: bool = True) -> str:
        '''Highlight the string with

        Args:
            do_highlight: Default is :const:`True`.
        '''
        x, y = string_align_complement(
            self.fwd_str,
            self.rev_str,
            self.alignment,
            self.rc_rev_str,
            do_highlight=do_highlight,
        )
        return DSEQ_STR % (
            type(self).__name__,
            len(self),
            self.is_circular,
            x, y,
        )

    def __add__(self, b: DSeq) -> DSeq:
        '''(1) Concatenates the forward strand with forward strand
        and the reverse strand with the reverse strand and preserves order
        (2) Realligns the two :class:`DSeq` objects involved.

        Args:
            b: a :class:`DSeq` object

        Returns:
            a :class:`DSeq` object

        Raises:
            ValueError, TypeError
        '''
        if isinstance(b, DSeq):
            if self.is_circular or b.is_circular:
                raise TypeError(
                    f'Cannot concatenate circular DSeq: {self} + {b}',
                )

            type3, seq3 = self.three_prime_end()
            type5, seq5 = b.five_prime_end()
            if type3 == type5 and len(seq3) == len(seq5):
                if seq3 != reverseComplement(seq5):
                    raise TypeError('Ends not complimentary')
                fwd = self.fwd_str + b.fwd_str
                rev = self.rev_str + b.rev_str
                return DSeq(
                    fwd,
                    rev,
                    self.overhang,
                )
            else:
                raise TypeError('Ends not compatible')
        else:
            raise ValueError(f'{b} object not a DSeq')

    def __len__(self) -> int:
        return self.the_length

    def __getitem__(self, x: Union[int, slice]) -> 'DSeq':
        fwd: str = self.fwd_str
        rc_rev: str = self.rc_rev_str
        overhang: int = self.overhang

        if not isinstance(x, slice):
            slc = slice(x, x + 1, 1)
        else:
            slc = x
        start_idx: int = slc.start or 0
        if overhang > 0:    # fwd shifted in 3' direction
            if slc.stop is not None:
                fwd_sl_stop_idx = (slc.stop - overhang)
            else:
                fwd_sl_stop_idx = None

            fwd_slc = slice(
                start_idx - overhang,
                fwd_sl_stop_idx,
                slc.step,
            )
            fwd_out_str: str = fwd[fwd_slc]
            rev_out_str: str = reverseComplement(rc_rev[slc])
            overhang_out: int = max(
                overhang - start_idx,
                0,
            )
        # fwd shited in the 5' direction relative to the reverse
        # (negative overhang)
        else:
            if slc.stop is not None:
                rev_slc_stop = slc.stop + overhang
            else:
                rev_slc_stop = None

            rev_slc: slice = slice(
                start_idx + overhang,
                rev_slc_stop,
                slc.step,
            )
            rev_out_str = reverseComplement(
                rc_rev[rev_slc],
            )
            fwd_out_str = fwd[slc]
            overhang_out = min(
                start_idx + overhang,
                0,
            )
        return DSeq(
            fwd_out_str,
            rev_out_str,
            overhang_out,
            alphabet=self.alphabet,
        )

    def __copy__(self) -> DSeq:
        '''No need to copy strings in python because they are immutable
        '''
        return type(self)(
            self.fwd_str,
            self.rev_str,
            self.overhang,
            self.is_circular,
            self.alphabet,
        )

    def __eq__(self, other) -> bool:
        elements = (
            'the_length',
            'fwd',
            'rev',
            'overhang',
            'is_circular',
        )
        for x in elements:
            if getattr(other, x) != getattr(self, x):
                return False
        return True

    def __iter__(self):
        '''NOTE: This is required to be defined such the __getitem__ is not
        called as a fall back for iteration
        '''
        for i in range(self.the_length):
            yield self.__getitem__(i)
        return

    def repeat(self, count: int) -> DSeq:
        '''Generatre

        Args:
            count: must be 2 or greater

        Returns:
            repeated concatenation :class:`DSeq`
        '''
        if not isinstance(count, int):
            raise TypeError(f'Must repeat integral times not {type(count)}')
        if count < 2:
            raise ValueError(f'Count must be 2 or greater not {count}')
        x = self
        for i in range(count - 1):
            x += self
        return x

    def circularize(self) -> DSeq:
        '''
        Returns:
            new :class:`DSeq` either a copy if already circular or a new object
        '''
        if self.is_circular:
            return self.__copy__()
        if self.isCircularizable():
            return DSeq(
                self.fwd_str,
                self.rev_str[-self.overhang:] + self.rev_str[:-self.overhang],
                0,
                is_circular=True,
                alphabet=self.alphabet,
            )
        else:
            raise TypeError("This object can't be circularized")

    def isCircularizable(self) -> bool:
        '''
        Returns:
            False if already circular or if ends don't match
        '''
        if self.is_circular:
            return False
        else:
            type3, seq3 = self.three_prime_end()
            type5, seq5 = self.five_prime_end()
            return type3 == type5 and seq3 == reverseComplement(seq5)

    def toLinear(self) -> DSeq:
        '''NOTE: you are not able to spcify the break point.  Uses internal
        string beginning and end

        Returns:
            new :class:`DSeq` either a copy if already linear or a new object
        '''
        x = self.__copy__()
        x.is_circular = False
        return x

    def getReverseIdx(self, fwd_idx: int) -> int:
        '''Does not check for a valid index

        Args:
            fwd_idx: the index into the forward string

        Returns:
            the index into the reverse strand string.
        '''
        return fwd_idx + self.overhang

    def getForwardIdxFrom3PrimeIdx(self, rev_idx: int) -> int:
        '''Does not check for a valid index

        Args:
            rev_idx: the index into the reverse string

        Returns:
            the index into the forward strand string.
        '''
        return rev_idx - self.overhang

    def getForwardIdxFrom5PrimeIdx(self, rev_idx: int) -> int:
        '''Does not check for a valid index

        Args:
            rev_idx: the index into the reverse string

        Returns:
            the index into the forward strand string.
        '''
        len_rev: int = len(self.rev_str)
        out_rev_idx: int = len_rev - rev_idx  # invert the index
        return out_rev_idx - self.overhang

    def five_prime_end(self) -> Tuple[int, str]:
        '''
        Returns:
            what kind of end is overhanging the 5' end of the forward
            strand of form::

                PrimeEnum.<X>, <sequence string>
        '''
        res, val = five_prime_type(self.alignment, self.fwd_str, self.rev_str)
        return res, val

    def three_prime_end(self) -> Tuple[int, str]:
        '''
        Returns:
            what kind of end is overhanging the 3' end of the forward
            strand of form::

                PrimeEnum.<X>, <sequence string>
        '''
        res, val = three_prime_type(
            self.alignment,
            self.fwd_str,
            self.rev_str,
        )
        return res, val

    @staticmethod
    def _cut(
        dseq: DSeq,
        rs: RestrictionSearcher,
    ) -> List[DSeq]:
        '''Helper for :method:`cut`
        '''
        out: List['DSeq'] = []
        fwd: str = dseq.fwd_str
        rev: str = dseq.rev_str
        reverse_rev: str = reverse(rev)
        fwd_match_list: List[RestrictionMatch] = rs.findSites(fwd)[0]
        if len(fwd_match_list) == 0:
            return [dseq]
        else:
            for fwd_match in fwd_match_list:
                rev_regex = fwd_match.pair_regex
                rev_match_list: List[Tuple[int, int]] = [
                    (m.start(), m.end())
                    for m in re.finditer(rev_regex, rev)
                ]
                if len(rev_match_list) == 0:
                    continue
                else:
                    rev_match_idxs = rev_match_list[-1]
                    potential_5p_idx = dseq.getForwardIdxFrom5PrimeIdx(
                        rev_match_idxs[1],
                    )
                    assert (potential_5p_idx == fwd_match.start_idx)

                    fwd_cuts: Tuple[Tuple[int, int], ...] = fwd_match.cut_idxs

                    # traverse through the cutsites in reverse
                    fwd = copy.copy(dseq.fwd_str)
                    reverse_rev = reverse(dseq.rev_str)
                    self_overhang = dseq.overhang

                    for fwd_cut, rev_cut in zip(
                        fwd_cuts[0][::-1],
                        fwd_cuts[1][::-1],
                    ):

                        cut_overhang = fwd_cut - rev_cut

                        fwd_cut_end_idx: int = (
                            fwd_match.start_idx + fwd_cut
                        )
                        fwd_slice_5p = slice(
                            0,
                            fwd_cut_end_idx,
                            1,
                        )
                        fwd_slice_3p = slice(
                            fwd_cut_end_idx,
                            None,
                            1,
                        )

                        rev_cut_end_idx = (
                            fwd_match.start_idx + rev_cut + self_overhang
                        )
                        rev_slc_3p = slice(
                            0,
                            rev_cut_end_idx,
                        )
                        rev_slc_5p = slice(
                            rev_cut_end_idx,
                            None,
                            1,
                        )

                        out.append(
                            DSeq(
                                fwd[fwd_slice_3p],
                                reverse(reverse_rev[rev_slc_5p]),
                                overhang=cut_overhang,
                            ),
                        )

                        fwd = fwd[fwd_slice_5p]
                        reverse_rev = reverse_rev[rev_slc_3p]
                    # end for
                    out.append(
                        DSeq(
                            fwd,
                            reverse(reverse_rev),
                            dseq.overhang,
                        ),
                    )
                '''NOTE: break out and then check the new DSeq's for additional
                cutsites'''
                break
            # end for
        # recurse at the 3prime end of the forward strand
        len_out = len(out)
        if len_out > 1:
            return (
                out[len_out - 1:0:-1] +
                dseq._cut(out[0], rs)
            )
        else:
            return out

    def cut(
        self,
        enzyme: str,
        restriction_searcher: Optional[RestrictionSearcher] = None,
    ) -> List[DSeq]:
        '''Return a list of the cut :class:`DSeq` objects created from cutting
        a Dseq with the ``enzyme`` argument.  If the enzyme finds more than one
        cutsite the list will be longer.  Cuts in the 5' to the 3' direction
        of the forward strand

        Args:
            enzyme: enzyme to check to cut
            restriction_searcher: default is None.  Creates one if needed

        Returns:
            List of :class:`DSeq` cuts from this :class:`DSeq` object. The list
            will be empty if no cuts are found
        '''
        if restriction_searcher is None:
            restriction_searcher = RestrictionSearcher(enzyme)
        is_circular: bool = self.is_circular
        dseq: DSeq
        if is_circular:
            dseq = self.toLinear()
            dseq = dseq.repeat(3)
        else:
            dseq = self
        out_list: List[DSeq] = self._cut(dseq, restriction_searcher)
        if is_circular:
            len_fwd: int = len(self.fwd_str)
            len_total: int = 0
            filtered_out_list: List[DSeq] = []
            '''Since the 3 repeat has an arbitrary cut we
            need to go from 1 to the -1 entry in the list
            of items
            '''
            for item in out_list[1:-1]:
                filtered_out_list.append(item)
                len_total += len(item.fwd_str)
                if len_total >= len_fwd:
                    break
            return filtered_out_list
        else:
            return out_list
# end class


if __name__ == '__main__':
    pass
    # a = DSeq("GGATCCAAA", "TTTGGATC")
    # print(a + a)

    # from pydna.dseq import Dseq

    # def tester(i, fwd, rev, overhang=None):
    #     print('Test#%d:' % (i))
    #     a = Dseq(fwd, rev, ovhg=overhang)
    #     print(a.fig())
    #     print(a.five_prime_end())
    #     print(a.three_prime_end())

    #     a = DSeq(fwd, rev, overhang=overhang)
    #     print(a.alignment.reference_start, a.alignment.reference_end)
    #     print(a.highlight())
    #     print(a.five_prime_end())
    #     print(a.three_prime_end())
    #     print(a.alignment.read_start, a.alignment.read_end)
    #     print('')
    #     return i + 1
    # # end def

    # i = 0
    # # i = tester(i, "AAA", "TTT")
    # # i = tester(i, "AAA", "TTT", overhang=1)
    # # i = tester(i,"AAA", "TTT", overhang=-1)
    # # i = tester(i,"AAAG", "TTT", overhang=-1)
    # # i = tester(i, "AAAGCCC", "TTT", overhang=-2)
    # # i = tester(i, "TTT", "AAAGCCC", overhang=4)

    # fwd = 'GGTCTCGAATTCAAA'
    # rev = 'GAATTCGAGACCAAA'
    # BsaI_ds = DSeq(fwd, rev)
    # print(BsaI_ds, len(BsaI_ds))
    # b = BsaI_ds.cut('BsaI')
    # print(BsaI_ds)
    # print('The cut')
    # print(b)
    # for x in b:
    #     print(x)

    # # b = a.cut('EcoRI')
    # # print(a)
    # # print("The cut")
    # # for x in b[0]:
    # #     print(x)
    # # print("SLICE")
    # # print(a[5:])

    # from libnano.datasets.build_enzyme_dataset import REGEX_BASE_LUT
    # BUF = 'CCCCCC'
    # fwd = BUF + 'A'*10 + 'AC' + 'A'*4 + 'GTAYC' + 'A'*12 + BUF
    # rev_rev = BUF + 'T'*15 + 'TG' + 'T'*4 + 'CATRG' + 'T'*7 + BUF

    # fwd = fwd.replace('Y', 'C')
    # rev_rev = rev_rev.replace('R', 'G')
    # rev = reverse(rev_rev)

    # baei_ds = DSeq(fwd, rev, overhang=5)
    # print('BaeI')
    # print(baei_ds)
    # b = baei_ds.cut('BaeI')
    # print("The cut")
    # for x in b:
    #     print(x)
    # print("$$$$$")
    # pyb = Dseq(fwd, rev)
    # print(len(pyb))
    # print(pyb.fig())
    # # print(pyb[2:].fig(), len(pyb[2:]))
    # print(pyb[5:].fig(), len(pyb[5:]))
