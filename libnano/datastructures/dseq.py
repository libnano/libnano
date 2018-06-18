# -*- coding: utf-8 -*-
import sys
from typing import (
    Tuple,
    List,
    Union
)
from os.path import (
    abspath,
    dirname
)
import math
import re

try:
    import libnano
except:
    LIBNANO_PATH = dirname(dirname(dirname(dirname(__file__))))
    sys.path = [LIBNANO_PATH] + sys.path

from libnano.datastructures.seqrecord.seqrecordbase import (
    SeqRecord,
    AlphaEnum,
    ALPHABETS
)
from libnano.seqstr import (
    reverseComplement,
    complement,
    reverse
)
from libnano.fileio.naformat import (
    align_complement,
    string_align_complement,
    Alignment,
    five_prime_type,
    three_prime_type,
    PRIME_ENUM_MAP
)
from libnano.search.restriction import (
    RestrictionSearcher,
    RestrictionMatch,
    StrandDirEnum
)

DSEQ_STR: str = '%s\r\n%s' if sys.platform == 'win32' else '%s\n%s'

class DSeq(object):
    def __init__(self,
        fwd: str,
        rev: str = None,
        overhang: int = None,
        alphabet: int = AlphaEnum.DNA):
        '''
        Args:
            fwd: ``fwd`` maps to reference in ``ssw-py``
            rev: ``rev`` maps to read in ``ssw-py``
            overhang: Use to force an alignment.  + overhang means the 5' end
                of the fwd strand is shifted in the 3' direction of the rev strand
                - overhang means the 5' end of the fwd strand is shifted in the
                3' direction of the reverse strand
            alphabet: whether this is DNA or RNA
        '''

        assert(alphabet in ALPHABETS)
        self.alphabet: int = alphabet
        self.fwd: str = fwd
        self.overhang: int = 0


        if rev is None:
            if overhang is not None:
                raise ValueError("overhang can't be defined for without a reverse strand")
            else:
                self.rev: str = reverseComplement(fwd)
        else:
            self.rev: str = rev


        max_idx_fwd: int = len(fwd) - 1
        max_idx_rev: int = len(self.rev) - 1
        the_length: int = max(max_idx_fwd, max_idx_rev) + 1 # default

        self.alignment: Alignment
        self.rc_rev: str
        if overhang is None:
            alignment, self.rc_rev = align_complement(fwd, self.rev)

            if max_idx_fwd > alignment.reference_end: # posiive overhang
                self.overhang = max_idx_fwd - alignment.reference_end
                the_length: int = max_idx_fwd + alignment.read_start + 1

            elif max_idx_rev > alignment.read_end: # negative overhang
                self.overhang = alignment.read_end - max_idx_fwd
                the_length: int = max_idx_rev + alignment.reference_start + 1
            self.alignment = alignment

        else:
            self.overhang = overhang
            if overhang < 0:
                reference_start: int = -overhang
                read_start: int = 0
                delta = min(max_idx_fwd + overhang, max_idx_rev)
                reference_end: int = reference_start + delta
                read_end: int = delta
                the_length: int = max_idx_rev + reference_start + 1
            elif overhang > 0:
                reference_start: int = 0
                read_start: int = overhang
                delta = min(max_idx_fwd, max_idx_rev - read_start)
                reference_end: int = delta
                read_end: int = read_start + delta
                the_length: int = max_idx_fwd + read_start + 1
            else:
                reference_start: int = 0
                read_start: int = 0
                delta: int = min(max_idx_fwd, max_idx_rev)
                reference_end: int = delta
                read_end: int = delta
                the_length: int = max(max_idx_fwd, max_idx_rev) + 1
            self.alignment = Alignment(
                '', # null
                0,  # null
                0,  # null
                reference_start,
                reference_end,
                read_start,
                read_end
            )
            self.rc_rev = reverseComplement(rev)
        self.the_length: int = the_length
    # end def

    def __str__(self) -> str:
        return self.highlight(do_highlight=False)

    def highlight(self, do_highlight=True):
        x, y = string_align_complement(
            self.fwd,
            self.rev,
            self.alignment,
            self.rc_rev,
            do_highlight=do_highlight)
        return DSEQ_STR % (x, y)

    def __add__(self, b: 'DSeq') -> 'DSeq':
        '''(1) Concatenates the forward strand with forward strand
        and the reverse strand with the reverse strand and preserves order
        (2) Realligns the two :class:`DSeq` objects involved

        Args:
            b: a :class:`DSeq` object

        Returns:
            a :class:`DSeq` object

        Raises:
            ValueError
        '''
        if isinstance(b, DSeq):
            fwd = self.fwd + b.fwd
            rev = self.rev + b.rev
            return DSeq(fwd, rev)
        else:
            raise ValueError("{} object not a DSeq".format(b))
    # end def

    def __len__(self) -> int:
        return self.the_length

    def __getitem__(self, x: Union[int, slice]) -> 'DSeq':
        fwd: str = self.fwd
        rc_rev: str = self.rc_rev
        overhang: int = self.overhang
        if not isinstance(x, slice):
            sl: slice = slice(x, x + 1, 1)
        else:
            sl: slice = x
        start_idx: int = sl.start or 0
        if overhang > 0:    # fwd shifted in 3' direction
            fwd_sl_stop: int = sl.stop - overhang if sl.stop is not None else None
            fwd_sl: slice = slice(start_idx - overhang, fwd_sl_stop, sl.step)
            fwd_out: str = fwd[fwd_sl]
            rev_out: str = reverseComplement(rc_rev[sl])
            overhang_out: int = max(overhang - start_idx, 0)
        else:   # fwd shited in the 5' direction relative to the reverse (negative overhang)
            rev_sl_stop: int = sl.stop + overhang if sl.stop is not None else None
            rev_sl: slice = slice(start_idx + overhang, rev_sl_stop, sl.step)
            rev_out: str = reverseComplement(rc_rev[rev_sl])
            fwd_out: str = fwd[sl]
            overhang_out: int = min(start_idx + overhang, 0)
        return DSeq(fwd_out,
                    rev_out,
                    overhang_out,
                    alphabet=self.alphabet)
    # end def

    def getReverseIdx(self, fwd_idx: int) -> int:
        '''Does not check for a valid index

        Args:
            fwd_idx: the index into the forward string

        Returns:
            the index into the reverse strand string.
        '''
        return fwd_idx + overhang
    # end def

    def getForwardIdxFrom3PrimeIdx(self, rev_idx: int) -> int:
        '''Does not check for a valid index

        Args:
            rev_idx: the index into the reverse string

        Returns:
            the index into the forward strand string.
        '''
        return rev_idx - self.overhang
    # end def

    def getForwardIdxFrom5PrimeIdx(self, rev_idx: int) -> int:
        '''Does not check for a valid index

        Args:
            rev_idx: the index into the reverse string

        Returns:
            the index into the forward strand string.
        '''
        len_rev: int = len(self.rev)
        out_rev_idx: int = len_rev - rev_idx # invert the index
        return out_rev_idx - self.overhang
    # end def


    def five_prime_end(self) -> Tuple[str, str]:
        '''Return what kind of end is overhanging the 5' end of the
        forward strand
        '''
        res, val  = five_prime_type(self.alignment, self.fwd, self.rev)
        return PRIME_ENUM_MAP.get(res), val
    # end def

    def three_prime_end(self) -> Tuple[str, str]:
        '''Return what kind of end is overhanging the 3' end of the forward
        strand
        '''
        res, val  = three_prime_type(self.alignment, self.fwd, self.rev)
        return PRIME_ENUM_MAP.get(res), val
    # end def

    # for m in regex_fwd_comp.finditer(seq)

    def cut(self,
            enzyme: str,
            restriction_searcher: RestrictionSearcher = None) -> List['DSeq']:
        if restriction_searcher is None:
            restriction_searcher = RestrictionSearcher(enzyme)
        rs: RestrictionSearcher = restriction_searcher
        fwd: str = self.fwd
        rev: str = self.rev
        reverse_rev: str = reverse(rev)
        fwd_match_list: List[RestrictionMatch] = rs.findSites(fwd)[0]
        # print(fwd_match_list)
        out: List[DSeq] = []
        if len(fwd_match_list) == 0:
            return out
        else:
            out = []
            for fwd_match in fwd_match_list:
                # if fwd_match.strand_dir == StrandDirEnum.Forward:
                #     print("A forward match!")
                # else:
                #     print("A reverse match!")
                rev_regex = fwd_match.pair_regex
                rev_match_list: List[Tuple[int, int]] = [(m.start(), m.end())for m in re.finditer(rev_regex, rev)]
                if len(rev_match_list) == 0:
                    return out
                else:
                    # print("rev strand has it!")
                    # print(fwd_match)
                    # print(rev_match_list)
                    rev_match_idxs = rev_match_list[0]
                    potential_5p_idx  = self.getForwardIdxFrom5PrimeIdx(rev_match_idxs[1])
                    assert potential_5p_idx == fwd_match.start_idx
                    overhang: int = self.overhang
                    rev_offset: int = overhang if overhang > 0 else 0

                    fwd_cuts: Tuple[Tuple[int, int]] = fwd_match.cut_idxs

                    # NOTE: THIS IS SETUP TO HANDLE ONLY SINGLE CUT SITE AT THE MOMENT
                    cut_delta: int =  fwd_cuts[0][1] - fwd_cuts[1][1]

                    fwd_slice_1 = slice(0,
                                        fwd_match.start_idx+fwd_cuts[0][1],
                                        1)
                    fwd_slice_2 = slice(fwd_match.start_idx+fwd_cuts[0][1],
                                        None,
                                        1)
                    fwd_out1, fwd_out2 = fwd[fwd_slice_1], fwd[fwd_slice_2]
                    # print(fwd_out1, fwd_out2)

                    rev_cut_end_idx: int = fwd_match.start_idx + fwd_cuts[1][1] + overhang
                    rev_slice_1 = slice(0,
                                        rev_cut_end_idx,
                                        1)
                    rev_slice_2 = slice(rev_cut_end_idx,
                                        None,
                                        1)
                    rev_rev_out1, rev_rev_out2 = reverse_rev[rev_slice_1], reverse_rev[rev_slice_2]
                    # print(rev_rev_out1, rev_rev_out2)

                    out.append((DSeq(   fwd_out1,
                                        reverse(rev_rev_out1)
                                ),
                                DSeq(   fwd_out2,
                                        reverse(rev_rev_out2),
                                        overhang=cut_delta
                                ) )
                    )
        return out
    # end def
# end class



if __name__ == '__main__':
    # a = DSeq("GGATCCAAA", "TTTGGATC")
    # print(a + a)

    from pydna.dseq import Dseq
    def tester(i, fwd, rev, overhang=None):
        print("Test#%d:" % (i))
        a = Dseq(fwd, rev, ovhg=overhang)
        print(a.fig())
        print(a.five_prime_end())
        print(a.three_prime_end())

        a = DSeq(fwd, rev, overhang=overhang)
        print(a.alignment.reference_start, a.alignment.reference_end)
        print(a.highlight())
        print(a.five_prime_end())
        print(a.three_prime_end())
        print(a.alignment.read_start, a.alignment.read_end)
        print('')
        return i + 1
    # end def

    i = 0
    # i = tester(i, "AAA", "TTT")
    # i = tester(i, "AAA", "TTT", overhang=1)
    # i = tester(i,"AAA", "TTT", overhang=-1)
    # i = tester(i,"AAAG", "TTT", overhang=-1)
    # i = tester(i, "AAAGCCC", "TTT", overhang=-2)
    # i = tester(i, "TTT", "AAAGCCC", overhang=4)

    'EcoRI'
    fwd = 'GGTCTCGAATTCAAA'
    rev = 'GAATTCGAGACCAAA'
    a = DSeq(fwd, rev)
    print(a, len(a))
    b = a.cut('BsaI')
    print(a)
    print("The cut")
    for x in b[0]:
        print(x)

    b = a.cut('EcoRI')
    print(a)
    print("The cut")
    for x in b[0]:
        print(x)
    # print("SLICE")
    # print(a[5:])



    # print("$$$$$")
    # pyb = Dseq(fwd, rev)
    # print(len(pyb))
    # print(pyb.fig())
    # # print(pyb[2:].fig(), len(pyb[2:]))
    # print(pyb[5:].fig(), len(pyb[5:]))