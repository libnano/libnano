# -*- coding: utf-8 -*-
import sys
from typing import Tuple
from os.path import abspath, dirname
import math

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
        self.alignment: Alignment
        self.rc_rev: str

        if overhang is None:
            self.alignment, self.rc_rev = align_complement(fwd, self.rev)

            if max_idx_fwd > alignment.reference_end:
                self.overhang = max_idx_fwd - alignment.reference_end
            elif max_idx_rev > alignment.read_end:
                self.overhang = alignment.read_end - max_idx_fwd
        else:
            self.overhang = overhang
            if overhang < 0:
                reference_start: int = -overhang
                read_start: int = 0
                delta = min(max_idx_fwd + overhang, max_idx_rev)
                reference_end: int = reference_start + delta
                read_end: int = delta
            elif overhang > 0:
                reference_start: int = 0
                read_start: int = overhang
                delta = min(max_idx_fwd, max_idx_rev - read_start)
                reference_end: int = delta
                read_end: int = read_start + delta
            else:
                reference_start: int = 0
                read_start: int = 0
                delta: int = min(max_idx_fwd, max_idx_rev)
                reference_end: int = delta
                read_end: int = delta
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

# end class



if __name__ == '__main__':
    a = DSeq("GGATCCAAA", "TTTGGATC")
    print(a + a)

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
    i = tester(i, "AAA", "TTT")
    i = tester(i, "AAA", "TTT", overhang=1)
    i = tester(i,"AAA", "TTT", overhang=-1)
    i = tester(i,"AAAG", "TTT", overhang=-1)
    i = tester(i, "AAAGCCC", "TTT", overhang=-2)
    i = tester(i, "TTT", "AAAGCCC", overhang=4)