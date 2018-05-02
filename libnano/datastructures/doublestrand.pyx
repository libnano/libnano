cimport numpy as cnp
from copy import copy, deepcopy
from libnano.seqstr import (
    reverseComplement,
    complement
)

cdef class DoubleSequence:
    """ forward strand index is always
    """
    cdef strand_fwd
    cdef strand_rev
    cdef int fwd_idx
    cdef int fwd_max
    cdef int fwd_len
    cdef int rev_idx
    cdef int rev_max
    cdef int rev_len

    cdef bool is_circular

    def __init__(self, strand_fwd=None, strand_rev=None, rev_idx=0, is_circular=False):
        cdef int overlap_low
        def int overlap_high = -1

        if strand_fwd is None and strand_rev is None:
            raise ValueError("can't be both None")
        if strand_rev is None:
            strand_rev = reverseComplement(strand_fwd)
        elif strand_fwd is None:
            strand_fwd = reverseComplement(strand_rev)
        self.strand_fwd = strand_fwd
        self.strand_rev = strand_rev
        self.is_circular = is_circular
        if rev_idx < 0:
            fwd_idx = rev_idx
            fwd_max = fwd_idx + len(strand_fwd)
            rev_idx = 0
            rev_max = len(strand_fwd)
        else:
            fwd_idx = 0
            fwd_max = len(strand_fwd)

            rev_max = rev_idx + len(strand_rev)

        if fwd_idx == rev_idx:
            overlap_low = 0
        elif fwd_idx > rev_idx:
            if fwd_idx < rev_max:
                overlap_low = fwd_idx
                if fwd_max <= rev_max:
                    overlap_high = fwd_max
                else:
                    overlap_high = rev_max
            # elif fwd_idx > rev_max:
            #     "no overlap"
            #     break
            elif fwd_idx == rev_max:
                overlap_low = rev_max
                overlap_high = rev_max
            # else:
            #     "no overlap"
            #     break
        elif rev_idx > fwd_idx:
            if rev_idx < fwd_max:
                overlap_low = rev_idx
                if rev_max <= fwd_max:
                    overlap_high = rev_max
                else:
                    overlap_high = fwd_max
            # elif rev_idx > fwd_max:
            #     "no overlap"
            #     break
            elif rev_idx == fwd_max:
                overlap_low = fwd_max
                overlap_high = fwd_max
            # else:
            #     "no overlap"
            #     break

        if overlap_high > 0 and
            strand_fwd[overlap_low:overlap_high] != complement(strand_rev[overlap_low:overlap_high]):
            print("no hybrid")
            return
    # end def

    def __add__(self, x, y):
        """ don't copy info
        """
        if not isinstance(x, DoubleSequence):
            raise ValueError("Can't concatenate non-DoubleSequence of type {}", type(x))
        if not isinstance(y, DoubleSequence):
            raise ValueError("Can't concatenate non-DoubleSequence of type {}", type(y))
        if x.is_circular or y.is_circular:
            raise ValueError("can't add a circular DoubleSequence x:{}, y:{}",
                                x.is_circular, y.is_circular)
        return DoubleSequence(x.strand_fwd + y.strand_fwd, x.strand_rev + y.strand_rev)
    # end def

    def tail(self):
        return fw

    def __copy__(self):
        """
        """
        ds_copy = type(self)(copy(self.strand_fwd),
                        copy(self.strand_rev),
                        is_circular=self.is_circular)
        return ds_copy
    # end def

# end cdef

