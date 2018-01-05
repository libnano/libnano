"""
seqrepeat.pxd
~~~~~~~~~~~~~~~~~~

Cython header file for seqrepeat.pyx -- allows for cross-project Cython.

This header mainly serves to help wrap functions written in raw c
that are wrapped in seqrepeat.pyx such that those functions need not
be recompiled nor the associated headers included
"""
from libc.stdint cimport uint64_t

cdef extern from "sr_seqrepeat.h":
    ctypedef struct repeatroot_t:
        int count
        int start
        int last

    ctypedef struct repeatcheck_t:
        repeatroot_t *repeatroot_arr
        int word_size
        int seq_length
        int count_max
        int *repeat_idx_arr

#cdef class RepeatCheck:
#    cdef inline repeatcheck_t* data(RepeatCheck self)