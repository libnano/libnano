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
seqrepeat.pxd
~~~~~~~~~~~~~~~~~~

Cython header file for seqrepeat.pyx -- allows for cross-project Cython.

This header mainly serves to help wrap functions written in raw c
that are wrapped in seqrepeat.pyx such that those functions need not
be recompiled nor the associated headers included

'''
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
