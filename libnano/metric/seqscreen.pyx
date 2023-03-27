
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
libnano.metric.seqscreen
~~~~~~~~~~~~~~~~~~~~~~~~

Simple filtering of DNA sequences

'''
cimport c_util

__doc__ = "Simple filtering of DNA sequences"

cdef extern from "sf_seqscreen.h":
    int sf_containsRun(
        char* seq, int seq_length,
        int maxA, int maxT, int maxG,
        int maxC, int maxAT, int maxGC
    )

    int sf_gcWindow(
        char* seq, int length,
        int gc_min_percent, int gc_max_percent,
        int window_size, int* window_start, int* gc_count
    )

from typing import (
    Tuple,
    Union,
)


def containsRun(
        seq_obj: Union[str, bytes],
        maxA: int,
        maxT: int,
        maxG: int,
        maxC: int,
        maxAT: int,
        maxGC: int,
) -> bool:
    '''
    Check a DNA sequence for A/T/G/C/AT/GC runs of a max length

    Args:
        seq: DNA sequence
        maxA: Maximum allowable run of A
        maxT: Maximum allowable run of T
        maxG: Maximum allowable run of G
        maxC: Maximum allowable run of C
        maxAT: Maximum allowable run of AT
        maxGC: Maximum allowable run of GC

    Returns:
        True if contains a run, False otherwise
    '''
    cdef:
        int     rc
        char*   seq
        Py_ssize_t  length

    seq = c_util.obj_to_cstr_len(
        seq_obj,
        &length,
    )

    rc = sf_containsRun(
        seq, length,
        maxA, maxT,
        maxG, maxC,
        maxAT, maxGC,
    )

    return True if rc else False


def gcWindow(
        seq_obj : Union[str, bytes],
        gc_min_percent: int,
        gc_max_percent: int,
        window_size: int
) -> Tuple[bool, int, int]:
    '''Check a DNA sequence for GC content in a sliding window of fixed size

    Args:
        seq: DNA sequence
        min_gc_percent: minimum GC percent within window (integer from 0-99)
        max_gc_percent: maximum GC percent within window (integer from 1-100)
        window_size: size of the sliding window in bases

    Returns:
        Tuple of the form::

        <has window>, <window start index>, <gc count>

    Raises:
        ValueError: gc_min_percent must be within 0-100
        ValueError: gc_max_percent must be within 0-100
        ValueError: gc_min_percent cannot be greater than gc_max_percent
        ValueError: window_size is larger than sequence length
    '''

    cdef:
        int         rc
        int         window_start = 0
        int         gc_count = 0
        char*       seq
        Py_ssize_t  length

    seq = c_util.obj_to_cstr_len(
        seq_obj,
        &length,
    )

    if gc_min_percent < 0 or gc_min_percent > 100:
        raise ValueError('gc_min_percent must be within 0-100')

    if gc_max_percent < 0 or gc_max_percent > 100:
        ValueError('gc_max_percent must be within 0-100')

    if gc_min_percent > gc_max_percent:
        ValueError('gc_min_percent cannot be greater than gc_max_percent')

    rc = sf_gcWindow(
        seq, length,
        gc_min_percent,
        gc_max_percent, window_size,
        &window_start, &gc_count
    )

    if rc == -1:
        raise ValueError('window_size is larger than sequence length')

    ret = True if rc else False
    return ret, window_start, gc_count
