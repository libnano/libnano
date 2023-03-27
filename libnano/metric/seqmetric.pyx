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
libnano.metric.seqmetric
~~~~~~~~~~~~~~~~~~~~~~~~

Simple DNA sequence metrics
'''
# from cpython cimport array

from libc.stdlib cimport (
    calloc,
    free,
)

from libnano.helpers cimport c_util

# import array


__doc__ = 'Simple DNA sequence metrics'


cdef extern from 'sm_seqmetric.h':
    int sm_maxRuns(char* seq, int seq_length, int* ret_arr)
    int sm_gcContent(char* seq, int length)
# end cdef

def maxRuns(object seq_obj):

    cdef char*          seq
    cdef Py_ssize_t     length
    cdef int*           runs_arr = NULL
    cdef list           out = []
    seq = c_util.obj_to_cstr_len(seq_obj, &length)
    runs_arr = <int*> calloc(6, sizeof(int))
    if runs_arr == NULL:
        raise OSError()
    sm_maxRuns(seq, length, runs_arr)
    for i in range(6):
        out.append(runs_arr[i])
    free(runs_arr)
    return out
# end def

def gcContent(object seq_obj):

    cdef char*          seq
    cdef Py_ssize_t     length
    cdef float          gc_content

    seq = c_util.obj_to_cstr_len(seq_obj, &length)
    gc_content = sm_gcContent(seq, length)

    return gc_content
# end def
