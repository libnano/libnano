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
libnano.metric.seqrepeat
~~~~~~~~~~~~~~~~~~~~~~~~

Find repeats in a DNA sequence

'''
from libnano.helpers cimport c_util

import numpy as np

cimport numpy as cnp
from cpython.int cimport PyInt_FromLong
from cpython.mem cimport (
    PyMem_Free,
    PyMem_Malloc,
)
from cpython.ref cimport PyObject
from libc.string cimport memset

__doc__ = 'Find repeats in a DNA sequence'

from typing import (
    Dict,
    List,
    Tuple,
    Union,
)


cdef extern from "sr_seqrepeat.h":
    int sr_buildRepeatData(
        char* seq,
        int seq_length,
        int word_size,
        repeatcheck_t** rcheck,
    )
    int sr_freeRepeatCheck(repeatcheck_t* rcheck)
    int sr_queryRepeatIdx(
        repeatcheck_t* rcheck,
        char* seq,
        char* query,
        int query_length,
    )
    int sr_repeatWindow(
        repeatcheck_t* rcheck,
        int window_size,
        int max_repeat_count,
        int** repeat_violation_idxs,
        int *violation_count,
        int** repeat_counts
    )


cdef extern from "Python.h":
    # we need this instead of from cpython.tuple cimport PyTuple_New,
    # PyTuple_SET_ITEM because object refcounting is broken with these guys
    PyObject* PyTuple_New(Py_ssize_t len)
    PyObject* PyLong_FromLong(long ival)
    int PyTuple_SET_ITEM(PyObject *p, Py_ssize_t pos, PyObject *o)


cdef class RepeatCheck:
    '''
    '''
    cdef:
        repeatcheck_t *rcheck
        bint is_extended
        seq_obj
        char* seq

    def __cinit__(
            self,
            seq_obj: Union[str, bytes],
            word_size: int,
    ):
        '''

        Args:
            seq_obj (str): string sequence
            word_size (int): the smallest considered repeat length
        '''
        cdef:
            Py_ssize_t seq_length
            char* seq = c_util.obj_to_cstr_len(
                seq_obj,
                &seq_length,
            )
            repeatcheck_t *rck = NULL

        self.seq_obj = seq_obj
        self.seq = seq

        if word_size > 32:
            raise ValueError(
                f'word_size {word_size} must be less than or equal to 32'
            )
        if word_size > seq_length:
            raise ValueError(
                f'word_size {word_size} must be less than or equal to '
                f'seq_length {seq_length}'
            )
        rck = <repeatcheck_t*> PyMem_Malloc(
            sizeof(repeatcheck_t)
        )
        if rck == NULL:
            raise OSError('RepeatCheck_init: problem allocating repeatcheck_t')
        else:
            self.rcheck = rck
            rck.repeatroot_arr = NULL
            rck.repeat_idx_arr = NULL
            if sr_buildRepeatData(
                    seq,
                    <int> seq_length,
                    word_size,
                    &rck,
            ) < 0:
                PyMem_Free(rck)
                self.rcheck = NULL
                raise OSError('RepeatCheck_init: problem building')

    def __dealloc__(RepeatCheck self):
        if self.rcheck != NULL:
            sr_freeRepeatCheck(self.rcheck)
            PyMem_Free(self.rcheck)
            self.rcheck = NULL

    cdef inline repeatcheck_t* data(RepeatCheck self):
        return self.rcheck

    def word_size(self):
        return self.rcheck.word_size

    def indices_of(
            self,
            query_obj,
    ) -> List[int]:
        ''' return all indices of a given query string
        root_word_size <= len(query_obj) <= max_word_size

        useful for validating the RepeatCheck code

        Returns:
            object (list) of indices
        '''
        cdef:
            Py_ssize_t query_length
            char* query = c_util.obj_to_cstr_len(
                query_obj,
                &query_length,
            )
            repeatcheck_t* rcheck = self.rcheck
            int idx = sr_queryRepeatIdx(
                rcheck,
                self.seq,
                query,
                query_length,
            )
            int* arr_query = rcheck.repeat_idx_arr
            list out

        if idx > -1:    # is there at least one hit?
            out = [idx]
            idx = arr_query[idx]
            while idx > 0:
                if self.seq_obj[idx:idx+query_length] == query_obj:
                    out.append(idx)
                idx = arr_query[idx]
            return out
        else:
            return []

    def screen_ends(
            self,
            num_bases: int,
    ) -> Tuple[bool, bool]:
        '''Screen the ends of the sequence for repeats on the order of
        self.word_size or greater

        Args:
            num_bases (int): Number of bases from either end to look for repeats

        Returns:
            Tuple[bool, bool]: (True/False, True/False) whether the 5' or 3'
            end have repeats

        Raises:
            ValueError: num_bases must be less than sequence length
        '''
        cdef:
            repeatcheck_t* rcheck = self.rcheck
            int* arr_query = rcheck.repeat_idx_arr
            Py_ssize_t seq_length = len(self.seq_obj)
            Py_ssize_t i, j

        is_5prime, is_3prime = False, False
        if num_bases > seq_length:
            raise ValueError(
                f'num_bases {num_bases} must be less than sequence length '
                f'{seq_length}'
            )
        j = seq_length - 1
        for i in range(num_bases):
            if arr_query[i] != 0:
                is_5prime = True
            if arr_query[j - i] < 0:
                is_3prime = True
        return (is_5prime, is_3prime)

    def gap_check(
            self,
            gap_length: int,
    ) -> List[Tuple[int, int, int]]:
        '''Finds all repeat pairs smaller than gap_length

        Args:
            gap_length (int): length that gap must be greater than or equal to
                to pass.  less than this length and it will be flagged

        Returns:
            List[Tuple[int, int, int]]: List of tuple of index pairs and the
                length of the repeat at those indices of the form:::

                    [(i, j, repeat length), ...]
        '''
        cdef:
            int i, j, k, next_idx
            repeatcheck_t* rcheck = self.rcheck
            int word_size = rcheck.word_size
            int* query_arr = rcheck.repeat_idx_arr
            list out = []
            int range_lim = rcheck.seq_length

        i = 0
        while i < range_lim:
            j = 0
            next_idx = query_arr[i]
            if next_idx > 0:
                if (next_idx - i) < gap_length:
                    k = i
                    for k in range(i, range_lim):
                        if query_arr[k] - next_idx != j:
                            break
                        j += 1
                    out.append((i, next_idx, word_size + (j - 1)))
                    i = k
                    continue
            i += 1
        return out

    def all_repeats(self) -> Dict[str, List[int]]:
        '''For each root repeat, return a dictionary of indices
        of the all repeats for which the root is the 5' end.  Filters
        out overlapping repeats returning only the list of the longest repeat
        string

        Returns:
            Dict: Key is the repeat string, value is a list of indices

                {key: [i, j, ...],
                        ...
                }
        '''
        cdef:
            int z, i, j, k, next_idx
            repeatcheck_t* rcheck = self.rcheck
            int word_size = rcheck.word_size
            int* query_arr = rcheck.repeat_idx_arr
            char* seqc = self.seq
            object seq_obj = self.seq_obj

            int range_lim = rcheck.seq_length
            dict res0 = {}

        res = {}
        ''' 1. Get all of the the base size index lists
        '''
        for i in range(range_lim):
            next_idx = query_arr[i]
            if next_idx > 0:
                key = seq_obj[i:i + word_size]
                if key not in res0:
                    res0[key] = [i, next_idx]
                else:
                    res0[key].append(next_idx)

        # 2. capture all extended repeats pairwise within
        # a word_size root and store in an output

        base = list(res0.items())
        for base_key, indices in base:
            for z, i in enumerate(indices[:-1]):
                for j in indices[z + 1:]:
                    k = word_size
                    for k in range(word_size, range_lim):
                        if seqc[i + k] != seqc[j + k]:
                            break
                    key = seq_obj[i:i + k] if k > word_size else base_key
                    if key not in res:
                        res[key] = [i, j]
                    else:
                        res[key].append(j)

        # pair things down only delete lists/sets of indices whose
        #  a. keys vary in length by 1
        #  b. the size of the list is the same
        #  c. the pairwise difference between all members is 1

        base = list(res.items())
        delete_list = []
        for i, keyval in enumerate(base[:-1]):
            key, indices = keyval
            for key_test, indices_test in base[i + 1:]:
                if len(indices) != len(indices_test) or abs(len(key) - len(key_test)) != 1:
                    continue
                if indices[0] > indices_test[0] and len(key_test) > len(key):
                    if any(x - y != 1 for x, y in zip(indices, indices_test)):
                        continue
                    else:
                        #print(key, key_test)
                        delete_list.append(key)
                elif indices[0] < indices_test[0] and len(key) > len(key_test):
                    if any(y - x != 1 for x, y in zip(indices, indices_test)):
                        continue
                    else:
                        #print(key, key_test)
                        delete_list.append(key_test)

        for x in delete_list:
            #print("deleting", x)
            del res[x]
        return res

    def window_check(
            self,
            window_size: int,
            max_repeat_count: int,
            get_counts: bool = True,
    ) -> Tuple[List, List]:
        '''
        Args:
            window_size: the length of the window.  Must be less than
                or equal to the sequence length
            max_repeat_count: max repeats of one kind allowed in a window
                before flagged into violation bin.  Example: 0 means no more than
                1 instance of a word can occur within a window
            get_counts: default True, whether to return the
                per window start index counts of repeats

        Returns:

            Tuple[List, List]::
                    (   List of violating indices of the window start,
                        the per window start index counts of repeat_counts)

        Raises:
            ValueError:
        '''
        cdef:
            int violation_count = 0
            int* repeat_violation_idxs = NULL
            int* repeat_counts = NULL
            repeatcheck_t* rcheck = self.rcheck
            size_t i
            list out = []
            list counts = []
            int ret
        if window_size > rcheck.seq_length:
            raise ValueError(
                 f'window_size: {window_size}, must be less than or equal '
                 f'to the sequence length: {rcheck.seq_length}'
            )
        try:
            ret = sr_repeatWindow(
                rcheck,
                window_size,
                max_repeat_count,
                &repeat_violation_idxs,
                &violation_count,
                &repeat_counts,
            )
            if ret < 0:
                raise ValueError('Failure in window calculation')

            for i in range(violation_count):
                out.append(repeat_violation_idxs[i])
            if get_counts:
                for i in range(rcheck.seq_length - window_size + 1):
                    counts.append(
                        repeat_counts[i],
                    )
        finally:
            if repeat_counts != NULL:
                PyMem_Free(repeat_counts)
            if repeat_violation_idxs != NULL:
                PyMem_Free(repeat_violation_idxs)
        return out, counts

    def get_repeat_window(
            self,
            idx: int,
            window_size: int,
    ) -> Dict[str, List[int]]:
        ''' Used for retreaving repeats on a specific window
        Args:
            idx (int): Index to search at
            window_size: (int): Window sizew

        Returns:
            Dict: key by word sequence with values of list of repeat indices of
                that word that occur in the window

        Raises:
            IndexError: Out of range
        '''
        cdef:
            repeatcheck_t* rcheck = self.rcheck
            int *repeat_idx_array = rcheck.repeat_idx_arr
            int word_size = rcheck.word_size
            int lim = rcheck.seq_length - word_size + 1
            size_t i
            int ilim = idx + window_size - word_size + 1
            int test_val
            dict out = {}
            str seq_obj = self.seq_obj
        if idx < 0 or idx > lim:
            raise IndexError(f'index {idx} out of range {lim}')
        if idx + window_size - word_size >= lim:
            raise IndexError(
                f'index {idx} and window_size {window_size} out of '
                f'range {lim}'
            )
        for i in range(idx, ilim):
            test_val = repeat_idx_array[i]
            if test_val > 0 and test_val < ilim:
                seq = seq_obj[i:i + word_size]
                if seq not in out:
                    out[seq] = ilist = [i]
                else:
                    ilist = out[seq]
                ilist.append(test_val)
        return out

    def fraction_repeats(self) -> float:
        '''Compute the total fraction of bases participating in repeats
        over the entire sequence

        Returns:
            float: fraction in range [0, 1]

        Raises:
            OSError: Out of memory
        '''
        cdef:
            repeatcheck_t* rcheck = self.rcheck
            int* repeat_idx_array = rcheck.repeat_idx_arr
            int word_size = rcheck.word_size
            int seq_length = rcheck.seq_length
            int lim = seq_length - word_size + 1
            cnp.uint8_t* count_buf
            size_t i, block_size, j
            int count = 0
            float out


        block_size = seq_length*sizeof(cnp.uint8_t)
        count_buf = <cnp.uint8_t*> PyMem_Malloc(block_size)
        if count_buf == NULL:
            raise OSError('Out of memory')
        memset(count_buf, 0, block_size)
        block_size = word_size*sizeof(cnp.uint8_t)
        for i in range(lim):
            if repeat_idx_array[i] != 0:
                memset(&count_buf[i], 1, block_size)
        for i in range(seq_length):
            count += count_buf[i]
        PyMem_Free(count_buf)
        out = count
        return out / seq_length

    def window_fraction_limit(
            self,
            window_size: int,
            cutoff: float,
    ) -> Tuple[List[Tuple[int, float]], Tuple[int, float]]:
        ''' Compute the fraction of bases in the window of window_size
        participating in repeats over all windows in the sequence
        Args:
            window_size (int): window_size to look for repeats
            cutoff (float): fraction in range [0, 1] at which to capture results

        Returns:
            Tuple(List, Tuple)::

                <List of fractions greater than cutoff>,
                ( <index of max violation>, <max_fraction> )

        Raises:
            ValueError: cutoff out of range 0 to 1
            OSError: Out of memory
        '''
        cdef:
            repeatcheck_t* rcheck = self.rcheck
            int* repeat_idx_array = rcheck.repeat_idx_arr
            int word_size = rcheck.word_size
            int seq_length = rcheck.seq_length
            int check, test_val, count
            size_t max_idx = 0
            float max_fraction = 0.
            float this_fraction
            list out = []
            size_t i, j, jl, block_size1, block_size2
            cnp.uint8_t* count_buf
            int il = seq_length - window_size + 1

        if cutoff > 1. or cutoff < 0.:
            raise ValueError(
                f'cutoff {cutoff} must be in range [0, 1]'
            )
        check = int(window_size*cutoff)

        block_size1 = window_size*sizeof(cnp.uint8_t)
        count_buf = <cnp.uint8_t*> PyMem_Malloc(block_size1)
        if count_buf == NULL:
            raise OSError('Out of memory')
        memset(
            count_buf,
            0,
            block_size1,
        )
        block_size2 = word_size * sizeof(cnp.uint8_t)

        for i in range(il):
            memset(count_buf, 0, block_size1)
            count = 0
            for j in range(i, i + window_size - word_size + 1):
                test_val = repeat_idx_array[j]
                if test_val != 0:
                    memset(&count_buf[j - i], 1, block_size2)
            for j in range(window_size):
                count += count_buf[j]
            if count > check:
                this_fraction = count / window_size
                out.append((i, this_fraction))
                if max_fraction < this_fraction:
                    max_fraction = this_fraction
                    max_idx = i
        PyMem_Free(count_buf)
        return out, (max_idx, max_fraction)
