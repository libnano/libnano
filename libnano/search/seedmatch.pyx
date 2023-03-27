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
libnano.search.seedmatch
~~~~~~~~~~~~~~~~~~~~~~~~

'''
from typing import (
    List,
    Tuple,
    Union,
)

from cython.operator cimport postincrement as inc
from libc.stdio cimport printf
from libc.string cimport memcpy

from libnano.cymem.cymem cimport Pool
from libnano.helpers cimport c_util

STR_T = Union[str, bytes]

cdef extern from "shl_seedhashlist.h":
    # leave out a few members
    ctypedef struct tuple_t:
        int idx
        int offset

    ctypedef struct seed_search_obj_t:
        pass

    int shl_seed_weight(char*, int)
    int shl_parse_seed(char*,int,int*)
    int shl_free_hash_obj(seed_search_obj_t*)
    int shl_buildSeedTable(
        int*, int, int,
        char*, int, int*,
        int*, int, seed_search_obj_t*
    )
    int shl_findSingleMatches(
        int, char*, int,
        char*, int, int*,
        int*, seed_search_obj_t*, tuple_t**,
        int,
    )
    int checkht(seed_search_obj_t*)

cdef class SeedMatcher:
    '''Generate the hash table from a list of
    sequences for submers of lenght m with k mismatches

    uses seed.pyx to generate seeds or seed pairs

    '''
    cdef:
        list reference_list
        Py_ssize_t reference_list_length

        char* ref_cstr_buffer  # buffer to put all c-strings of
        Py_ssize_t ref_cstr_buf_length

        int* reference_lengths     # length of each string in library
        int* reference_start_idxs  # start indices of each string in the library

        int* seed_idxs
        int seed_idxs_len

        object seed
        Py_ssize_t seed_len

        char* seed_cstr
        seed_search_obj_t* search_obj

        Pool mem   # Memory manager

    def __cinit__(
            self,
            seed: STR_T,
            reference_list: List[STR_T],
            debug: bool = False,
        ):
        self.mem = None
        self.reference_list = None
        self.reference_list_length = 0
        self.seed = None
        self.seed_len = 0

        # the concatenated c string version of the reference, null terminated
        self.ref_cstr_buffer = NULL

        # this is length of each string in reference_list
        self.reference_lengths = NULL

        # this is the start index of each string in reference_list in
        #       self.ref_cstr_buffer
        self.reference_start_idxs = NULL


        self.seed_idxs = NULL
        self.seed_cstr = NULL
        self.search_obj = NULL

    def __init__(
            self,
            seed: STR_T,
            list reference_list,
            debug: bool = False,
    ):
        '''
        Args:
            seed: Seed string
            reference_list: List of references
            debug:  If True, priner debug messages
        '''
        cdef:
            Py_ssize_t i, j, k
            Py_ssize_t reference_list_length, total_seq_length
            Py_ssize_t seed_len
            int offset_last, check

            char* temp_ref_cstr_buff = NULL
            int* temp_ref_lengths = NULL
            int* temp_ref_start_idxs = NULL
            char* seq_cstr = NULL
            Py_ssize_t seq_cstr_len

        self.mem = mem = Pool()
        self.seed = seed

        # 1 Setup the seed
        self.seed_cstr = c_util.obj_to_cstr_len(
            seed,
            &seed_len,
        )
        self.seed_idxs_len = shl_seed_weight(
            self.seed_cstr,
            seed_len,
        )
        self.seed_idxs = <int*> mem.malloc(
            self.seed_idxs_len,
            sizeof(int),
        )

        shl_parse_seed(
            self.seed_cstr,
            <int> seed_len,
            self.seed_idxs,
        )

        self.seed_len = seed_len

        # 2. Get constants from the list of sequences
        self.reference_list = reference_list
        reference_list_length = len(reference_list)
        total_seq_length = sum(
            len(x) for x in reference_list
        )


        # 3. Memory allocations
        self.ref_cstr_buf_length = total_seq_length + reference_list_length
        temp_ref_cstr_buff = <char*> mem.malloc(
            self.ref_cstr_buf_length,
            sizeof(char),
        )
        temp_ref_lengths = <int*> mem.malloc(
            reference_list_length,
            sizeof(int),
        )
        temp_ref_start_idxs = <int*> mem.malloc(
            reference_list_length,
            sizeof(int),
        )

        self.search_obj = <seed_search_obj_t *> mem.calloc(
            1,
            sizeof(seed_search_obj_t)
        )

        # 4. Copy strings into one buffer
        i = 0
        k = 0
        offset_last = 0
        for seq in reference_list:
            seq_cstr = c_util.obj_to_cstr_len(
                seq,
                &seq_cstr_len,
            )

            temp_ref_lengths[k] = <int> inc(seq_cstr_len)  # Notice the inc
            temp_ref_start_idxs[inc(k)] = offset_last

            offset_last += <int> seq_cstr_len
            # Copy everything including the NULL
            memcpy(
                &temp_ref_cstr_buff[i],
                seq_cstr,
                seq_cstr_len,
            )
            i += seq_cstr_len

        # Assign pointers
        self.ref_cstr_buffer = temp_ref_cstr_buff
        self.reference_lengths = temp_ref_lengths
        self.reference_start_idxs = temp_ref_start_idxs

        if debug:
            print("building table", self.seed_idxs[0])
        check = shl_buildSeedTable(
            self.seed_idxs,
            self.seed_idxs_len,
            self.seed_len,
            self.ref_cstr_buffer,
            self.ref_cstr_buf_length,
            self.reference_lengths,
            self.reference_start_idxs,
            reference_list_length,
            self.search_obj,
        )
        if debug:
            print("table built")
        if check < 0:
            raise OSError("something happened")

    def __dealloc__(self):
        #print("freeing")
        shl_free_hash_obj(self.search_obj)


    cdef list toList(self, tuple_t* matches, int repeat_count):
        cdef:
            Py_ssize_t i
            list outlist = [None] * repeat_count

        for i in range(repeat_count):
            outlist[i] = (
                matches[i].idx,
                 matches[i].offset,
                )
        return outlist


    def _checkht(self):
        checkht(self.search_obj)

    def match(
            self,
            target: STR_T,
            mismatches: int,
    ) -> List[Tuple[int, int]]:
        cdef:
            int repeat_count
            char* target_cstr = NULL
            Py_ssize_t target_len
            tuple_t* matches = NULL

            int matches_buffer_size = self.ref_cstr_buf_length

        target_cstr = c_util.obj_to_cstr_len(
            target,
            &target_len,
        )

        repeat_count = shl_findSingleMatches(
            mismatches,
            target_cstr,
            target_len,
            self.ref_cstr_buffer,
            self.ref_cstr_buf_length,
            self.reference_lengths,
            self.reference_start_idxs,
            self.search_obj,
            &matches,
            matches_buffer_size,
        )
        if repeat_count > 0:
            #print("Found", repeat_count)
            self.mem.own(
                matches,
                repeat_count,
                sizeof(tuple_t),
            )
            outlist = self.toList(
                matches,
                repeat_count,
            )
            self.mem.free(matches)
            return outlist
        else:
            return []
