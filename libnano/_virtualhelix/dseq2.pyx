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
libnano._virtualhelix.dseq2
~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''
# cimport libnano.seqstr
from libnano.cymem.cymem cimport Pool
from libnano.helpers cimport c_util


cdef class DSeq2:
    cdef int* fwd_gaps
    cdef int* rev_gaps
    cdef char* fwd_cstr
    cdef Py_ssize_t fwd_cstr_len
    cdef char* rev_cstr
    cdef Py_ssize_t rev_cstr_len

    cdef Pool mem   # Memory manager

    def __cinit__(self, fwd: str, rev: str):
        self.mem = None
        self.fwd_gaps = NULL
        self.rev_gaps = NULL
        self.fwd_cstr = NULL
        self.fwd_cstr_len = 0
        self.rev_cstr = NULL
        self.rev_cstr_len = 0


    def __init__(self, fwd: str, rev: str):
        self.mem = mem = Pool()
        self.fwd_cstr = c_util.obj_to_cstr_len(
            fwd,
            &self.fwd_cstr_len,
        )
        self.rev_cstr = c_util.obj_to_cstr_len(
            rev,
            &self.rev_cstr_len,
        )

    def awesome(self):
        print("AWESOME")
# end class
