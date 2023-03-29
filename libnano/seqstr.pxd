# cython: language_level=3, boundscheck=False, wraparound=False
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
libnano.seqstr
~~~~~~~~~~~~~~

Cython header file for seqstr.pyx -- allows for cross-project Cython /
C integration of the low-level thermodynamic analysis bindings.

----

Basic operations of DNA, RNA, and amino acid (AA) sequences (e.g.,
complement, reverse, reverse complement, etc)

Includes hamming distance method and derivative heuristic methods.
'''

cdef:
    char* reverse_seq_c(char* seq, int length)
    void reverse_seq_cb(char* seq, char* out, int length)

    char* complement_c(char* seq, int length)
    void complement_cb(char* seq, char* out, int length)

    char* reverse_complement_c(char* seq, int length)
    void reverse_complement_cb(char* seq, char* out, int length)

    int hamming_distance_c(char* seq1, char* seq2, int len1)

    void rolling_hamming_distance_c(
        char* seq1, char* seq2,
        int len1, int len2,
        int overlap,
        int* hamming_distance_arr,
    )
