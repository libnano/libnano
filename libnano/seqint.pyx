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
libnano.seqint
~~~~~~~~~~~~~~

Compact integer representations of DNA sequences
'''
cimport c_util
from cpython.mem cimport (  # PyMem_Realloc,
    PyMem_Free,
    PyMem_Malloc,
)

from typing import Union

STR_T = Union[str, bytes]

cdef extern from "si_seqint.h":
    uint64_t si_seq2Int(char*, int)
    int si_int2Seq(uint64_t, char*, int)
    uint64_t si_revSeqInt(uint64_t, int)
    uint64_t si_compSeqInt(uint64_t, int)
    uint64_t si_revCompSeqInt(uint64_t, int)
    uint64_t si_addBase(uint64_t, char)
    uint64_t si_removeBase(uint64_t)
    uint64_t si_seqIntSubstring(uint64_t, int, int, int)
    uint64_t si_addToWindow(uint64_t, char, int)


def tester(s: STR_T) -> str:
    cdef:
        Py_ssize_t slen
        char* s_copy_cstr
    c_util.copy_obj_to_cstr(s, &slen, &s_copy_cstr)
    return c_util.cstr_to_obj(s_copy_cstr, slen, 0)


def seq2Int(
        seq_obj: STR_T,
) -> int:
    '''Convert a DNA sequence to an integer representation

    Args:
        seq: DNA sequence comprised of the bases A/a, T/t, G/g, and C/c

    Returns:
        Binary integer representation of sequence

    Raises:
        ValueError: Sequences over 30 bp cannot be converted
    '''
    cdef:
        char*           seq
        uint64_t        seqint
        int length

    length = <int> len(seq_obj)

    if length > 30:
        raise ValueError(
            f'seq2Int: Sequences over 30 bp cannot be converted (length)'
        )


    seq = c_util.obj_to_cstr(seq_obj)

    ``seqint`` = si_seq2Int(seq, length)
    return seqint


def reverseComplement(
        uint64_t seqint,
        int length,
) -> float:
    '''Return the reverse complement ``seqint`` representation of ``seqint``

    Args:
        seqint: integer representation of a DNA sequence
        length: length of the represented sequence in bases

    Returns:
        Reverse complement ``seqint`` representation of ``seqint``
    '''
    cdef uint64_t         seqintrc

    seqintrc = si_revCompSeqInt(seqint, length)

    return seqintrc


def complement(
        uint64_t seqint,
        int length,
) -> float:
    '''Return the complement ``seqint`` representation of ``seqint``

    Args:
        seqint: integer representation of a DNA sequence
        length: length of the represented sequence in bases

    Returns:
        The complement ``seqint`` representation of ``seqint``

    '''
    cdef uint64_t         seqintrc

    seqintrc = si_compSeqInt(seqint, length)

    return seqintrc


def reverse(
        uint64_t seqint,
        int length,
) -> float:
    '''Return the reverse complement ``seqint`` representation of ``seqint``

    Args:
        seqint: integer representation of a DNA sequence
        length: length of the represented sequence in bases

    Returns:
        The reverse complement ``seqint`` representation of ``seqint``
    '''
    cdef uint64_t         seqintrc

    seqintrc = si_revSeqInt(seqint, length)

    return seqintrc


def addBase(uint64_t seqint, base_obj: STR_T) -> int:
    '''Add a base to the right-hand side of a seqint

    Args:
        seqint: integer representation of a DNA sequence
        base_obj: base (A/a, T/t, G/g, or C/c) to be added

    Return:
        ``seqint`` with righthand appended base
    '''
    cdef:
        uint64_t    seqintmod
        bytes base_py = c_util._bytes(base_obj)
        char*       base = base_py

    seqintmod = si_addBase(seqint, base[0])

    return seqintmod


def removeBase(
        uint64_t seqint,
        int length,
) -> int:
    '''Remove a base from the right-hand side of a seqint

    Args:
        seqint: integer representation of a DNA sequence
        length: length of the represented sequence in bases

    Returns:
        seqint with base removed
    '''
    cdef uint64_t    seqintmod

    seqintmod = si_removeBase(seqint)

    return seqintmod


def addToWindow(
        uint64_t seqint,
        base_obj: STR_T,
        int length,
) -> int:
    '''Add a base to the right-hand side of a ``seqint`` (fixed length window)

    Args:
        seqint: integer representation of a DNA sequence
        base: base (A/a, T/t, G/g, or C/c) to be added
        length: length of the represented sequence in bases

    Returns:

    '''
    cdef:
        char*       base
        uint64_t    seqintmod

        bytes base_py = c_util._bytes(base_obj)

    base = base_py

    seqintmod = si_addToWindow(
        seqint,
        base[0],
        length,
    )

    return seqintmod


def getSubstring(
        uint64_t seqint,
        int start_idx,
        int end_idx,
        int length,
) -> int:
    '''Returns the ``seqint`` representing the defined substring

    Args:
        seqint: integer representation of a DNA sequence
        start_idx: substring start index (inclusive, 0-based)
        start_idx: substring end index (exclusive)
        length: length of the represented sequence in bases

    Returns:
        Integer representation of string
    '''
    cdef uint64_t    seqintmod

    seqintmod = si_seqIntSubstring(
        seqint,
        start_idx,
        end_idx,
        length,
    )

    return seqintmod


def int2Seq(
        uint64_t seqint,
        int length,
) -> str:
    '''Return the DNA sequence string represented by ```seqint```

    Args:
        seqint: integer representation of a DNA sequence
        length: length of the represented sequence in bases

    Returns:
        DNA sequence string represented by ```seqint```
    '''
    cdef object ret_obj
    cdef char* out_seq

    if length == 0:
        return ''

    out_seq = <char *> PyMem_Malloc((length+1)*sizeof(char))
    if out_seq == NULL:
        raise OSError('Could not allocate memory for sequence.')

    if si_int2Seq(seqint, out_seq, length) != 0:
        PyMem_Free(out_seq)
        raise OSError('Could not convert integer to sequence.')

    try:
        ret_obj = out_seq.decode('utf8')
    finally:
        PyMem_Free(out_seq)

    return ret_obj
# end def

cdef inline int2Seq_c(uint64_t seqint, char* out, int length):
    '''Cython header exposed wrapper for si_int2Seq
    '''
    return si_int2Seq(seqint, out, length)
