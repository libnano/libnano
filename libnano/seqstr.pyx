#cython: boundscheck=False, wraparound=False
"""libnano.seqstr
~~~~~~~~~~~~~~

Basic operations of DNA, RNA, and amino acid (AA) sequences (e.g.,
complement, reverse, reverse complement, etc)

Includes hamming distance method and derivative heuristic methods.
"""


import random
from typing import (
    Tuple,
    Union,
    Iterable,
    List
)

from cython.operator cimport postincrement as inc

import numpy as np
cimport numpy as cnp
cnp.import_array()
from cpython.mem cimport (
    PyMem_Malloc,
    PyMem_Free
)

cimport c_util

STR_T = Union[str, bytes]

__doc__ = "Basic operations on DNA, RNA, and amino acid sequences"

cdef extern from "ss_seqstr.h":
    int ss_revSeq(char*, int)
    int ss_compSeq(char*, int)
    int ss_revCompSeq(char*, int)
    int ss_revCompSeqCopy(char*, const char*, int)
    int ss_hamming(const char*, const char*, int)
    int ss_rollingHamming(const char*, const char*,
                        int, int, int, int*)
    int ss_minHammingThreshold(const char*, const char*, int, int, int)

def randomer(int length) -> str:
    """Create a randomer of length `length`

    Args:
        length (int): length of the randomer

    Returns:
        str
    """
    choices = ('A', 'C', 'G', 'T')
    f = random.choice
    return ''.join([f(choices) for i in range(length)])
# end def

def dsRandomer(int length, bint do_print = False) -> Tuple[str, str]:
    """Create a double stranded randomer

    Args:
        length (int): length of the randomer
        do_print (Optional, bool): False by default

    Returns:
        tuple of strings (str)
    """
    s = randomer(length)
    comp = complement(s)
    if do_print:
        print(s + '\n\r' + comp)
    return s, comp
# end def

def reverse(seq: STR_T) -> STR_T:
    """Compute the reverse of a DNA/RNA/amino acid sequence

    Args:
        seq: DNA, RNA, or amino acid sequence

    Returns:
        The reverse of the sequence (as a string)
    """
    cdef Py_ssize_t s_length
    cdef char* rev_seq_cstr
    cdef object rev_seq
    cdef int obj_type

    obj_type = c_util.copy_obj_to_cstr(seq, &s_length, &rev_seq_cstr)

    ss_revSeq(rev_seq_cstr, <int> s_length)

    rev_seq = c_util.cstr_to_obj(rev_seq_cstr, s_length, obj_type)

    return rev_seq
# end def


cdef inline char* reverse_c(char* seq, int length):
    cdef char* r_seq = c_util.copy_string(seq, length)
    ss_revSeq(r_seq, length)
    return r_seq
# end def


cdef inline void reverse_cb(char* seq, char* out, int length):
    c_util.copy_string_buffer(seq, out, length)
    ss_revSeq(out, length)
# end def


def complement(seq: STR_T) -> STR_T:
    """Compute the complement of a DNA/RNA/amino acid sequence

    Args:
        seq (str): DNA, RNA, or amino acid sequence

    Returns:
        The complement of the sequence (as a string)
    """

    cdef Py_ssize_t s_length
    cdef char* comp_seq_cstr
    cdef object comp_seq
    cdef int obj_type
    obj_type =  c_util.copy_obj_to_cstr(seq, &s_length, &comp_seq_cstr)

    ss_compSeq(comp_seq_cstr, <int> s_length)

    comp_seq = c_util.cstr_to_obj(comp_seq_cstr, s_length, obj_type)
    return comp_seq
# end def


cdef inline char* complement_c(char* seq, int length):
    cdef char* c_seq = c_util.copy_string(seq, length)
    ss_compSeq(c_seq, length)
    return c_seq
# end def


cdef inline void complement_cb(char* seq, char* out, int length):
    c_util.copy_string_buffer(seq, out, length)
    ss_compSeq(out, length)
# end def


def reverseComplement(seq: STR_T) -> STR_T:
    """Compute the reverse complement of a DNA/RNA/amino acid sequence

    Args:
        seq (str): DNA, RNA, or amino acid sequence

    Returns:
        The reverse complement of the sequence (as a string)
    """
    cdef Py_ssize_t s_length
    cdef char* rcomp_seq_cstr
    cdef object rcomp_seq
    cdef int obj_type

    obj_type = c_util.copy_obj_to_cstr(seq, &s_length, &rcomp_seq_cstr)

    ss_revCompSeq(rcomp_seq_cstr, <int> s_length)

    rcomp_seq = c_util.cstr_to_obj(rcomp_seq_cstr, s_length, obj_type)
    return rcomp_seq
# end def

cdef inline char* reverseComplement_c(char* seq, int length):
    cdef char* c_seq = c_util.copy_string(seq, length)
    ss_revCompSeq(c_seq, length)
    return c_seq
# end def


cdef inline void reverseComplement_cb(char* seq, char* out, int length):
    c_util.copy_string_buffer(seq, out, length)
    ss_revCompSeq(out, length)
# end def


def hammingDistance(seq1: STR_T, seq2: STR_T) -> int:
    """Compute the Hamming distance between two DNA/RNA/AA sequences

    Args:
        seq1 (str): DNA/RNA/AA sequence
        seq2 (str): DNA/RNA/AA sequence

    Returns:
        The hamming distance between the sequences as an integer
    """
    cdef char* seq1_cstr
    cdef char* seq2_cstr
    cdef Py_ssize_t len1

    seq2_cstr = c_util.obj_to_cstr_len(seq2, &len1)
    seq1_cstr = c_util.obj_to_cstr_len(seq1, &len1)

    return ss_hamming(  <const char *> seq1_cstr,
                        <const char *> seq2_cstr, <int>len1)
# end def

cdef inline int hammingDistance_c(char* seq1, char* seq2, int len1):
    """Compute the Hamming distance between two strings
    seq: DNA, RNA, or amino acid sequence
    """
    return ss_hamming(  <const char *> seq1,
                        <const char *> seq2, <int>len1)
# end def


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simple heuristics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def minHammingDistance(seq1: STR_T, seq_list: List[STR_T]) -> int:
    """Compute the min. hamming distance between a seq and as list of seqs.

    All of the sequences in `seq_list` must have the same length as
    `seq1` (raises ValueError).

    Args:
        seq1 (str): DNA/RNA/AA sequence for min. hamming calculations
        seq_list (list): List of DNA/RNA/AA sequences as targets for min.
                         hamming calculations

    Returns:
        Integer representing the minimum hamming distance between `seq1` and
        any sequence in `seq_list`

    Raises:
        ValueError
    """
    cdef char* seq1_cstr
    cdef char* seq2_cstr
    cdef Py_ssize_t len1
    cdef Py_ssize_t len2
    cdef object seq2

    seq1_cstr = c_util.obj_to_cstr_len(seq1, &len1)

    min_hd = len1

    for seq2 in seq_list:
        seq2_cstr = c_util.obj_to_cstr_len(seq2, &len2)
        if not len2 == len1:
            raise ValueError('All strings in `seq_list` must have the same'
                             ' length as `seq1`')
        min_hd = min(min_hd, ss_hamming(<const char *> seq1_cstr,
                                        <const char *> seq2_cstr, <int>len1))


    return min_hd
# end def

def thresholdRollingHammingList(target: STR_T,
                                seq_list: List[STR_T],
                                int threshold) -> int:
    """Compute the min. hamming distance between a seq and as list of seqs.

    All of the sequences in `seq_list` must have the same length as
    `target` (raises ValueError).

    Args:
        target (str): DNA/RNA/AA sequence for min. hamming calculations
        seq_list (list): List of DNA/RNA/AA sequences as targets for min.
                         hamming calculations
        threshold (int): threshold value to be greater than or equal to
            to not be a violation

    Returns:
        int: the minimum hamming distance between `target` and
            any sequence in `seq_list`

    Raises:
        ValueError
    """
    cdef char* target_cstr
    cdef char* ref_cstr
    cdef Py_ssize_t target_len
    cdef Py_ssize_t ref_len
    cdef object ref
    cdef int num_positions, hamming_distance
    cdef Py_ssize_t i, j, k
    target_cstr = c_util.obj_to_cstr_len(target, &target_len)
    cdef int tlim = threshold + 1

    cdef list out_list = []
    for i in range(len(seq_list)):
        ref = seq_list[i]
        ref_cstr = c_util.obj_to_cstr_len(ref, &ref_len)
        #1. Do the left end of reference (requires for threshold >= 2 to do anything)
        for j in range(threshold - 1):
            hamming_distance = 0
            for k in range(target_len - j - 1):
                if target_cstr[j + 1 + k] != ref_cstr[k]:
                    inc(hamming_distance)
            if hamming_distance < tlim:
                out_list.append((i, -(j+1)))

        #2. Do the middle
        num_positions = ref_len - target_len + 1
        for j in range(num_positions):
            hamming_distance = 0
            with nogil:
                for k in range(target_len):
                    if target_cstr[k] != ref_cstr[j + k]:
                        inc(hamming_distance)
            if hamming_distance < tlim:
                out_list.append((i, j))
        #3. Do the right end
        num_positions = ref_len - target_len
        for j in range(threshold - 1):
            hamming_distance = 0
            for k in range(target_len - j - 1):
                if target_cstr[k] != ref_cstr[num_positions + j + 1 + k]:
                    inc(hamming_distance)
            if hamming_distance < tlim:
                out_list.append((i, num_positions + j + 1))
    return out_list
# end def

def rollingHammingDistance(seq1: STR_T, seq2: STR_T, int overlap=0) -> np.ndarray:
    """Compute the Hamming distance between seq1 and seq2 at each idx of seq2

    Args:
        seq1 (str): DNA/RNA/AA sequence for rolling hamming comparison
        seq2 (str): DNA/RNA/AA target for rolling hamming comparison
        overlap (int): the number of bases to overlap seq1 onto seq2 at the 3'
            and 5' ends. So a 10-mer seq1 with overlap of 8 would overhang seq2
            by 2 bases

    Returns:
        A numpy array of the hamming distance between seq1 and seq2 computed
        at each index of seq2

    Raises:
        ValueError


    For example:
        seq1: ATGCCG
        seq2: TAGCCGGACCGTTAGGACCACGTA

        ATGCCG
        TAGCCGGACCGTTAGGACCACGTA
        idx 0 hamming distance: 2

         ATGCCG
        TAGCCGGACCGTTAGGACCACGTA
        idx 1 hamming distance: 3
                   ...
                          ATGCCG
        TAGCCGGACCGTTAGGACCACGTA
        idx 18 hamming distance: 6

        returned numpy array: [2, 3, ..., 6]
            length: len(seq2)-len(seq1) + 1
    """
    cdef char* seq1_cstr
    cdef char* seq2_cstr
    cdef Py_ssize_t len1, len2

    cdef cnp.ndarray[cnp.int32_t, ndim=1] hamming_distance_np_arr

    cdef int* hamming_distance_arr
    cdef int num_positions

    seq1_cstr = c_util.obj_to_cstr_len(seq1, &len1)
    seq2_cstr = c_util.obj_to_cstr_len(seq2, &len2)

    num_positions = <int> (len2 - len1 + 1 +2*overlap)

    hamming_distance_np_arr = np.zeros(num_positions, dtype=np.int32)

    hamming_distance_arr = <int*> cnp.PyArray_DATA(hamming_distance_np_arr)

    ss_rollingHamming(<const char *> seq1_cstr, <const char *> seq2_cstr,
                      <int> len1, <int>len2, overlap, hamming_distance_arr)

    # get rid of PyObject* to get cython objects
    return hamming_distance_np_arr
# end def

cdef inline void rollingHammingDistance_c(char* seq1, char* seq2,
                                          int len1, int len2, int overlap, int*
                                          hamming_distance_arr):
    ss_rollingHamming(<const char *> seq1, <const char *> seq2,
                      len1, len2, overlap, hamming_distance_arr)
# end def

def misprimeCheck(  putative_seq: STR_T,
                    sequences: Iterable[STR_T],
                    int hamming_threshold) -> bool:
    """Calculate the heterodimer formation thermodynamics of a DNA
    sequence, ``putative_seq`` with a list of sequences relative to
    a melting temperature threshold

    Args:
        putative_seq (str):
        sequences (iterable of str): each sequence should be of length
                                    greater than or equal to putative_seq
        hamming_threshold (int): value if less than will be equivalent to a misprime

    Returns:
        bool: True if will misprime False otherwise.
    """
    cdef bint is_offtarget = False
    cdef int res
    cdef char* seq1_cstr    # the shorter sequence
    cdef char* seq2_cstr    # the longer sequence
    cdef Py_ssize_t len1, len2

    seq1_cstr = c_util.obj_to_cstr_len(putative_seq, &len1)
    for seq in sequences:
        seq2_cstr = c_util.obj_to_cstr_len(seq, &len2)
        res = ss_minHammingThreshold(<const char *> seq1_cstr,
                                    <const char *> seq2_cstr,
                                    <int> len1, <int>len2, hamming_threshold)
        if res < hamming_threshold:
            return True
    return False
# end def