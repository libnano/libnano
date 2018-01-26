"""
libnano.seqstr
~~~~~~~~~~~~~~

Cython header file for seqstr.pyx -- allows for cross-project Cython /
C integration of the low-level thermodynamic analysis bindings.

----

Basic operations of DNA, RNA, and amino acid (AA) sequences (e.g.,
complement, reverse, reverse complement, etc)

Includes hamming distance method and derivative heuristic methods.
"""

cdef char* reverse_c(char* seq, int length)
cdef void reverse_cb(char* seq, char* out, int length)

cdef char* complement_c(char* seq, int length)
cdef void complement_cb(char* seq, char* out, int length)

cdef char* reverseComplement_c(char* seq, int length)
cdef void reverseComplement_cb(char* seq, char* out, int length)

cdef int hammingDistance_c(char* seq1, char* seq2, int len1)

cdef void rollingHammingDistance_c(char* seq1, char* seq2,
                                        int len1, int len2,
                                        int overlap,
                                        int* hamming_distance_arr)
