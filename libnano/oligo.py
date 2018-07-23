import sys
from typing import (
    List,
    Generator
)

import __init__
from libnano.seqrecord.seqrecordbase import (
    AlphaEnum,
    ALPHABETS
)

from strand import Strand
# from libnano.cymem.cymem cimport Pool

# cdef char* cseq(self):
#     cdef char* offset = self.offset
#     return self.oligo.seq[offset]

# cdef struct Strand_t:
#     int start   # start index in the Oligo
#     int length  # length of the Strand
#     int vh_id   # Virtual Helix ID Number


# cdef int OLIGO_DEFAULT_LEN = 4

# cdef class Oligo:
#     cdef:
#         Pool mem
#         char* cseq  # pointer to the c string reference in the _
#         int len_seq
#         Strand_t* strands
#         int len_strands

#     def __cinit__(Self, assembly, strand5p):
#         self.mem = None
#         self.strands = NULL
#         self.cseq = NULL
#         self._seq = None
#         self.len_seq = 0

#     def __init__(self, assembly, strand5p: Strand = None):
#         self.mem = mem = Pool()
#         self.assembly = assembly
#         cdef:
#             Strand_t* strands = <Strand_t*>mem.malloc(OLIGO_DEFAULT_LEN, sizeof(Strand_t)) # order 5' to 3'
#             self.strands = strands
#         if strand5p is not None:
#             strands[0].start = 0
#             strands[0].length = strand5p.length
#             strands[0].vh_id = strand5p.id_num
#             self.len_seq = strand5p.length
#     # end def

#     @property
#     def seq(self) -> str:
#         return self._seq = seq

#     @seq.setter
#     def seq(self, value: str):
#         self._seq = value
#         self.cseq =  c_util.obj_to_cstr_len(value, &self.len_seq)

#     @seq.deleter
#     def seq(self):
#         self.cseq = NULL
#         del self._seq = None
#         self.len_seq = 0

#     cdef char* cseq_index(self, int index):
#         cdef int offset = self.strands[index].start
#         return &self.cseq[offset]

#     cdef append3p(self, strand: Strand):
#         '''Append a :class:`Strand` to the 3 prime end of the oligo
#         '''
#         strand5p: Strand = self.strand5p
#         strand3p: Strand = strand5p.strand3p
#         while strand3p is not None:
#             strand5p = strand3p
#             strand3p = strand5p.strand3p
#         strand5p.strand3p = strand
#     # end def

class  Oligo(object):
    """An Oligo"""
    def __init__(self, seq: str = None):
        self.seq: str = seq
        strand: Strand = None
        if seq is not None:
            strand = Strand(seq, self)
        self.strand5p: Strand = strand

    def __len__(self) -> int:
        return sum(len(x) for x in self.strand5p.gen3p())
    # end def

    def __str__(self) -> str:
        return "Oligo %d: %s" % (len(self), self.seq)

    def __repr__(self) -> str:
        return  "Oligo <%d> %d: %s" % (id(self), len(self), self.seq)

    def prepend5p(self, strand: Strand):
        '''Append a :class:`Strand` to the 3 prime end of the oligo
        '''
        old_strand5p: Strand = self.strand5p
        assert(old_strand5p.strand5p is None)
        old_strand5p.strand5p = strand
        self.strand5p = strand
        strand.strand3p = old_strand5p
    # end def
# end class