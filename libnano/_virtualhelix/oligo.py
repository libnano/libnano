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
libnano._virtualhelix.oligo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''
# import sys
# from typing import (
#     Optional,
# )

# import __init__
# from strand import Strand

# from libnano.seqrecord.seqrecordbase import (
#     ALPHABETS,
#     AlphaEnum,
# )

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
#             Strand_t* strands = <Strand_t*> mem.malloc(
#                  OLIGO_DEFAULT_LEN,
#                   sizeof(Strand_t),
#             ) # order 5' to 3'
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

# class Oligo(object):
#     '''An Oligo'''

#     def __init__(self, seq: str = ''):
#         self.seq: str = seq
#         self.strand5p: Optional[Strand] = Strand(seq, self) if seq else None

#     def __len__(self) -> int:
#         if self.strand5p:
#             return sum(
#                 len(x) for x in self.strand5p.gen3p()
#             )
#         else:
#             return 0

#     def __str__(self) -> str:
#         return f'Oligo {len(self)}: {self.seq}'

#     def __repr__(self) -> str:
#         return f'Oligo <{id(self)}> {len(self)}: {self.seq}'

#     def prepend5p(self, strand: Strand):
#         '''Append a :class:`Strand` to the 3 prime end of the oligo
#         '''
#         old_strand5p: Strand = self.strand5p
#         assert (old_strand5p.strand5p is None)
#         old_strand5p.strand5p = strand
#         self.strand5p = strand
#         strand.strand3p = old_strand5p
