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
# '''
# libnano._virtualhelix.virtualhelix
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# '''
# from __future__ import annotations

# import sys
# from enum import IntEnum
# from typing import (
#     List,
#     NamedTuple,
#     Optional,
#     Tuple,
#     Union,
# )

# # from dseq import DSeq
# from oligo import Oligo
# from strand import Strand

# # from libnano.fileio.naformat import (
# #     five_prime_type,
# #     three_prime_type,
# # )
# from libnano.seqrecord.seqrecordbase import (  # type: ignore
#     AlphaEnum,
# )
# from libnano.seqstr import (  # type: ignore
#     reverse,
#     reverseComplement,
# )

# NEWLINE_STR: str = '\r\n' if sys.platform == 'win32' else '\n'
# STR_OR_STRAND_T = Union[str, Strand]


# def str2Oligo(x: str) -> Tuple[bool, Oligo, Strand]:
#     if isinstance(x, Strand):
#         return False, x.oligo, x
#     else:
#         oligo = Oligo(x)
#         return True, oligo, oligo.strand5p


# class StrandArray(NamedTuple):
#     strands: List[Strand]
#     idx_offsets: List[int]
#     seq: List[str]


# class VHDirEnum(IntEnum):
#     Forward: int = 0
#     Reverse: int = 1

#     @staticmethod
#     def check(x: int):
#         if x not in (0, 1):
#             raise ValueError(f'{x} not in value in VHDirEnum')


# class VirtualHelix:
#     FWD: int = 0
#     REV: int = 1

#     def __init__(
#             self,
#             fwd_strands: List[Strand],
#             fwd_idx_offsets: List[int],
#             rev_strands: List[Strand],
#             rev_idx_offsets: List[int],
#             alphabet: int = AlphaEnum.DNA,
#     ):
#         '''
#         Args:
#             fwd_strands: List of strands in the 5' to 3' direction
#             rev_strands: List of strands in the 3' to 5' direction
#             fwd_idx_offsets: List of integer offsets from the 5' end of the
#                 forward strand
#             rev_idx_offsets: List of integer offsets from the 5' end of the
#                 forward strand
#             alphabet: DNA or RNA
#         '''
#         self.strand_arrays: Tuple[StrandArray, ...] = (
#             StrandArray(
#                 fwd_strands,
#                 fwd_idx_offsets,
#                 [''],
#             ),
#             StrandArray(
#                 rev_strands,
#                 rev_idx_offsets,
#                 [''],
#             ),
#         )
#         self.alphabet: int = alphabet

#         # fwd_gaps: list = []  # track the free space in the VirtualHelix
#         # fwd_endpoints: list = []
#         # rev_gaps: list = []  # track the free space in the VirtualHelix
#         # rev_endpoints: list = []

#     def oligos(self) -> Tuple[List[Oligo], List[Oligo]]:
#         strand_arrays: Tuple[StrandArray, ...] = self.strand_arrays
#         strands_f: List[Strand] = strand_arrays[0].strands
#         strands_r: List[Strand] = strand_arrays[0].strands
#         out = (
#             [x.oligo for x in strands_f],
#             [x.oligo for x in strands_r],
#         )
#         return out

#     @staticmethod
#     def reverseStrandArray(
#             strands: List[Strand],
#             idx_offsets: List[int],
#     ) -> StrandArray:
#         '''
#         Args:
#             strands: List of :class:`Strand`s
#             idx_offsets: offsets per strand TODO: clarify this

#         Returns:
#             Reversed :class:`StrandArray`
#         '''
#         strands_out: List[Strand] = strands[::-1]
#         # seq_len: int = sum(len(x) for x in strands)
#         x0: int = idx_offsets[0] + len(strands[0])
#         for x1, strand1 in zip(idx_offsets[1:], strands[1:]):
#             check = x1 + len(strand1)
#             assert (check >= x0)
#             x0 = check
#         total_length: int = x0
#         gener = zip(idx_offsets[::-1], strands_out)
#         idx_offsets_out: List[int] = [
#             total_length - (y + len(z) - 1) for y, z in gener
#         ]
#         return StrandArray(
#             strands_out,
#             idx_offsets_out,
#             [''],
#         )

#     def get_seq(
#         self,
#         dir_idx: int,
#         do_cache: bool = False,
#     ) -> str:
#         '''
#         Args:
#             dir_idx: Directiom of interest:
#                       {VirtualHelix.FWD, VirtualHelix.REV}
#             do_cache: Cache the results

#         Returns:
#             Concatenated sequence in the direction of interests
#         '''
#         strand_array: StrandArray = self.strand_arrays[dir_idx]
#         the_seq: str = strand_array.seq[0]
#         if the_seq:
#             return the_seq
#         else:
#             strand_array.strands
#             strand_seq_list: List[str] = []
#             idx_last: int = 0
#             for strand, idx in zip(
#                 strand_array.strands,
#                 strand_array.idx_offsets,
#             ):
#                 strand_seq_list.append(' ' * (idx - idx_last))
#                 seq = strand.seq
#                 strand_seq_list.append(seq)
#                 idx_last = idx + len(seq)
#             the_seq = ''.join(strand_seq_list)
#             if do_cache:
#                 strand_array.seq[0] = the_seq
#             return the_seq

#     @property
#     def fwd_strands(self) -> List[Strand]:
#         return self.strand_arrays[self.FWD].strands

#     @property
#     def rev_strands(self) -> List[Strand]:
#         return reverse(self.strand_arrays[self.REV].strands)

#     def fwd_seq(self, do_cache=True) -> str:
#         return self.get_seq(
#             self.FWD,
#             do_cache=do_cache,
#         )

#     def rev_seq(self, do_cache=True) -> str:
#         return self.get_seq(
#             self.REV,
#             do_cache=do_cache,
#         )

#     def len_strands(
#         self,
#         dir_idx: int,
#     ) -> int:
#         '''
#         Args:
#             dir_idx:

#         Returns:
#             The length of all the strands including the offsets
#             while checking to make sure strands do not overlap
#         '''
#         strand_array: StrandArray = self.strand_arrays[dir_idx]
#         strands: List[Strand] = strand_array.strands
#         idx_offsets: List[int] = strand_array.idx_offsets

#         # seq_len: int = sum(len(x) for x in strands)
#         x0: int = idx_offsets[0] + len(strands[0])
#         for x1, strand1 in zip(idx_offsets[1:], strands[1:]):
#             check = x1 + len(strand1)
#             assert (check >= x0)
#             x0 = check
#         return x0

#     def len_fwd(self) -> int:
#         return self.len_strands(self.FWD)

#     def len_rev(self) -> int:
#         return self.len_strands(self.REV)

#     def __str__(self) -> str:
#         return '%s%s%s' % (
#             self.fwd_seq(),
#             NEWLINE_STR,
#             reverse(self.rev_seq()),
#         )

#     def addForwardStrand(self, strand: Strand, offset: int):
#         pass

#     # def addSequence(self, seq: str) -> Oligo:
#     #     pass

#     def breakStrand(
#         self,
#         dir_idx: int,
#         strand: Strand,
#         idx: int,
#     ) -> Tuple[Oligo, Oligo]:
#         '''Break a Strand in two and create two new Oligos to assign the all
#         of the strands in the pre-existing Oligo to

#         Args:
#             dir_idx: is this on the forward [0] or reverse [1] direction of
#                 the :class:`VirtualHelix`.  Also use :enum:`VHDirEnum` to get
#                 these idxs
#             strand: :class:`Strand` object to break
#             idx: index to break the strand at in terms of it's sequence

#         Returns:
#             Two new :class:`Oligo` objects of form::

#                  Oligo_5p,  Oligo_3p
#         '''
#         VHDirEnum.check(dir_idx)
#         strand_array = self.strand_arrays[dir_idx]
#         vh_strands: List[Strand] = strand_array.strands

#         if strand not in vh_strands:
#             raise ValueError(
#                 f'Strand {strand} not in the {VHDirEnum(dir_idx).name} '
#                 'StrandArray of the VirtualHelix',
#             )

#         idx_offsets: List[int] = strand_array.idx_offsets

#         seq = strand.seq
#         # oligo_old = strand.oligo

#         # 1. Do the 5' portion of the break
#         oligo_break5p = Oligo(seq[0:idx])
#         strand_break5p = oligo_break5p.strand5p

#         neighbor_strand_5p = strand.strand5p
#         if neighbor_strand_5p is not None:  # update existing neighbor oligos
#             strand_break5p.strand5p = neighbor_strand_5p
#             neighbor_strand_5p.strand3p = strand_break5p
#             for seg in neighbor_strand_5p.gen5p():
#                 seg.oligo = oligo_break5p

#         # 2. Do the 3' portion of the break
#         oligo_break3p = Oligo(seq[idx:])
#         strand_break3p = oligo_break3p.strand5p

#         neighbor_3p: Strand = strand.strand3p
#         if neighbor_3p is not None:  # update existing neighbor oligos
#             strand_break3p.strand3p = neighbor_3p
#             neighbor_3p.strand5p = strand_break3p
#             for seg in neighbor_3p.gen3p():
#                 seg.oligo = oligo_break3p

#         # 3. Update the strands

#         list_idx: int = vh_strands.index(strand)
#         offset_5p: int = idx_offsets[list_idx]
#         list_idx_plus_1: int = list_idx + 1

#         vh_strands.insert(
#             list_idx_plus_1,
#             strand_break3p,
#         )
#         vh_strands.insert(
#             list_idx_plus_1,
#             strand_break5p,
#         )

#         idx_offsets.insert(
#             list_idx_plus_1,
#             offset_5p + len(strand_break5p),
#         )

#         vh_strands.pop(list_idx)  # pop out the original strand

#         return oligo_break5p, oligo_break3p

#     def __add__(self, b: VirtualHelix) -> VirtualHelix:
#         '''(1) Concatenates the forward strand with forward strand
#         and the reverse strand with the reverse strand and preserves order
#         (2) Realigns the two :class:`VirtualHelix` objects involved

#         Args:
#             b: a :class:`VirtualHelix` object

#         Returns:
#             a :class:`VirtualHelix` object

#         Raises:
#             ValueError, TypeError
#         '''
#         if isinstance(b, VirtualHelix):
#             type3, seq3 = self.three_prime_end()
#             type5, seq5 = b.five_prime_end()
#             if type3 == type5 and len(seq3) == len(seq5):
#                 if seq3 != reverseComplement(seq5):
#                     raise TypeError('Ends not complimentary')
#                 fwd_strands = self.fwd_strands + b.fwd_strands
#                 rev_strands = self.rev_strands + b.rev_strands

#                 fwd_idx_offsets =self.strand_arrays[self.FWD].idx_offsets
#                 rev_idx_offsets =self.strand_arrays[self.REV].idx_offsets

#                 return VirtualHelix(
#                     fwd_strands=fwd_strands,
#                     fwd_idx_offsets=fwd_idx_offsets,
#                     rev_strands=rev_strands,
#                     rev_idx_offsets=rev_idx_offsets,
#                 )
#             else:
#                 raise TypeError('Ends not compatible')
#         else:
#             raise ValueError(f'{b} object not a DSeq')

#     def five_prime_end(self) -> Tuple[int, str]:
#         '''
#         Returns:
#             What kind of end is overhanging the 5' end of the
#             forward strand
#         '''
#         # fwd_idx0: int = self.fwd_idx_offsets[0]
#         # rev_idx0: int = self.rev_idx_offsets[-1]
#         res, val = five_prime_type(
#             self.alignment,
#             self.fwd_seq(),
#             self.rev_seq(),
#         )
#         return res, val

#     def three_prime_end(self) -> Tuple[int, str]:
#         '''
#         Returns:
#             What kind of end is overhanging the 3' end of the forward
#             strand
#         '''
#         res, val = three_prime_type(
#             self.alignment,
#             self.fwd_seq(),
#             self.rev_seq(),
#         )
#         return res, val


# def DSeqVH(
#         fwd: str,
#         rev: str = '',
#         overhang: Optional[int] = None,
#         alphabet: int = AlphaEnum.DNA,
# ) -> VirtualHelix:
#     '''Helper function for creating :class:`VirtualHelix` in the style of
#     the :class:`DSeq` with strings
#     '''
#     dseq: DSeq = DSeq(fwd, rev, overhang, alphabet)
#     # overhang: int = dseq.overhang
#     if dseq.overhang > 0:
#         fwd_idx_offsets = [dseq.overhang]
#         rev_idx_offsets = [0]
#     else:
#         fwd_idx_offsets = [0]
#         rev_idx_offsets = [dseq.overhang]
#     oligo_fwd = Oligo(fwd)
#     if rev is None:
#         rev = reverseComplement(fwd)
#     oligo_rev = Oligo(rev)
#     return VirtualHelix(
#         [oligo_fwd.strand5p],
#         fwd_idx_offsets,
#         [oligo_rev.strand5p],
#         rev_idx_offsets,
#     )


# if __name__ == '__main__':
#     fwd = 'GGTCTCGAATTCAAA'
#     oligo_fwd = Oligo(fwd)
#     rev = 'TTTGAATTCGAGACC'
#     oligo_rev = Oligo(rev)

#     BsaI_vh = VirtualHelix(
#         [oligo_fwd.strand5p], [0],
#         [oligo_rev.strand5p], [0],
#     )
#     print('1.\n%s' % BsaI_vh)
#     BsaI_vh = DSeqVH(fwd, rev, 0)
#     print('2.\n%s' % BsaI_vh)
#     print(BsaI_vh.fwd_strands)
#     BsaI_vh = DSeqVH(fwd)
#     print('3.\n%s' % BsaI_vh)
#     print('Da Oligos', BsaI_vh.oligos())
#     strand0 = BsaI_vh.fwd_strands[0]
#     print(strand0.oligo)
#     broken_oligos = BsaI_vh.breakStrand(dir_idx=0, strand=strand0, idx=4)
#     import pprint
#     pprint.pprint(broken_oligos)
#     print(BsaI_vh.len_fwd())

#     bonus_oligo = Oligo('ACGT')
#     try:
#         BsaI_vh.breakStrand(dir_idx=1, strand=bonus_oligo.strand5p, idx=2)
#     except ValueError:
#         print('Handled bad strand successfully')
