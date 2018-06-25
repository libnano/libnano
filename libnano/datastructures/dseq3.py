
from typing import (
    List,
    Tuple,
    Generator
)
import uuid
from heapq import (
    heapify,
    heappush,
    nsmallest
)
import numpy as np
import pandas as pd

DEFAULT_SIZE: int = 256

class Strand:
    def __init__(self, seq: str, uid: str = None):
        self.uid: str = uuid.uuid4() uid is None else uid
        self.seq: str = seq
        self.strand5p: Strand = None
        self.strand3p: Strand = None
        self.oligo: Oligo = None

    def gen3p(self) -> Generator['Strand', None, None]:
        '''Iterate from self to the final `strand3p` of the :class:`Oligo` this
        :class:`Strand` is part of
        '''
        node0 = node = self
        while node in not None:
            yield node
            node = node.strand3p
            if node0 == node:
                break
    # end def

    def gen5p(self) -> Generator['Strand', None, None]:
        '''Iterate from self to the final `strand5p` of the :class:`Oligo` this
        :class:`Strand` is part of
        '''
        node0 = node = self
        while node in not None:
            yield node
            node = node.strand5p
            if node0 == node:
                break
    # end def

class Oligo:
    def __init__(self, strand5p: Strand = None):
        self.strand5p: Strand = strand5p

    def append3p(self, strand: Strand):
        '''Append a :class:`Strand` to the 3 prime end of the oligo
        '''
        strand5p: Strand = self.strand5p
        strand3p: Strand = strand5p.strand3p
        while strand3p is not None:
            strand5p = strand3p
            strand3p = strand5p.strand3p
        strand5p.strand3p = strand
    # end def

    def prepend5p(self, strand: Strand):
        '''Append a :class:`Strand` to the 3 prime end of the oligo
        '''
        old_strand5p: Strand = self.strand5p
        assert(old_strand5p.strand5p is None)
        old_strand5p.strand5p = strand
        self.strand5p = strand
        strand.strand3p = old_strand5p
    # end def

    def seq(self) -> str:
        '''Return the sequence string from 5' to 3'
        '''
        return ''.join(strand.seq for strand in strand5p.gen3p())
# end class

class OligoAssembly:
    def __init__(self, default_size: int = DEFAULT_SIZE):
        self.default_size: int = default_size
        self.oligos: List[Oligo] = []
        self.strands: List[str] = []
        self.id_nums: np.ndarray = np.full((DEFAULT_SIZE,), -1, dtype=int)
        self.fwd_strandsets: List[Strand] = [None] * default_size
        self.rev_strandsets: List[Strand] = [None] * default_size

        self.highest_id_num_used: int = -1
        self.recycle_bin: List[int] = []
        self.reserved_ids: set = set()
    # end def

    def getNewIdNum(self): -> int:
        '''Query the lowest available (unused) id_num. Internally id numbers are
        recycled when virtual helices are deleted.
        Returns:
            int: ID number
        '''
        if len(self.recycle_bin):
            return nsmallest(1, self.recycle_bin)[0]
        else:
            # use self.highest_id_num_used if the recycle bin is empty
            # and highest_id_num_used + 1 is not in the reserve bin
            return self.highest_id_num_used + 1
    # end def

    def reserveIdNum(self, requested_id_num: int):
        '''Reserves and returns a unique numerical id_num appropriate for a
        virtualhelix of a given parity. If a specific index is preferable
        (say, for undo/redo) it can be requested in num.
        Args:
            requested_id_num (int): virtual helix ID number
        '''
        num: int = requested_id_num
        assert num >= 0, int(num) == num
        # assert not num in self._number_to_virtual_helix
        if num in self.recycle_bin:
            self.recycle_bin.remove(num)
            # rebuild the heap since we removed a specific item
            heapify(self.recycle_bin)
        if self.highest_id_num_used < num:
            self.highest_id_num_used = num
        self.reserved_ids.add(num)
    # end def

    def recycleIdNum(self, id_num: int):
        '''The caller's contract is to ensure that id_num is not used in *any*
        helix at the time of the calling of this function (or afterwards, unless
        `reserveIdNumForHelix` returns the id_num again).
        Args:
            id_num (int): virtual helix ID number
        '''
        heappush(self.recycle_bin, id_num)
        self.reserved_ids.remove(id_num)
    # end def

    def getIdNumMax(self) -> int:
        '''The max id number
        Returns:
            int: max virtual helix ID number used
        '''
        return self.highest_id_num_used
    # end def
# end class

class Helix:
    def __init__(self,
        fwd_strands: List[str],
        fwd_idx_offsets: List[int]
        rev_strands: List[str]
        rev_idx_offsets: List[int],
        alphabet: int = AlphaEnum.DNA):

        self.fwd_strands = fwd_strands
        fwd_str: str = ''.join(' '*x + y for x, y in zip(fwd_strands, fwd_idx_offsets))
        self.rev_strands = rev_strands
        rev_str: str = ''.join(' '*x + y for x, y in zip(rev_strands, rev_idx_offsets))

        self.alphabet: int = alphabet
    # end def
# end class

def HelixHelper(fwd: str, rev: str, overhang: int) -> Helix:
    if overhang > 0:
        fwd_idx_offsets = [overhang]
        rev_idx_offsets = []
    else:
        fwd_idx_offsets = []
        rev_idx_offsets = [overhang]
    return Helix([fwd], fwd_idx_offsets, [rev], rev_idx_offsets)
# end def