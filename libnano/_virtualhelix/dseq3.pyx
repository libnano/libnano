
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
from libnano.cymem.cymem cimport Pool


DEFAULT_SIZE: int = 256

cdef class OligoAssembly:
    cdef Pool mem   # Memory manager

    def __cinit__(self):
        self.mem = None

    def __init__(self, default_size: int = DEFAULT_SIZE):
        self.mem = Pool()
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

    def seq(self) -> str:
        '''Return the sequence string from 5' to 3'
        '''
        return ''.join(x.seq for x in strand5p.gen3p())

    def setSeq(self, seq: str):
        if len(seq) != len(self):
            raise ValueError("sequence length does not match oligo length")
        len_x: int
        for x in self.strand5p().gen3p():
            len_x = len(x)
            x.seq = seq[:len_x]
            seq = seq[len_x:]
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