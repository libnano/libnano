import sys
from typing import (
    List,
    Generator
)
import uuid

import __init__
from libnano.seqrecord.seqrecordbase import (
    AlphaEnum,
    ALPHABETS
)

class Strand:
    def __init__(self,
            seq: str,
            oligo: 'Oligo',
            uid: str = None,
            ctx: dict = None):
        if uid is None:
            uid = uuid.uuid4()
        self.uid: str = uid
        self.seq: str = seq
        self.strand5p: Strand = None
        self.strand3p: Strand = None
        # if not isinstance(oligo, Oligo):
        #     raise TypeError('Every strand needs an Oligo')
        self.oligo: 'Oligo' = oligo
        if ctx is not None:
            ctx[uid] = self

    def __len__(self) -> int:
        return len(self.seq)

    def __str__(self) -> str:
        return self.seq

    def __repr__(self) -> str:
        return '%s: %s' % (self.uid, self.seq)

    def gen3p(self) -> Generator['Strand', None, None]:
        '''Iterate from self to the final `strand3p` of the :class:`Oligo` this
        :class:`Strand` is part of
        '''
        node0 = node = self
        while node is not None:
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
        while node is not None:
            yield node
            node = node.strand5p
            if node0 == node:
                break
    # end def
# end class

def UnknownStrand(length: int) -> Strand:
    seq: str = 'N'*length
    return Strand(seq)
# end def

if __name__ == '__main__':
    fwd = 'GGTCTCGAATTCAAA'
    strand = Strand(fwd, None)
    print(strand)