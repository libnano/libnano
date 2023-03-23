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
libnano._virtualhelix.strand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
'''
from __future__ import annotations

import uuid
from typing import (
    Dict,
    Iterator,
    Optional,
)

# from libnano._virtualhelix import oligo


class Strand:
    def __init__(
            self,
            seq: str,
            oligo,
            uid: str = '',
            ctx: Optional[Dict] = None,
    ):
        if not uid:
            uid = str(uuid.uuid4())
        self.uid: str = uid
        self.seq: str = seq
        self.strand5p: Optional[Strand] = None
        self.strand3p: Optional[Strand] = None
        # if not isinstance(oligo, Oligo):
        #     raise TypeError('Every strand needs an Oligo')
        self.oligo = oligo
        if ctx is not None:
            ctx[uid] = self

    def __len__(self) -> int:
        return len(self.seq)

    def __str__(self) -> str:
        return self.seq

    def __repr__(self) -> str:
        return f'{self.uid}: {self.seq}'

    def gen3p(self) -> Iterator[Strand]:
        '''Iterate from self to the final `strand3p` of the :class:`Oligo` this
        :class:`Strand` is part of
        '''
        node0 = node = self
        while node is not None:
            yield node
            node = node.strand3p  # type: ignore
            if node is None or node0 == node:
                break

    def gen5p(self) -> Iterator[Strand]:
        '''Iterate from self to the final `strand5p` of the :class:`Oligo` this
        :class:`Strand` is part of
        '''
        node0 = node = self
        while node is not None:
            yield node
            node = node.strand5p  # type: ignore
            if node is None or node0 == node:
                break


def UnknownStrand(length: int) -> Strand:
    seq: str = 'N' * length
    return Strand(seq, None)
