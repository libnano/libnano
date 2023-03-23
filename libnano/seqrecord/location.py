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
libnano.seqrecord.location
~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
from __future__ import annotations

from typing import (
    List,
    NamedTuple,
    Tuple,
    Union,
)

# Compound_location = (operator, (remote, idx5p, idx3p, flags))
# flags = <IS_PARTIAL_5 | IS_PARTIAL_3 | IS_INTERNAL]
IS_INTERNAL = 1
IS_PARTIAL_3 = 2
IS_PARTIAL_5 = 4


# base 'simple' location
class Location(NamedTuple):
    remote: str
    idx5p: int
    idx3p: int
    flags: int


class CompoundLocation(NamedTuple):
    operator: str
    location_like: Union[
        List[Union[CompoundLocation, Location]],
        CompoundLocation,
        Location,
    ]


ULocation_T = Union[
    List[Union[CompoundLocation, Location]],
    CompoundLocation,
    Location,
]


def createLocationStr(
        location_like: ULocation_T,
) -> str:
    '''
    Args:
        location_like: :class:`Location` instance or tuple of form:

            <operator>, <CompoundLocation instance>

    Returns:
        Location string
    '''
    if isinstance(location_like, CompoundLocation):
        return (
            f'{location_like.operator}'
            f'({createLocationStr(location_like.location_like)})'
        )
    elif isinstance(location_like, list):
        out_list = []
        for item in location_like:
            out_list.append(createLocationStr(item))
        return ','.join(out_list)
    else:
        return createLocationStrBases(
            location_like,
        )


def createLocationStrBases(
        loc: Location,
) -> str:
    '''
    flags = <IS_PARTIAL_5 | IS_PARTIAL_3 | IS_INTERNAL>

    Args:
        location: :class:`Location` instance

    Returns:
        string form of location
    '''
    out_list = []
    if loc.remote:
        out_list += [loc.remote, ':']
    if loc.flags & IS_PARTIAL_5:
        out_list.append('<')
    out_list.append(str(loc.idx5p + 1))
    if loc.idx3p != loc.idx5p:  # is not None:
        if loc.flags & IS_INTERNAL:
            out_list.append('^')
        else:
            out_list.append('..')
            if loc.flags & IS_PARTIAL_3:
                out_list.append('>')
        out_list.append(str(loc.idx3p + 1))
    return ''.join(out_list)


def parseLocation(
        loc_str: str,
) -> Union[CompoundLocation, Location]:
    '''Assumes last character is a matching parenthesis

    Args:
        loc_str: Location string form

    Returns:
        :class:`CompoundLocation` instance or :class:`Location` instance
    '''
    if loc_str[0:5] == 'join(':
        return CompoundLocation(
            'join',
            parseJoinOrder(loc_str[5:-1]),
        )
    elif loc_str[0:6] == 'order(':
        return CompoundLocation(
            'order',
            parseJoinOrder(loc_str[6:-1]),
        )
    elif loc_str[0:11] == 'complement(':
        return CompoundLocation(
            'complement',
            parseLocation(loc_str[11:-1]),
        )
    else:
        return parseLocationBases(loc_str)


def parseJoinOrder(
        loc_str: str,
) -> List[Union[CompoundLocation, Location]]:
    '''Find outer most commas

    Args:
        loc_str: Location string form

    Returns:
        List of tuples of the form:

            <remote string>, <parsed location string>
    '''
    comma_idx_list = []
    ignore_flag = 0
    for i, c in enumerate(loc_str):
        if c == ',' and ignore_flag == 0:
            comma_idx_list.append(i)
        elif c == '(':
            ignore_flag += 1
        elif c == ')':
            ignore_flag -= 1
    j = 0
    out_list = []
    for i in comma_idx_list:
        out_list.append(
            parseLocation(loc_str[j:i]),
        )
        j = i + 1
    out_list.append(
        parseLocation(loc_str[j:]),
    )
    return out_list


def parseLocationBases(
        bases_str: str,
) -> Location:
    '''Find remote if it exists

    Args:
        bases_str: String form of bases

    Returns:
        :class:`Location` instance of ``bases_str``
    '''
    if ':' in bases_str:
        bases_list = bases_str.split(':')
        remote = bases_list[0]
        bases_str = bases_list[1]
    else:
        remote = ''

    if '^' in bases_str:
        bl = bases_str.split('^')
        return Location(
            remote,
            int(bl[0]) - 1,
            int(bl[1]) - 1,
            IS_INTERNAL,
        )
    elif '..' in bases_str:
        flags = 0
        bl = bases_str.split('..')
        if bl[0][0] == '<':
            bl[0] = bl[0][1:]
            flags += IS_PARTIAL_5
        if bl[1][0] == '<':
            bl[1] = bl[1][1:]
            flags += IS_PARTIAL_3
        return Location(
            remote,
            int(bl[0]) - 1,
            int(bl[1]) - 1,
            flags,
        )
    else:  # Single base
        base_idx = int(bases_str) - 1
        return Location(
            remote,
            base_idx,
            base_idx,
            0,
        )
