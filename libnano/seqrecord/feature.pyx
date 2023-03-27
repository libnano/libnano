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
libnano.seqrecord.feature
~~~~~~~~~~~~~~~~~~~~~~~

'''
from copy import (
    copy,
    deepcopy,
)
from typing import (
    Any,
    Dict,
    Optional,
)

from libnano.seqrecord.location import (
    CompoundLocation,
    Location,
    createLocationStr,
    parseLocation,
)

from cpython.object cimport (  # type ignore
    Py_EQ,
    Py_LT,
    Py_NE,
)


cdef class DummyFeature:
    cdef:
        public int start
        public int end

    def __cinit__(
            self,
            start: int,
            end: int,
    ):
        '''
        Args:
            start: Start index
            end: End index
        '''
        self.start = start
        self.end = end

    def __richcmp__(self, other, int op):
        '''
        Raises:
            NotImplementedError: Only __lt__ implemented
        '''
        cdef:
            int start = self.start
            int o_start = other.start

        if op == 0: # __lt__
            if start < o_start:
                return True
            elif start == o_start and self.end < other.end:
                return True
            else:
                return False
        else:
            raise NotImplementedError('Only __lt__ implemented')


def locationStr2Feature(
        feature_type: str,
        location_str: str,
        qualifiers: Optional[Dict] = None,
) -> Feature:
    '''
    Args:
        feature_type: Type or name of feature
        location_str: Location string
        qualifiers: Optional qualifiers dictionary

    Returns:
        :class:`Feature` instance

    Raises:
        NotImplementedError: Compound locations not supported
        ValueError: Need to pass this a location
    '''
    cdef:
        bint is_fwd
        object out_location

    location = parseLocation(
        location_str,
    )

    if len(location) == 2:
        if (
            location[0] == 'complement' and
            isinstance(location[1], Location)
        ):
            out_location = location[1]
            is_fwd = False
        else:
            raise NotImplementedError('Compound locations not supported')
    elif isinstance(location, Location):
        out_location = location
        is_fwd = True
    else:
        raise ValueError('Need to pass this a location')

    return Feature(
        feature_type,
        out_location,
        is_fwd,
        qualifiers=qualifiers,
    )

cdef class Feature:
    '''
    Attributes:
        feature_type: Type or name of feature:
        is_fwd: If True, is a forward strand
        location: :class:`Location` of the feature
        start: Start index
        end: End index
        qualifiers: Qualifiers dictionary
    '''
    cdef:
        readonly object feature_type
        readonly bint is_fwd
        readonly object location
        readonly int start
        readonly int end
        readonly object qualifiers

    def __cinit__(
            self,
            feature_type: str,
            location: Location,
            bint is_fwd,
            qualifiers=None,
    ):
        '''
        Args:
            feature_type: Type or name of feature
            location: Location
            is_fwd: If True, is a forward strand
            qualifiers: Optional qualifiers dictionary
        '''
        self.feature_type = feature_type
        self.is_fwd = is_fwd
        self.location = location
        self.start = location[1]
        self.end = location[2]
        self.qualifiers = qualifiers

    def name(self) -> str:
        '''
        Returns:
            :attr:`feature_type` string
        '''
        return self.feature_type

    def __copy__(self):
        return type(self)(
            copy(self.feature_type),
            copy(self.location),
            self.is_fwd,
            deepcopy(self.qualifiers),
        )

    def copyForSlice(
            self,
            offset: int,
     ) -> Feature:
        remote, idx5p, idx3p, flags = self.location
        loc = Location(
            remote,
            idx5p + offset,
            idx3p + offset,
            flags,
        )
        return type(self)(
            copy(self.feature_type),
            loc,
            self.is_fwd,
            deepcopy(self.qualifiers),
        )

    def __richcmp__(self, other, op: int) -> bool:
        '''

        Raises:
            NotImplementedError: Only __lt__, __eq__, and __ne__ implemented
        '''
        cdef:
            int start = self.start
            int o_start = other.start

        if op == Py_LT: # __lt__
            if start < o_start:
                return True
            elif (
                start == o_start and
                self.end < other.end
            ):
                return True
            else:
                return False
        elif op == Py_EQ: # __eq__
            return (
                other.location == self.location and
                other.feature_type == self.feature_type and
                other.qualifiers == self.qualifiers
            )
        elif op == Py_NE: # __ne__
            return (
                other.location != self.location or
                other.feature_type != self.feature_type or
                other.qualifiers != self.qualifiers
            )
        else:
            raise NotImplementedError(
                'Only __lt__, __eq__, and __ne__ implemented'
            )

    def dump(self, bint clone):
        '''
        Args:
            clone: If True, clone the feature
        '''
        if self.is_fwd:
            loc_str = createLocationStr(
                CompoundLocation('', self.location)
            )
        else:
            loc_str = createLocationStr(
                CompoundLocation('complement', self.location),
            )
        if clone:
            return {
                'type': copy(self.feature_type),
                'location': loc_str,
                'qualifiers': deepcopy(self.qualifiers),
            }
        else:
            return {
                'type': self.feature_type,
                'location': loc_str,
                'qualifiers': self.qualifiers,
            }

    def addQualifier(
            self,
            qualifier: str,
            value: Any,
    ) -> None:
        '''
        Args:
            qualifier: Qualifier name
            value: Qualifier value
        '''
        self.qualifiers[qualifier] = value

    def removeQualifier(self, qualifier):
        val = self.qualifiers[qualifier]
        del self.qualifiers[qualifier]
        return qualifier, val
