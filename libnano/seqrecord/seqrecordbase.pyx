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
libnano.seqrecord.seqrecordbase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
from __future__ import annotations

import io
import json
import os.path
from copy import copy  # deepcopy,
from enum import IntEnum
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
    Union,
)

from cpython.list cimport (  # type: ignore
    PyList_New,
    PyList_SET_ITEM,
    PyList_Size,
)
from cpython.ref cimport Py_INCREF  # type: ignore; Py_XINCREF,; PyObject,
from cpython.slice cimport (  # type: ignore
    PySlice_Check,
    PySlice_GetIndicesEx,
)

from libnano.fileio import (
    fasta,
    gb_reader,
    gb_writer,
)
from libnano.seqrecord.feature import (  # type: ignore
    DummyFeature,
    Feature,
    locationStr2Feature,
)
from libnano.seqrecord.featuretypes import FeatureTypes

from libnano.list_bisect cimport (
    bisect_left,
    bisect_right,
    insort_left,
)


class AlphaEnum(IntEnum):
    DNA: int = 0
    RNA: int = 1

ALPHABETS: Tuple[int, int] = tuple(a.value for a in AlphaEnum)

cdef class SeqRecord:
    '''SeqRecord

    Attributes:
        info: Information dictionary
    '''
    cdef:
        public object info
        object sequence
        object feature_types
        dict feature_instance_dict
        list feature_instance_list

    def __init__(
            self,
            sequence: str = '',
            feature_types: Optional[FeatureTypes ]= None,
            name: str = '',
            alphabet: int = AlphaEnum.DNA,
    ):
        '''
        Args:
            sequence: Optional sequence string
            feature_types: Optional :class:`FeatureTypes`
            name: Optional name of the record,
            alphabet: {AlphaEnum.DNA, AlphaEnum.RNA}
        '''
        cdef object obj
        self.feature_instance_dict: Dict[str, List[Feature]] = {}
        self.feature_instance_list: List[Feature] = []

        if feature_types is None:
            self.feature_types = FeatureTypes()
        else:
            self.feature_types = feature_types

        assert(alphabet in ALPHABETS)
        if name is None:
            name = ''
        assert(isinstance(name, str))

        self.info: Dict[str, Any] = {
            'alphabet': alphabet,
            'name': name
        }

        if sequence:
            # TODO:add alphabet linting for DNA and/or RNA
            self.sequence = sequence.lower()
        else:
            self.sequence = ''

    @property
    def name(self) -> str:
        return self.info['name']

    @name.setter
    def name(self, value: str):
        self.info['name'] = value

    @name.deleter
    def name(self):
        del self.info['name']

    @property
    def description(self) -> str:
        return self.info.get('definition', '')

    @description.setter
    def description(self, value: str):
        self.info['definition'] = value

    @description.deleter
    def description(self):
        del self.info['definition']

    @property
    def ID(self) -> str:
        return self.info.get('version', '')

    @ID.setter
    def ID(self, value: str):
        self.info['version'] = value

    @ID.deleter
    def ID(self):
        del self.info['version']

    @property
    def alphabet(self) -> int:
        return self.info['alphabet']


    @property
    def seq(self) -> str:
        return self.sequence

    @seq.setter
    def seq(self, value: str):
        self.sequence = value.lower()

    def __contains__(
            self,
            feature_type: Feature,
    ) -> bool:
        '''Allow for ``feature_type in SeqRecord()`` construction
        '''
        return feature_type in self.feature_instance_dict

    def __getitem__(self, in_slice):
        '''Slicing a :class:`SeqRecord`.  Does not support step in slice
        but does support negative indexing

        Information is dropped


        '''
        cdef:
            Py_ssize_t start
            Py_ssize_t stop
            Py_ssize_t step
            Py_ssize_t slicelength
            Py_ssize_t seq_length = len(self.sequence)
            list feature_list
            int i

        if PySlice_Check(in_slice):
            PySlice_GetIndicesEx(
                in_slice,
                seq_length,
                &start,
                &stop,
                &step,
                &slicelength,
            )
        else:
            start = in_slice
            if seq_length < start:
                raise IndexError(f'Index {in_slice} out of bounds')
            else:
                if start < 0:
                    start = seq_length + start
                stop = start + 1

        sr = type(self)(
            feature_types=self.feature_types,
        )
        feature_list = self.getFeatures(
            start,
            stop,
        )
        for i in range(len(feature_list)):
            f = feature_list[i]
            fcopy = f.copyForSlice(start)
            sr.addFeature(fcopy)
        sr.seq = self.sequence[start:stop]
        return sr

    def __copy__(self):
        '''Do not copy info
        '''
        sr = type(self)(
            feature_types=copy(self.feature_types),
            name=self.name,
            alphabet=self.alphabet,
            sequence=copy(self.sequence),
        )
        for f in self.feature_instance_list:
            sr.addFeature(f.__copy__()) # copy the feature too?
        return sr

    def __add__(self, x: SeqRecord):
        '''Do not copy info

        Raises:
            ValueError: Cannot concatenate non-SeqRecord of type
        '''
        cdef Py_ssize_t start = len(self.sequence)

        if not isinstance(x, SeqRecord):
            raise ValueError(
                f'Cannot concatenate non-SeqRecord of type {type(x)}'
            )
        # Unsure about copy of feature_types
        sr = self.__copy__()
        for f in x.iterFeatures():
            fcopy = f.copyForSlice(start)
            sr.addFeature(fcopy)
        sr.seq = sr.seq + x.seq
        return sr

    def __eq__(self, x: SeqRecord):
        return True if self.seq == x.seq else False

    def __hash__(self):
        return hash(self.seq.lower())

    cdef constructFeaturesFromGenbankLike(
            self,
            list feature_list,
    ):
        '''Add :class:`Feautures` from a list of Genbank like items.
        Does not assume input feature_list is sorted
        Could speed this up with that assumption

        Args:
            feature_list: List of :class:`Features`
        '''
        cdef:
            Py_ssize_t i
            Py_ssize_t len_features = PyList_Size(feature_list)

        for i in range(len_features):
            fobj = feature_list[i]
            f = locationStr2Feature(
                fobj['type'],
                fobj['location'],
                fobj['qualifiers'],
            )
            self.addFeature(f)

    def fromGenbankLike(
            self,
            fn_string: str,
    ) -> None:
        '''Create record from Genbank-like string

        Args:
            fn_string: string

        Raises:
            NotImplementedError: Only json and gb are supported
        '''
        extension = os.path.splitext(fn_string)[1]
        if extension == '.json':
            with io.open(fn_string, 'r', encoding='utf-8') as fd_json:
                obj = json.load(fd_json)
        elif extension == '.gb':
            obj = gb_reader.parse(fn_string)
        else:
            raise NotImplementedError('Only json and gb are supported')

        self.seq = obj['seq']
        self.constructFeaturesFromGenbankLike(obj['features'])
        self.info = obj['info']

    cdef list serializeFeatures(
            self,
            bint clone,
    ):
        '''
        Args:
            clone: If True, clone the out put

        Returns:
            List of serialized Features
        '''
        cdef:
            object item
            Py_ssize_t i
            list fi_list = self.feature_instance_list
            Py_ssize_t len_features = PyList_Size(fi_list)
            list out = PyList_New(len_features)    # Preallocate known size

        for i in range(len_features):
            item = fi_list[i].dump(clone)
            Py_INCREF(item) # Cython needs this due to SET_ITEM ref steal
            PyList_SET_ITEM(
                out,
                i,
                item,
            )
        return out

    def dump(
            self,
            clone: bool = True,
            to_json: bool = False,
            to_json_file: str = None,
    ) -> Union[str, Dict]:
        '''Create a dictionary or string of the :class:`SeqRecord`

        Args:
            clone: whether to clone the features or not
            to_json: True returns a string.  False returns a dictionary
            to_json_file: whether to write to a file

        Returns:
            dictionary or a string
        '''
        d = {
            'seq': self.sequence,
            'features': self.serializeFeatures(clone),
            'info': self.info
        }
        if to_json:
            return json.dumps(d)
        elif to_json_file is not None:
            with io.open(to_json_file, 'w', encoding='utf-8') as fd:
                fd.write(json.dumps(d, indent=4))
        return d

    def toGenBank(
            self,
            fn_string: str,
            order_qualifiers: bool = False,
    ) -> None:
        '''Write record to a Genbank file

        Args:
            fn_string: filepath string
            order_qualifiers: If True, write features in order
        '''
        with open(fn_string, 'w') as fd:
            gb_writer.write(
                fd,
                self.dump(clone=False),
                order_qualifiers=order_qualifiers,
            )

    def toFasta(
            self,
            fn_string: str,
    ) -> None:
        '''Write record to a FASTA file

        Args:
            fn_string: filepath string
        '''
        description = self.info['name']

        fasta.write(
            fn_string,
            [(description, self.sequence)],
        )

    cpdef addFeature(
            self,
            feature: Feature,
            description: str = None,
    ):
        '''Add :class:`Feature` to the record

        Args:
            feature: :class:`Feature` instance to add
            description: Desctription of the feature
        '''
        # location = feature.location
        feature_name = feature.name()
        # feature_types = self.feature_types

        fi_list: List[Feature] = self.feature_instance_list
        fi_dict: Dict[str, List[Feature]] = self.feature_instance_dict

        # if feature_name not in feature_types:
        #     ft_id = feature_types.addFeatureType(
        #         feature_name,
        #         description=description,
        #     )
        # else:
        #     ft_id = feature_types.getFTID(
        #         feature_name,
        #     )

        insort_left(
            fi_list,
            feature,
            0,
            -1,
        )

        if feature_name in fi_dict:
            insort_left(
                fi_dict[feature_name],
                feature,
                0,
                -1,
            )
        else:
            fi_dict[feature_name] = [feature]

    cpdef removeFeature(
            self,
            feature: Feature,
    ):
        '''Remove :class:`Feature` from the record

        Args:
            feature: :class:`Feature` instance to add

        Raises:
            ValueError: removeFeature: Feature not in SeqRecord
        '''
        cdef Py_ssize_t idx
        feature_name = feature.name()

        fi_list: List[Feature] =              self.feature_instance_list
        fi_dict: Dict[str, List[Feature]] =   self.feature_instance_dict

        idx = bisect_left(fi_list, feature, 0, -1)
        if idx == len(fi_list) or \
            fi_list[idx] != feature:
            raise ValueError(
                'removeFeature: Feature not in SeqRecord'
            )
        fi_list.pop(idx)

        ilist = fi_dict[feature_name]
        ilist.remove(feature)
        if len(ilist) == 0:
            del fi_dict[feature_name]

    cpdef list getFeatures(
            self,
            Py_ssize_t start,
            Py_ssize_t end,
    ):
        '''Get :class:`Feature`s between indices

        Args:
            start: Start index
            end: End index

        Returns:
            List of :class:`Features` between the indices

        '''
        cdef:
            Py_ssize_t idx_lo, idx_hi
            object fi_list = self.feature_instance_list
            object temp = DummyFeature(start, start)

        idx_lo = bisect_left(
            fi_list,
            temp,
            0,
            -1,
        )
        temp = DummyFeature(
            end,
            end,
        )
        idx_hi = bisect_right(
            fi_list,
            temp,
            0,
            -1,
        )
        return fi_list[idx_lo:idx_hi]

    def iterFeatures(self) -> Iterator[Feature]:
        '''Iterator and not the list itself since we must not
        expose the feature_instance_list object itself

        Returns:
            Iterator over all :class:`Features` in the list

        '''
        return iter(self.feature_instance_list)


def fromFasta(
        fasta_fn: str,
        alphabet: str = 'DNA',
        not_allowed: str = None,
) -> List[SeqRecord]:
    '''Generate a list of :class:`SeqRecord`s from a FASTA file

    Args:
        fasta_fn: FASTA file path
        alphabet: {'DNA', 'RNA', 'AMINO_ACID'}
        not_allowed: String of characters not allow in the sequence

    Returns:
        List of :class:`SeqRecord`s

    '''
    recs = fasta.parseFasta(
        fasta_fn,
        alphabet,
        not_allowed,
    )
    out = []
    for rec in recs:
        sr = SeqRecord()
        sr.name = rec[0]
        sr.seq = rec[1]
        out.append(sr)
    return out


def fromGenbankLike(
        fn_string: str,
        feature_types: Optional[FeatureTypes] = None,
) -> SeqRecord:
    '''Create :class:`SeqRecord`s from a Genbank file

    Args:
        fn_string: Genbank file path
        not_allowed: String of characters not allow in the sequence

    Returns:
        :class:`SeqRecord` instance

    '''
    sr = SeqRecord(
        feature_types=feature_types,
    )
    sr.fromGenbankLike(fn_string)
    return sr
