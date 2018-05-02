import json
import io
import os.path
from copy import (
    copy,
    deepcopy
)

from cpython.list cimport (
    PyList_New,
    PyList_SET_ITEM,
    PyList_Size
)
from cpython.ref cimport (
    PyObject,
    Py_XINCREF,
    Py_INCREF
)
from cpython.slice cimport (
    PySlice_GetIndicesEx,
    PySlice_Check
)

from libnano.fileio import (
    gb_reader,
    gb_writer
)
from libnano.fileio import fasta
from libnano.datastructures.seqrecord.feature import (
    Feature,
    DummyFeature,
    locationStr2Feature
)
from libnano.datastructures.seqrecord.featuretypes import FeatureTypes
from libnano.datastructures.list_bisect cimport (
    bisect_left,
    bisect_right,
    insort_left
)

def fromFasta(fasta_fn, alphabet='DNA', not_allowed=None):
    recs = fasta.parseFasta(fasta_fn, alphabet, not_allowed)
    out = []
    for rec in recs:
        sr = SeqRecord()
        sr.name = rec[0]
        sr.seq = rec[1]
        out.append(sr)
    return out
# end def

def fromGenbankLike(fn_string, feature_types=None):
    sr = SeqRecord(feature_types=feature_types)
    sr.fromGenbankLike(fn_string)
    return sr
# end def

cdef class SeqRecord:
    """
    def __getitem__(self, slice):
            return self.data[slice]

    slice.start slice.step slice.stop
    """
    cdef public object info
    cdef object sequence, feature_types
    cdef dict feature_instance_dict
    cdef list feature_instance_list

    def __init__(self, feature_types=None):
        cdef object obj
        self.feature_instance_dict = {}
        self.feature_instance_list = []

        if feature_types is None:
            self.feature_types = FeatureTypes()
        else:
            self.feature_types = feature_types

        self.sequence = ''
        self.info = {}
    # end def

    property name:
        def __get__(self):
            return self.info['name']
        def __set__(self, value):
            self.info['name'] = value
        def __del__(self):
            del self.info['name']

    property description:
        def __get__(self):
            return self.info['definition']
        def __set__(self, value):
            self.info['definition'] = value
        def __del__(self):
            del self.info['definition']

    property ID:
        def __get__(self):
            return self.info['version']
        def __set__(self, value):
            self.info['version'] = value
        def __del__(self):
            del self.info['version']

    property seq:
        def __get__(self):
            return self.sequence
        def __set__(self, value):
            self.sequence = value

    def __contains__(self, feature_type):
        """ allow for
            feature_type in SeqRecord()
        construction
        """
        return feature_type in self.feature_instance_dict

    def __getitem__(self, in_slice):
        """ slicing a SeqRecord.  Does not support step in slice
        but does support negative indexing

        information is dropped
        """
        cdef Py_ssize_t start, stop, step, slicelength
        cdef Py_ssize_t seq_length = len(self.sequence)
        cdef list feature_list
        cdef int i
        if PySlice_Check(in_slice):
            PySlice_GetIndicesEx(in_slice, seq_length, &start, &stop, &step, &slicelength)
        else:
            start = in_slice
            if seq_length < start:
                raise IndexError("Index {} out of bounds".format(in_slice))
            else:
                if start < 0:
                    start = seq_length + start
                stop = start + 1

        sr = type(self)(feature_types=self.feature_types)
        feature_list = self.getFeatures(start, stop)
        for i in range(len(feature_list)):
            f = feature_list[i]
            fcopy = f.copyForSlice(start)
            sr.addFeature(f)
        sr.seq = self.sequence[start:stop]
        return sr
    # end def

    def __copy__(self):
        """ don't copy info
        """
        sr = type(self)(feature_types=copy(self.feature_types))
        sr.seq = copy(self.sequence)
        for f in self.feature_instance_list:
            sr.addFeature(f.__copy__()) # copy the feature too?
        return sr

    def __add__(self, x):
        """ don't copy info
        """
        cdef Py_ssize_t start = len(self.sequence)
        if not isinstance(x, SeqRecord):
            raise ValueError("Can't concatenate non-SeqRecord of type {}", type(x))
        # unsure about copy of feature_types
        sr = self.__copy__()
        for f in x.iterFeatures():
            fcopy = f.copyForSlice(start)
            sr.addFeature(fcopy)
        sr.seq = sr.seq + x.seq
    # end def

    cdef constructFeaturesFromGenbankLike(self, list feature_list):
        """ Does not assume input feature_list is sorted
        Could speed this up with that assumption
        """
        cdef Py_ssize_t i
        cdef Py_ssize_t len_features = PyList_Size(feature_list)
        for i in range(len_features):
            fobj = feature_list[i]
            f = locationStr2Feature(fobj['type'],
                                    fobj['location'],
                                    fobj['qualifiers'])
            self.addFeature(f)
        # end for
    # end def

    def fromGenbankLike(self, fn_string):
        extension = os.path.splitext(fn_string)[1]
        if extension == '.json':
            with io.open(fn_string, 'r', encoding='utf-8') as fd_json:
                obj = json.load(fd_json)
        elif extension == '.gb':
            obj = gb_reader.parse(fn_string)
        else:
            raise NotImplementedError("Only json and gb are supported")
        self.seq = obj['seq']
        self.constructFeaturesFromGenbankLike(obj['features'])
        self.info = obj['info']
    # end def

    cdef serializeFeatures(self, clone):
        cdef object item
        cdef Py_ssize_t i
        cdef list feature_instance_list = self.feature_instance_list
        cdef Py_ssize_t len_features = PyList_Size(feature_instance_list)
        cdef list out = PyList_New(len_features)    # preallocate known size

        for i in range(len_features):
            item = feature_instance_list[i].dump(clone)
            Py_INCREF(item) # cython needs this due to SET_ITEM ref steal
            PyList_SET_ITEM(out, i, item)
        # end for
        return out
    # end def

    def dump(self, clone=True, to_json=False, to_json_file=None):
        """
        Create a dictionary of the SeqRecord
        """
        d = { 'seq': self.sequence,
                'features': self.serializeFeatures(clone),
                'info': self.info
            }
        if to_json:
            return json.dumps(d)
        elif to_json_file is not None:
            with io.open(to_json_file, 'w', encoding='utf-8') as fd:
                fd.write(json.dumps(d, indent=4))
            return None
        else:
            return d
    # end def

    def toGenBank(self, fn_string, order_qualifiers=False):
        return gb_writer.write(fn_string,
                                self.dump(clone=False),
                                order_qualifiers=order_qualifiers)

    def toFasta(self, fn_string):
        description = self.info['name']
        return fasta.write(fn_string, [(description, self.sequence)])

    cpdef addFeature(self, feature, description=None):
        location = feature.location
        feature_name = feature.name()
        feature_types = self.feature_types

        feature_instance_list = self.feature_instance_list
        feature_instance_dict = self.feature_instance_dict

        if feature_name not in feature_types:
            ft_id = feature_types.addFeatureType(feature_name, description=description)
        else:
            ft_id = feature_types.getFTID(feature_name)

        insort_left(self.feature_instance_list, feature, 0, -1)

        if feature_name in self.feature_instance_dict:
            insort_left(feature_instance_dict[feature_name], feature, 0, -1)
        else:
            feature_instance_dict[feature_name] = [feature]
    # end def

    cpdef removeFeature(self, feature):
        cdef Py_ssize_t idx
        feature_name = feature.name()

        feature_instance_list = self.feature_instance_list
        feature_instance_dict = self.feature_instance_dict

        idx = bisect_left(feature_instance_list, feature, 0, -1)
        if idx == len(feature_instance_list) or \
            feature_instance_list[idx] != feature:
            raise ValueError("removeFeature: feature not in SeqRecord")
        feature_instance_list.pop(idx)


        ilist = feature_instance_dict[feature_name]
        ilist.remove(feature)
        if len(ilist) == 0:
            del feature_instance_dict[feature_name]
    # end def

    cpdef list getFeatures(self, Py_ssize_t start, Py_ssize_t end):
        cdef Py_ssize_t idx_lo, idx_hi
        cdef object feature_instance_list = self.feature_instance_list
        cdef object temp = DummyFeature(start, start)
        idx_lo = bisect_left(feature_instance_list, temp, 0 , -1)
        temp = DummyFeature(end, end)
        idx_hi = bisect_right(feature_instance_list, temp, 0, -1)
        return feature_instance_list[idx_lo:idx_hi]
    # end def

    def iterFeatures(self):
        """ iterator and not the list itself since we must not
        expose the feature_instance_list object itself
        """
        return iter(self.feature_instance_list)
# end class

