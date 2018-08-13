from copy import (
    copy,
    deepcopy
)

from libnano.seqrecord.location import (
    Location,
    parseLocation,
    createLocationStr
)

from cpython.object cimport (
    Py_LT,
    Py_EQ,
    Py_NE
)

cdef class DummyFeature:
    cdef public int start
    cdef public int end
    def __cinit__(self, int start, int end):
        self.start = start
        self.end = end

    def __richcmp__(self, other, int op):
        cdef int start = self.start
        cdef int o_start = other.start
        if op == 0: # __lt__
            if start < o_start:
                return True
            elif start == o_start and self.end < other.end:
                return True
            else:
                return False
        else:
            raise NotImplementedError("Only __lt__ implemented")
# end def

def locationStr2Feature(feature_type, location_str, qualifiers=None):
    cdef bint is_fwd
    cdef object out_location

    location = parseLocation(location_str)

    if len(location) == 2:
        if location[0] == "complement" and isinstance(location[1], Location):
            out_location = location[1]
            is_fwd = False
        else:
            raise NotImplementedError("compound locations not supported")
    elif isinstance(location, Location):
        out_location = location
        is_fwd = True
    else:
        raise ValueError("need to pass this a location")

    return Feature(feature_type, out_location, is_fwd, qualifiers=qualifiers)

cdef class Feature:
    """
    """
    cdef readonly object feature_type
    cdef readonly bint is_fwd
    cdef readonly object location
    cdef readonly int start
    cdef readonly int end
    cdef readonly object qualifiers

    def __cinit__(self, feature_type, location, bint is_fwd, qualifiers=None):
        self.feature_type = feature_type
        self.is_fwd = is_fwd
        self.location = location
        self.start = location[1]
        self.end = location[2]
        self.qualifiers = qualifiers
    # end def

    def name(self):
        return self.feature_type

    def __copy__(self):
        return type(self)(copy(self.feature_type),
                            copy(self.location),
                            self.is_fwd,
                            deepcopy(self.qualifiers))
    # end def

    def copyForSlice(self, int offset):
        remote, idx5p, idx3p, flags = self.location
        loc = Location(remote, idx5p + offset, idx3p + offset, flags)
        return type(self)(copy(self.feature_type),
                            loc,
                            self.is_fwd,
                            deepcopy(self.qualifiers))

    def __richcmp__(self, other, int op):
        cdef int start = self.start
        cdef int o_start = other.start
        if op == Py_LT: # __lt__
            if start < o_start:
                return True
            elif start == o_start and self.end < other.end:
                return True
            else:
                return False
        elif op == Py_EQ: # __eq__
            return other.location == self.location and \
                    other.feature_type == self.feature_type and \
                    other.qualifiers == self.qualifiers
        elif op == Py_NE: # __ne__
            return other.location != self.location or \
                    other.feature_type != self.feature_type or \
                    other.qualifiers != self.qualifiers
        else:
            raise NotImplementedError("Only __lt__, __eq__, and __ne__ implemented")
    # end def

    def dump(self, bint clone):
        if self.is_fwd:
            loc = createLocationStr(self.location)
        else:
            loc = createLocationStr("complement", self.location)
        if clone:
            return { 'type': copy(self.feature_type),
                    'location': loc,
                    'qualifiers': deepcopy(self.qualifiers)
            }
        else:
            return { 'type': self.feature_type,
                    'location': loc,
                    'qualifiers': self.qualifiers
            }
    # end def


    def addQualifier(self, qualifier, value):
        self.qualifiers[qualifier] = value
    # end def

    def removeQualifier(self, qualifier):
        val = self.qualifiers[qualifier]
        del self.qualifiers[qualifier]
        return qualifier, val
    # end def