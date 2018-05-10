from cpython.ref cimport PyObject
from cpython.list cimport (
    PyList_CheckExact,
    PyList_Insert,
    PyList_GET_ITEM
)
from cpython.list cimport (
    PyList_GetItem,
    PyList_Size
)
from cpython.object cimport (
    PyObject_RichCompareBool,
    Py_LT
)

cdef extern from "Python.h":
    PyObject* PyExc_ValueError
    void PyErr_SetString(PyObject *type, char *msg) nogil

"""
 Derived from cpython's _bisectmodule.c

PySequence_GetItem returns a new reference
PyList_GetItem returns a borrowed reference
PyList_GET_ITEM has no error checking
"""
cdef inline Py_ssize_t bisect_right(object inlist,
                                    object item,
                                    Py_ssize_t lo, Py_ssize_t hi) except -1:
    cdef object litem_cy
    cdef Py_ssize_t mid, res

    if lo < 0:
        PyErr_SetString(PyExc_ValueError, "lo must be non-negative")
        return -1
    if hi == -1:
        hi = PyList_Size(inlist)
        if hi < 0:
            return -1
    while lo < hi:
        """
        The (size_t)cast ensures that the addition and subsequent division
           are performed as unsigned operations, avoiding difficulties from
           signed overflow.  (See issue 13496.)
        """
        mid = (<size_t>lo + hi) // 2
        litem_cy = <object> PyList_GetItem(inlist, mid)
        res = PyObject_RichCompareBool(item, litem_cy, Py_LT)
        if res < 0:
            return -1
        if res:
            hi = mid
        else:
            lo = mid + 1

    return lo
# end def

cdef inline Py_ssize_t insort_right(object inlist,
                                    object item,
                                    Py_ssize_t lo, Py_ssize_t hi) except -1:
    cdef Py_ssize_t index

    index = bisect_right(inlist, item, lo, hi)
    if index < 0:
        return -1
    if PyList_CheckExact(inlist):
        if PyList_Insert(inlist, index, item) < 0:
            return -1
    else:
        PyErr_SetString(PyExc_ValueError, "inlist must be a list")
        return -1
    return 0
# end def

cdef inline Py_ssize_t bisect_left(object inlist,
                                   object item,
                                    Py_ssize_t lo, Py_ssize_t hi) except -1:
    cdef object litem_cy
    cdef Py_ssize_t mid, res

    if lo < 0:
        PyErr_SetString(PyExc_ValueError, "lo must be non-negative")
        return -1

    if hi == -1:
        hi = PyList_Size(inlist)
        if hi < 0:
            return -1;

    while lo < hi:
        """The (size_t)cast ensures that the addition and subsequent division
            are performed as unsigned operations, avoiding difficulties from
            signed overflow.  (See issue 13496.)"""
        mid = (<size_t>lo + hi) // 2
        litem_cy = <object> PyList_GetItem(inlist, mid)
        res = PyObject_RichCompareBool(litem_cy, item, Py_LT)
        if res < 0:
            PyErr_SetString(PyExc_ValueError, "lo must be non-negative")
            return -1
        if res:
            lo = mid + 1
        else:
            hi = mid

    return lo
# end def

cdef inline Py_ssize_t insort_left(object inlist,
                                    object item,
                                    Py_ssize_t lo, Py_ssize_t hi) except -1:
    cdef Py_ssize_t index

    index = bisect_left(inlist, item, lo, hi)
    if index < 0:
        return -1
    if PyList_CheckExact(inlist):
        if PyList_Insert(inlist, index, item) < 0:
            return -1
    else:
        PyErr_SetString(PyExc_ValueError, "inlist must be a list")
        return -1

    return 0
# end def


