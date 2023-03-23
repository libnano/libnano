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
libnano.list_bisect.pxd
~~~~~~~~~~~~~~~~~~~~~~~

 Derived from cpython's _bisectmodule.c

PySequence_GetItem returns a new reference
PyList_GetItem returns a borrowed reference
PyList_GET_ITEM has no error checking
'''
from cpython.list cimport (
    PyList_CheckExact,
    PyList_GetItem,
    PyList_Insert,
    PyList_Size,
)
from cpython.object cimport (
    Py_LT,
    PyObject_RichCompareBool,
)
from cpython.ref cimport PyObject


cdef extern from 'Python.h':
    PyObject* PyExc_ValueError
    void PyErr_SetString(PyObject *type, char *msg) nogil


cdef inline Py_ssize_t bisect_right(
        object inlist,
        object item,
        Py_ssize_t lo,
        Py_ssize_t hi,
) except -1:
    cdef:
        object litem_cy
        Py_ssize_t mid, res

    if lo < 0:
        PyErr_SetString(
            PyExc_ValueError,
            'lo must be non-negative',
        )
        return -1
    if hi == -1:
        hi = PyList_Size(inlist)
        if hi < 0:
            return -1
    while lo < hi:
        # The (size_t)cast ensures that the addition and subsequent division
        #    are performed as unsigned operations, avoiding difficulties from
        #    signed overflow.  (See issue 13496.)

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


cdef inline Py_ssize_t insort_right(
        object inlist,
        object item,
        Py_ssize_t lo,
        Py_ssize_t hi,
) except -1:
    cdef Py_ssize_t index

    index = bisect_right(
        inlist,
        item,
        lo,
        hi,
    )
    if index < 0:
        return -1
    if PyList_CheckExact(inlist):
        if PyList_Insert(
            inlist,
            index,
            item,
        ) < 0:
            return -1
    else:
        PyErr_SetString(
            PyExc_ValueError,
            'inlist must be a list',
        )
        return -1
    return 0


cdef inline Py_ssize_t bisect_left(
        object inlist,
        object item,
        Py_ssize_t lo,
        Py_ssize_t hi,
) except -1:
    cdef:
        object litem_cy
        Py_ssize_t mid, res

    if lo < 0:
        PyErr_SetString(
            PyExc_ValueError,
            'lo must be non-negative'
        )
        return -1

    if hi == -1:
        hi = PyList_Size(inlist)
        if hi < 0:
            return -1;

    while lo < hi:
        # The (size_t)cast ensures that the addition and subsequent division
        #     are performed as unsigned operations, avoiding difficulties from
        #     signed overflow.  (See issue 13496.)
        mid = (<size_t>lo + hi) // 2
        litem_cy = <object> PyList_GetItem(inlist, mid)
        res = PyObject_RichCompareBool(
            litem_cy,
            item,
            Py_LT,
        )
        if res < 0:
            PyErr_SetString(
                PyExc_ValueError,
                'lo must be non-negative'
            )
            return -1
        if res:
            lo = mid + 1
        else:
            hi = mid

    return lo


cdef inline Py_ssize_t insort_left(
        object inlist,
        object item,
        Py_ssize_t lo,
        Py_ssize_t hi,
) except -1:
    cdef Py_ssize_t index

    index = bisect_left(
        inlist,
        item,
        lo,
        hi,
    )
    if index < 0:
        return -1
    if PyList_CheckExact(inlist):
        if PyList_Insert(
            inlist,
            index,
            item,
        ) < 0:
            return -1
    else:
        PyErr_SetString(
            PyExc_ValueError,
            'inlist must be a list',
        )
        return -1

    return 0
