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
libnano.helpers.c_util.pxd
~~~~~~~~~~~~~~~~~~~~~~~~~~

Raw free is not available in Cython; PyMem_RawMalloc, PyMem_RawFree

'''


from cpython.mem cimport (
    PyMem_Free,
    PyMem_Malloc,
)
from libc.string cimport memcpy


cdef inline bytes _bytes(s):
    if isinstance(s, str):
        # Encode to the specific encoding used inside of the module
        return (<str> s).encode('utf8')
    else:
        return s


# Not all of these methods are available in Cython pxd file so we are importing
# them here
cdef extern from 'Python.h':
    const char* PyUnicode_AsUTF8AndSize(object, Py_ssize_t *)
    object PyUnicode_FromString(const char *)
    object PyUnicode_FromStringAndSize(const char *, Py_ssize_t)
    const char* PyUnicode_AsUTF8(object)


from cpython.bytes cimport (
    PyBytes_AsString,
    PyBytes_AsStringAndSize,
    PyBytes_Check,
    PyBytes_FromString,
    PyBytes_FromStringAndSize,
)


cdef inline object copy_obj_to_cstr_unsafe(
        object      o1,
        Py_ssize_t* length,
        char**      c_str2,
):
    '''For Python 3 it can take a bytes object or a unicode object

    This is unsafe for single characters in cpython due to object reuse
    https://github.com/python/cpython/blob/master/Objects/unicodeobject.c#L4688

    Args:
        o1 - python string object
        length - a pointer that will store the length of the string
        c_str2 - a pointer to the C NULL-terminated character array of the
            copy of o1

    Returns:
        o2 - a python string object

    Raises:
        ValueError: PyBytes_AsStringAndSize issue
        OSError: PyUnicode_AsUTF8AndSize memory error

    '''
    cdef:
        object o2
        char* c_str1
        size_t b_length

    if PyBytes_Check(o1):
        if PyBytes_AsStringAndSize(o1, &(c_str1), length) == -1:
            raise ValueError('PyBytes_AsStringAndSize issue')
        b_length = length[0]
        o2 = PyBytes_FromStringAndSize(c_str1, b_length)
        c_str2[0] = PyBytes_AsString(o2)
    else:
        c_str1 = <char*> PyUnicode_AsUTF8AndSize(o1, length)
        if c_str1 == NULL:
            raise OSError('PyUnicode_AsUTF8AndSize memory error')
        b_length = length[0]
        o2 = PyUnicode_FromStringAndSize(<const char *>c_str1, b_length)
        c_str2[0] = <char*> PyUnicode_AsUTF8(o2)

    return o2


cdef inline int copy_obj_to_cstr(
        object      o1,
        Py_ssize_t* length,
        char**      c_str2,
) except -1:
    '''Copy object string to a C NULL-terminated character array

    Args:
        o1: Python string-like object
        length: A pointer that will store the length of the string
        c_str2: A pointer to the C NULL-terminated character array of the
            copy of o1

    Returns:
        obj_type - type a python string object for Python 3
            1 == bytes
            0 == str (unicode)

    '''
    cdef:
        char* c_str1 = NULL
        char* temp = NULL
        size_t b_length
        int obj_type = 0

    if PyBytes_Check(o1):
        if PyBytes_AsStringAndSize(
            o1,
            &(c_str1),
            length,
        ) == -1:
            return -1
        obj_type = 1
    else:
        obj_type = 0
        c_str1 = <char*> PyUnicode_AsUTF8AndSize(
            o1,
            length,
        )
        if c_str1 == NULL:
            return -1
    b_length = length[0] + 1 # add 1 byte for the 0 or NULL
    temp = <char *> PyMem_Malloc(
        b_length * sizeof(char)
    )
    if temp == NULL:
        return -1
    memcpy(
        temp,
        c_str1,
        b_length
    )
    c_str2[0] = temp
    return obj_type


cdef inline cstr_to_obj(
        char*       c_str,
        Py_ssize_t  length,
        int         obj_type,
):
    '''Convert C NULL-terminated character array to a Python string-like

    Frees the memory associated with c_str

    Args:
        c_str: A pointer to the C NULL-terminated character array of the copy
            of o1
        length: The length of the string
        obj_type: for Python 3
                1 == bytes
                0 == str (unicode)

    Returns:
        A python string object
    '''
    cdef object obj
    if obj_type:
        obj = PyBytes_FromStringAndSize(
            c_str,
            length,
        )
    else:
        obj = PyUnicode_FromStringAndSize(
            <const char *> c_str,
            length,
        )
    PyMem_Free(c_str)
    return obj


cdef inline cstr_to_obj_nofree(
        char*       c_str,
        Py_ssize_t  length,
        int         obj_type,
):
    '''Convert C NULL-terminated character array to a Python string-like

    Does not free the memory associated with C character array

    Args:
        c_str: C NULL-terminated character array
        length: Length of character array to convert
        obj_type: Type to coerce to (1 bytes, 0 unicode string)

    Returns:
        Python object from the C string

    Raises:
        OSError: Memory error

    '''
    cdef object obj
    if obj_type:
        obj = PyBytes_FromStringAndSize(
            c_str,
            length,
        )
    else:
        obj = PyUnicode_FromStringAndSize(
            <const char*> c_str,
            length,
        )
    return obj


cdef inline cstr_to_obj_nolength(
        char*   c_str,
        int     obj_type,
):
    '''Convert C NULL-terminated character array to a Python string-like

    Args:
        c_str: C NULL-terminated character array
        obj_type: Type to coerce to (1 bytes, 0 unicode string)

    Returns:
        Python object from the C character array

    Raises:
        OSError: Memory error

    '''
    cdef object obj
    if obj_type:
        obj = PyBytes_FromString(c_str)
    else:
        obj = PyUnicode_FromString(
            <const char*> c_str,
        )
    return obj


cdef inline char* obj_to_cstr(
        object o1,
):
    '''Convert python bytes or unicode string to a pointer to a C
        NULL-terminated character array

    Args:
        o1: Python bytes object or a unicode string object

    Returns:
        c_str1: string pointer to the internal string of the o1

    Raises:
        OSError: Memory error

    '''
    cdef:
        char* c_str1
        Py_ssize_t length

    if PyBytes_Check(o1):
        if PyBytes_AsStringAndSize(
            o1,
            &(c_str1),
            &length,
        ) == -1:
            raise OSError('PyBytes_AsStringAndSize: memory error')
        return c_str1
    else:
        c_str1 = <char*> PyUnicode_AsUTF8AndSize(
            o1,
            &length,
        )
        if c_str1 == NULL:
            raise OSError('PyUnicode_AsUTF8AndSize: memory error')
    return c_str1


cdef inline char* obj_to_cstr_len(
        object      o1,
        Py_ssize_t* length,
):
    '''Same as above but fetches the string length too

    Args:
        o1: Python bytes object or a unicode string object
        length: pointer to length term to set

    Returns:
        pointer to C string of o1

    Raises:
        OSError: Memory error

    '''
    cdef:
        char* c_str1

    if PyBytes_Check(o1):
        if PyBytes_AsStringAndSize(
            o1,
            &(c_str1),
            length,
        ) == -1:
            raise OSError('PyBytes_AsStringAndSize: bytes error')
        return c_str1
    else:
        c_str1 = <char*> PyUnicode_AsUTF8AndSize(
            o1,
            length,
        )
        if c_str1 == NULL:
            raise OSError('PyUnicode_AsUTF8AndSize: unicode error')
    return c_str1


cdef inline char* copy_string(
        char*   src_str,
        int     length,
) except NULL:
    '''Copy C strings.

    Args:
        src_str: source string
        length: length of source string

    Returns:
        Pointer to newly allocated C string copy.  Caller is responsible for
        freeing the memory
    '''
    cdef:
        char* dest_str = <char *> PyMem_Malloc(
            (length + 1) * sizeof(char)
        )
    if dest_str == NULL:
        raise OSError('Could not allocate memory for sequence.')
    memcpy(dest_str, src_str, length + 1)
    return dest_str


cdef inline void copy_string_buffer(
        char*   src_str,
        char*   dest_str,
        int     length,
):
    '''Copy a C NULL-terminated character array.  Assumes memory is allocated
    correctly

    Args:
        src_str: source string
        dest_str: destination string
        length: length of source string
    '''
    memcpy(dest_str, src_str, length + 1)
