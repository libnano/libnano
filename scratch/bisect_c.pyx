cdef extern from "bisectmodule.h":
    cdef Py_ssize_t internal_bisect_right(PyObject* inlist, 
                        PyObject* item, Py_ssize_t lo, Py_ssize_t hi)
    cdef Py_ssize_t internal_insort_right(PyObject* inlist, 
                        PyObject* item, Py_ssize_t lo, Py_ssize_t hi)
    cdef Py_ssize_t internal_bisect_left(PyObject* inlist, 
                    PyObject* item, Py_ssize_t lo, Py_ssize_t hi)
    cdef Py_ssize_t internal_insort_left(PyObject* inlist, 
                        PyObject* item, Py_ssize_t lo, Py_ssize_t hi)

"""
This wraps the C code.  Could easily be done in cython code alone and not
C if we dropped support for objects not PyLists
"""
cdef Py_ssize_t bisect_right(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1:
    return internal_bisect_right(inlist, item, lo, hi)

cdef Py_ssize_t insort_right(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1:
    return internal_insort_right(inlist, item, lo, hi)

cdef Py_ssize_t bisect_left(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1:
    return internal_bisect_left(inlist, item, lo, hi)

cdef Py_ssize_t insort_left(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1:
    return internal_insort_left(inlist, item, lo, hi)