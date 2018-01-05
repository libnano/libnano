from cpython.ref cimport PyObject

cdef Py_ssize_t bisect_right(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1

cdef Py_ssize_t insort_right(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1

cdef Py_ssize_t bisect_left(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1

cdef Py_ssize_t insort_left(PyObject* inlist, PyObject* item, Py_ssize_t lo, Py_ssize_t hi) except -1
