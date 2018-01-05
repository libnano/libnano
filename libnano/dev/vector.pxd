cdef extern from "Python.h":
    ctypedef int Py_intptr_t
    ctypedef struct PyTypeObject
    ctypedef struct PyObject:
        Py_ssize_t ob_refcnt
        PyTypeObject *ob_type

    void Py_XDECREF(PyObject *)
    void Py_XINCREF(PyObject *)

from helpers.inttypes cimport uint64_t, int32_t, int64_t

cdef struct pyobject_vec_t:
    size_t n, capacity
    PyObject **data

cdef struct int32_vec_t:
    size_t n, capacity
    int32_t *data

cdef struct int64_vec_t:
    size_t n, capacity
    int64_t *data

cdef struct size_vec_t:
    size_t n, capacity
    size_t *data

ctypedef fused fused_t:
    PyObject *
    int32_t
    int64_t
    size_t


ctypedef fused fused_vec_t:
    pyobject_vec_t
    int32_vec_t
    int64_vec_t
    size_vec_t

cdef inline fused_t init(fused_vec_t v, size_t i)
cdef inline fused_t take(fused_vec_t v, size_t i)
cdef inline resize(fused_vec_t v, size_t s, size_t tsize)
cdef inline fused_t pop(fused_vec_t v)
cdef inline size_t size(fused_vec_t v)
cdef inline push(fused_vec_t v, fused_t x, size_t tsize)
cdef inline destroy(fused_vec_t v)
cdef inline insert(fused_vec_t v, fused_t x, size_t index, size_t tsize)

cdef class ObjectVector
