#include "Python.h"

Py_ssize_t
internal_bisect_right(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi);

Py_ssize_t
internal_insort_right(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi);


Py_ssize_t
internal_bisect_left(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi);

Py_ssize_t
internal_insort_left(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi);

