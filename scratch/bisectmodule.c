 /* Derived from cpython _bisectmodule.c

 Bisection algorithms. Drop in replacement for bisect.py
Converted to C by Dmitry Vasiliev (dima at hlabs.spb.ru).
*/

#define PY_SSIZE_T_CLEAN
#include "Python.h"

_Py_IDENTIFIER(insert); // for weird non list edge case for custom types

Py_ssize_t
internal_bisect_right(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi)
{
    PyObject *litem;
    Py_ssize_t mid, res;

    if (lo < 0) {
        PyErr_SetString(PyExc_ValueError, "lo must be non-negative");
        return -1;
    }
    if (hi == -1) {
        hi = PySequence_Size(list);
        if (hi < 0)
            return -1;
    }
    while (lo < hi) {
        /* The (size_t)cast ensures that the addition and subsequent division
           are performed as unsigned operations, avoiding difficulties from
           signed overflow.  (See issue 13496.) */
        mid = ((size_t)lo + hi) / 2;
        litem = PySequence_GetItem(list, mid);
        if (litem == NULL)
            return -1;
        res = PyObject_RichCompareBool(item, litem, Py_LT);
        Py_DECREF(litem);
        if (res < 0)
            return -1;
        if (res)
            hi = mid;
        else
            lo = mid + 1;
    }
    return lo;
}

Py_ssize_t
internal_insort_right(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi)
{
    PyObject *result;
    Py_ssize_t index;

    index = internal_bisect_right(list, item, lo, hi);
    if (index < 0)
        return -1;
    if (PyList_CheckExact(list)) {
        if (PyList_Insert(list, index, item) < 0)
            return -1;
    } else {
        result = _PyObject_CallMethodId(list, &PyId_insert, "nO", index, item);
        if (result == NULL)
            return -1;
        Py_DECREF(result);
    }
    return 0;
}


Py_ssize_t
internal_bisect_left(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi)
{
    PyObject *litem;
    Py_ssize_t mid, res;

    if (lo < 0) {
        PyErr_SetString(PyExc_ValueError, "lo must be non-negative");
        return -1;
    }
    if (hi == -1) {
        hi = PySequence_Size(list);
        if (hi < 0)
            return -1;
    }
    while (lo < hi) {
        /* The (size_t)cast ensures that the addition and subsequent division
           are performed as unsigned operations, avoiding difficulties from
           signed overflow.  (See issue 13496.) */
        mid = ((size_t)lo + hi) / 2;
        litem = PySequence_GetItem(list, mid);
        if (litem == NULL)
            return -1;
        res = PyObject_RichCompareBool(litem, item, Py_LT);
        Py_DECREF(litem);
        if (res < 0)
            return -1;
        if (res)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

Py_ssize_t
internal_insort_left(PyObject *list, PyObject *item, Py_ssize_t lo, Py_ssize_t hi)
{
    PyObject *result;
    Py_ssize_t index;

    index = internal_bisect_left(list, item, lo, hi);
    if (index < 0)
        return -1;
    if (PyList_CheckExact(list)) {
        if (PyList_Insert(list, index, item) < 0)
            return -1;
    } else {
        result = _PyObject_CallMethodId(list, &PyId_insert, "nO", index, item);
        if (result == NULL)
            return -1;
        Py_DECREF(result);
    }

    return 0;
}

