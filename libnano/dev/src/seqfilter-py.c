#include <Python.h>
#include "sf_seqscreen.h"

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif


PyDoc_STRVAR(seqscreen__doc__,
              "Simple filtering of DNA sequences\n");

PyDoc_STRVAR(containsRun__doc__,
    "Check a DNA sequence for A/T/G/C/AT/GC runs of a max length\n\n"
    "seq: DNA sequence\n"
    "maxA: Maximum allowable run of A\n"
    "maxT: Maximum allowable run of T\n"
    "maxG: Maximum allowable run of G\n"
    "maxC: Maximum allowable run of C\n"
    "maxAT: Maximum allowable run of AT\n"
    "maxGC: Maximum allowable run of GC\n"
);

PyObject*
containsRun(PyObject* self, PyObject* args) {
    int     seq_length, maxA, maxT, maxG, maxC, maxAT, maxGC, rc;
    char*   seq;

    if (!PyArg_ParseTuple(args, "s#iiiiii", &seq, &seq_length, &maxA, &maxT,
                          &maxG, &maxC, &maxAT, &maxGC)) {
        return NULL;
    }

    rc = sf_containsRun(seq, seq_length, maxA, maxT, maxG, maxC, maxAT, maxGC);

    return Py_BuildValue("O", rc ? Py_True : Py_False);

}


PyDoc_STRVAR(gcWindow__doc__,
    "Check a DNA sequence for GC content in a sliding window of fixed size\n\n"
    "Returns a tuple of (Passed Filter [True/False], start of window,\n"
    "                    gc count of window\n"
    "seq: DNA sequence\n"
    "min_gc_percent: minimum GC percent within window (integer from 0-99)\n"
    "max_gc_percent: maximum GC percent within window (integer from 1-100)\n"
    "window_size: size of the sliding window in bases\n"
);

PyObject*
gcWindow(PyObject* self, PyObject* args) {
    int         length, gc_min_percent, gc_max_percent, window_size, rc;
    int         window_start = 0, gc_count = 0;
    char        *seq;

    if (!PyArg_ParseTuple(args, "s#iii", &seq, &length, &gc_min_percent,
                          &gc_max_percent, &window_size)) {
        return NULL;
    }

    if (gc_min_percent < 0 || gc_min_percent > 100) {
        PyErr_SetString(PyExc_ValueError, "gc_min_percent must be within" \
                        " 0-100");
        return NULL;
    }

    if (gc_max_percent < 0 || gc_max_percent > 100) {
        PyErr_SetString(PyExc_ValueError, "gc_max_percent must be within" \
                        " 0-100");
        return NULL;
    }

    if (gc_min_percent > gc_max_percent) {
        PyErr_SetString(PyExc_ValueError, "gc_min_percent cannot be greater" \
                        " than gc_max_percent");
        return NULL;
    }

    rc = sf_gcWindow(seq, length, gc_min_percent, gc_max_percent, window_size,
                   &window_start, &gc_count);

    if (rc == -1) {
        PyErr_SetString(PyExc_ValueError, "window_size is larger than " \
                                         "sequence length");
        return NULL;
    }

    return Py_BuildValue("Oii", rc ? Py_True : Py_False, window_start,
                         gc_count);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static PyMethodDef seqscreen_methods[] = {
    { "containsRun", containsRun, METH_VARARGS, containsRun__doc__ },
    { "gcWindow", gcWindow, METH_VARARGS, gcWindow__doc__ },
    { NULL, NULL}
};

MOD_INIT(seqscreen) {
    #if PY_MAJOR_VERSION >= 3
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "seqscreen",           /* m_name */
            seqscreen__doc__,      /* m_doc */
            -1,                 /* m_size */
            seqscreen_methods,     /* m_methods */
            NULL,               /* m_reload */
            NULL,               /* m_traverse */
            NULL,               /* m_clear */
            NULL,               /* m_free */
        };
        PyObject* m = PyModule_Create(&moduledef);
        return m;
    #else
        Py_InitModule3("seqscreen", seqscreen_methods, seqscreen__doc__);
    #endif
};