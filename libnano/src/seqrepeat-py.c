#include <Python.h>
#include "structmember.h"
#include <string.h> /* for NULL pointers */
#include "sr_seqrepeat.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include    <numpy/arrayobject.h>

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C API functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

PyDoc_STRVAR(seqrepeat__doc__,
              "Find repeats in a DNA sequence\n");


typedef struct {
    PyObject_HEAD
    repeatcheck_t *rcheck;
    PyObject* is_extended;
} RepeatCheck;

static void
RepeatCheck_dealloc(RepeatCheck* self) {
    sr_freeRepeatCheck(self->rcheck);
    PyMem_Free(self->rcheck);
    self->rcheck = NULL;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
RepeatCheck_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    RepeatCheck *self;
    self = (RepeatCheck *)type->tp_alloc(type, 0);
    self->rcheck = NULL;
    self->is_extended = Py_False;
    return (PyObject *)self;
}

repeatcheck_t * RepeatCheck_data(RepeatCheck* self) {
    return self->rcheck;
};

static int
RepeatCheck_init(RepeatCheck *self, PyObject *args, PyObject *kwds) {
    char * seq;
    int seq_length;
    int seed_size;
    int max_seed_size;

    if (self == NULL) {
        return -1;
    }

    if (!PyArg_ParseTuple(args, "s#ii", &seq, &seq_length,
            &seed_size, &max_seed_size))  {
        PyErr_SetString(PyExc_TypeError, "RepeatCheck_init: arg problem\n");
        return -1;
    }
    repeatcheck_t *rck = PyMem_Malloc(sizeof(repeatcheck_t));
    if (rck == NULL) {
        PyErr_SetString(PyExc_TypeError, "RepeatCheck_init: problem allocating repeatcheck_t\n");
        return -1;
    } else {
        self->rcheck = rck;
        rck->seqint_arr = NULL;
        rck->repeatroot_arr = NULL;
        rck->repeat_idx_arrs = NULL;
    }
    if (sr_buildRepeatData(seq, seq_length, seed_size, max_seed_size, &rck) < 0) {
        PyErr_SetString(PyExc_TypeError, "RepeatCheck_init: problem building\n");
        PyMem_Free(rck);
        return -1;
    }

    return 0;
};

static PyObject *
RepeatCheck_buildExtendedRepeats(RepeatCheck* self, PyObject *args) {
    if (sr_buildExtendedRepeats(self->rcheck) < 0) {
        PyErr_SetString(PyExc_TypeError, "RepeatCheck_sr_buildExtendedRepeats: problem with repeats\n");
        return NULL;
    }
    self->is_extended = Py_True;
    Py_RETURN_NONE;
};

static PyObject *
RepeatCheck_window(RepeatCheck* self, PyObject *args) {
    int test_word_size, i;
    int window_size, max_repeat_percent, max_repeat_count;
    int *repeat_counts = NULL;
    int *repeat_violation_idxs = NULL;
    int lim_violation_count;

    if (!PyArg_ParseTuple(args, "iii", &test_word_size, &window_size, &max_repeat_percent))  {
        return NULL;
    }

    max_repeat_count = max_repeat_percent * window_size / 100;
    int counts_size = self->rcheck->seq_length - window_size + 1;

    // allocated as a 1D array zeroed C-Contiguous
    npy_intp counts_dims[1] = {counts_size};
    PyArrayObject *counts_arr = NULL;
    counts_arr = (PyArrayObject *) PyArray_ZEROS(1, counts_dims, NPY_INT, 0);
    if (counts_arr == NULL) {
        PyErr_Format(PyExc_ValueError, "RepeatCheck_window: failed to create Numpy array");
        return NULL;
    }
    repeat_counts = (int *) PyArray_DATA(counts_arr);

    if (sr_repeatWindow(self->rcheck, test_word_size, window_size,
            max_repeat_count, &lim_violation_count, &repeat_violation_idxs, &repeat_counts) < 0) {
        PyErr_SetString(PyExc_TypeError, "RepeatCheck_window: problem with sr_repeatWindow\n");
        return NULL;
    }

    PyObject* outtuple_violations = NULL;
    outtuple_violations = PyTuple_New(lim_violation_count);
    if (outtuple_violations == NULL) {
        goto rcw_fail;
    }
    for (i=0; i < lim_violation_count; i++) {
        PyTuple_SET_ITEM(outtuple_violations, i,
            PyLong_FromLong((long) repeat_violation_idxs[i]));
    }
    PyMem_Free(repeat_violation_idxs);
    return Py_BuildValue("NO", counts_arr, outtuple_violations);
rcw_fail:
    PyErr_SetString(PyExc_MemoryError, "RepeatCheck_window: memory error issue\n");
    PyMem_Free(repeat_violation_idxs);
    return NULL;
};

static PyObject*
RepeatCheck_pileup(RepeatCheck* self, PyObject *args) {
    int min_word_size, max_word_size;
    int *pileup_array;
    PyArrayObject *pileup_npy_array=NULL;

    if (self->is_extended == Py_False) {
        PyErr_SetString(PyExc_ValueError, "RepeatCheck_pileup: call buildExtendedRepeats first\n");
        return NULL;
    }
    if (!PyArg_ParseTuple(args, "ii", &min_word_size, &max_word_size))  {
        return NULL;
    }
    npy_intp pileup_array_len = self->rcheck->seq_length;
    pileup_npy_array = (PyArrayObject *) PyArray_ZEROS(1,
                                 &pileup_array_len, NPY_INT, 0);
    pileup_array = (int *) PyArray_DATA(pileup_npy_array);
    sr_repeatPileup(self->rcheck, min_word_size, max_word_size, pileup_array);

    return (PyObject*) pileup_npy_array; //Py_BuildValue("N", pileup_npy_array);

}


static PyMemberDef RepeatCheck_members[] = {
    {NULL}  /* Sentinel */
};

static PyGetSetDef RepeatCheck_getsetters[] = {
    {NULL}  /* Sentinel */
};

static PyMethodDef RepeatCheck_methods[] = {
    {"buildExtendedRepeats", (PyCFunction)RepeatCheck_buildExtendedRepeats, METH_VARARGS,
     "build extended repeats\n"
    },
    {"window", (PyCFunction)RepeatCheck_window, METH_VARARGS,
     "get the repeat window\n"
    },
    {"pileup", (PyCFunction)RepeatCheck_pileup, METH_VARARGS,
     "generate pileup of repeat counts per base index\n"},
    {NULL}  /* Sentinel */
};

static PyTypeObject RepeatCheckType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "seqrepeat.RepeatCheck",        /*tp_name*/
    sizeof(RepeatCheck),            /*tp_basicsize*/
    0,                              /*tp_itemsize*/
    (destructor)RepeatCheck_dealloc,/*tp_dealloc*/
    0,                              /*tp_print*/
    0,                              /*tp_getattr*/
    0,                              /*tp_setattr*/
    0,                              /*tp_compare*/
    0,                              /*tp_repr*/
    0,                              /*tp_as_number*/
    0,                              /*tp_as_sequence*/
    0,                              /*tp_as_mapping*/
    0,                              /*tp_hash */
    0,                              /*tp_call*/
    0,                              /*tp_str*/
    0,                              /*tp_getattro*/
    0,                              /*tp_setattro*/
    0,                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "RepeatCheck objects",          /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    RepeatCheck_methods,            /* tp_methods */
    RepeatCheck_members,            /* tp_members */
    RepeatCheck_getsetters,          /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)RepeatCheck_init,     /* tp_init */
    0,                              /* tp_alloc */
    RepeatCheck_new,                /* tp_new */
};

static PyMethodDef seqrepeat_mod_methods[] = {
    {NULL}
};

MOD_INIT(seqrepeat) {
    if (PyType_Ready(&RepeatCheckType) < 0) {
    #if PY_MAJOR_VERSION >= 3
        return NULL;
    #else
        return;
    #endif
    }
    #if PY_MAJOR_VERSION >= 3
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "seqrepeat",            /* m_name */
            seqrepeat__doc__,       /* m_doc */
            -1,                     /* m_size */
            seqrepeat_mod_methods,  /* m_methods */
            NULL,                   /* m_reload */
            NULL,                   /* m_traverse */
            NULL,                   /* m_clear */
            NULL,                   /* m_free */
        };
        PyObject* m = PyModule_Create(&moduledef);
        import_array();
        if (m == NULL) { return NULL; }

        Py_INCREF(&RepeatCheckType);
        PyModule_AddObject(m, "RepeatCheck", (PyObject *)&RepeatCheckType);
        return m;
    #else
        PyObject* m = Py_InitModule3("seqrepeat", seqrepeat_mod_methods, seqrepeat__doc__);
        if (m == NULL) { return; }
        import_array();
        Py_INCREF(&RepeatCheckType);
        PyModule_AddObject(m, "RepeatCheck", (PyObject *)&RepeatCheckType);
    #endif
};