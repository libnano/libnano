/******************************************************************************
** seqstr.c
**
** Python C API functions for general manipulation of DNA, RNA, and amino acid
** (i.e., computing complement / reverse / reverse complement etc).
**
******************************************************************************/

// Necessary for PyArg_ParseTuple to inject Py_ssize_t-type values into
// string size pointers
#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "ss_seqstr.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include    <numpy/arrayobject.h>


// GET_STR_AND_COPY_OBJ_3 and GET_STR_AND_COPY_STR_2
// macros create new 2/3 string or bytes objects and manipulate the strings
// in place proper encoding as utf-8 makes this possible in Python 3

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
    #define NEW_STRING_BUFFER(py_obj, char_ptr)

// support for bytes or unicode input for Python 3
// if not bytes assumes unicode
    #define GET_STR_AND_COPY_OBJ_3(o1, c_str1, c_str1_s, o2, c_str2)           \
    if (PyBytes_Check(o1)) {                                                   \
        if (PyBytes_AsStringAndSize(o1, &(c_str1), &(c_str1_s)) == -1) {       \
            return NULL;                                                       \
        }                                                                      \
        if ((o2 = PyBytes_FromStringAndSize(c_str1, c_str1_s)) == NULL) {      \
            return NULL;                                                       \
        }                                                                      \
        c_str2 = PyBytes_AsString(o2);                                         \
    } else {                                                                   \
        if ((c_str1 = PyUnicode_AsUTF8AndSize(o1, &(c_str1_s))) == NULL) {     \
            return NULL;                                                       \
        }                                                                      \
        if ((o2 = PyUnicode_FromStringAndSize(c_str1, c_str1_s)) == NULL) {    \
            return NULL;                                                       \
        }                                                                      \
        c_str2 = PyUnicode_AsUTF8(o2);                                         \
    }
#else
    // Python 2 PyString_FromStringAndSize returns a new ref so no Py_INCREF needed
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
    #define GET_STR_AND_COPY_STR_2(c_str1, c_str1_s, o2, c_str2)               \
    if( (o2 = PyString_FromStringAndSize(c_str1, c_str1_s)) == NULL) {         \
        return NULL;                                                           \
    }                                                                          \
    c_str2 = PyString_AsString(o2);
#endif


PyDoc_STRVAR(seqstr__doc__,
              "Basic operations on DNA, RNA, and amino acid sequences\n");

PyDoc_STRVAR(reverse__doc__,
    "Reverse a DNA/RNA/amino acid sequence\n\n"
    "seq: DNA, RNA, or amino acid sequence\n"
);

PyObject*
reverse(PyObject* self, PyObject* args) {

    char        *seq=NULL, *r_seq=NULL;
    Py_ssize_t  length;
    PyObject    *return_value=NULL;

    #if PY_MAJOR_VERSION >= 3
    PyObject * seq_obj;
    if (!PyArg_ParseTuple(args, "O", &seq_obj)) {
        return NULL;
    }
    GET_STR_AND_COPY_OBJ_3(seq_obj, seq, length, return_value, r_seq);
    #else
    if (!PyArg_ParseTuple(args, "s#", &seq, &length)) {
        return NULL;
    }
    GET_STR_AND_COPY_STR_2(seq, length, return_value, r_seq);
    #endif

    ss_revSeq(r_seq, (int)length);

    return return_value;

}


PyDoc_STRVAR(complement__doc__,
    "Compute the complement of a DNA/RNA/amino acid sequence\n\n"
    "seq: DNA, RNA, or amino acid sequence\n"
);

PyObject*
complement(PyObject* self, PyObject* args) {

    char        *seq=NULL, *c_seq=NULL;
    Py_ssize_t  length;
    PyObject    *return_value=NULL;

    #if PY_MAJOR_VERSION >= 3
    PyObject * seq_obj;
    if (!PyArg_ParseTuple(args, "O", &seq_obj)) {
        return NULL;
    }
    GET_STR_AND_COPY_OBJ_3(seq_obj, seq, length, return_value, c_seq);
    #else
    if (!PyArg_ParseTuple(args, "s#", &seq, &length)) {
        return NULL;
    }
    GET_STR_AND_COPY_STR_2(seq, length, return_value, c_seq);
    #endif

    ss_compSeq(c_seq, (int)length);

    return return_value;
}


PyDoc_STRVAR(reverseComplement__doc__,
    "Compute the reverse complement of a DNA/RNA/amino acid sequence\n\n"
    "seq: DNA, RNA, or amino acid sequence\n"
);

PyObject*
reverseComplement(PyObject* self, PyObject* args) {

    char        *seq=NULL, *rc_seq=NULL;
    Py_ssize_t  length;
    PyObject    *return_value=NULL;

    #if PY_MAJOR_VERSION >= 3
    PyObject * seq_obj;
    if (!PyArg_ParseTuple(args, "O", &seq_obj)) {
        return NULL;
    }
    GET_STR_AND_COPY_OBJ_3(seq_obj, seq, length, return_value, rc_seq);
    #else
    if (!PyArg_ParseTuple(args, "s#", &seq, &length)) {
        return NULL;
    }
    GET_STR_AND_COPY_STR_2(seq, length, return_value, rc_seq);
    #endif

    ss_revCompSeq(rc_seq, (int)length);

    return return_value;
}

// PyDoc_STRVAR(reverseComplementIP__doc__,
//     "Compute the reverse complement of a DNA/RNA/amino acid sequence\n\n"
//     "seq: DNA, RNA, or amino acid sequence\n"
// );

// PyObject*
// reverseComplementIP(PyObject* self, PyObject* args) {

//     char*       seq;
//     int         length;

//     if (!PyArg_ParseTuple(args, "s#", &seq, &length)) {
//         return NULL;
//     }

//     ss_revCompSeq(seq, length);

//     return Py_BuildValue("s#", seq, length);
// }

PyDoc_STRVAR(hammingDistance__doc__,
    "Compute the Hamming distance between two strings\n\n"
    "seq: DNA, RNA, or amino acid sequence\n"
);

PyObject*
hammingDistance(PyObject* self, PyObject* args) {

    char        *seq1, *seq2;
    Py_ssize_t  len1, len2;

    if (!PyArg_ParseTuple(args, "s#s#", &seq1, &len1, &seq2, &len2)) {
        return NULL;
    }

    int res = ss_hamming((const char *) seq1, (const char *) seq2, (int)len1);

    return Py_BuildValue("i", res);
}

PyDoc_STRVAR(rollingHammingDistance__doc__,
    "Compute the Hamming distance between seq1 and seq2 at each idx of seq2\n\n"
    "seq: DNA, RNA, or amino acid sequence\n"
);

PyObject*
rollingHammingDistance(PyObject* self, PyObject* args) {

    char       *seq1, *seq2;
    Py_ssize_t len1, len2;
    int        *hamming_distance_arr;


    if (!PyArg_ParseTuple(args, "s#s#", &seq1, &len1, &seq2, &len2)) {
        return NULL;
    }

    if (len1 > len2) {
        PyErr_Format(PyExc_ValueError, "rollingHammingDistance: length of seq 2 must be >= length of seq 1");
        return NULL;
    }

    // allocated as a 1D array zeroed C-Contiguous
    int num_positions = (int)len2 - (int)len1 + 1;
    npy_intp arr_dims[1] = {num_positions};
    PyArrayObject *hamming_distance_np_arr=NULL;
    hamming_distance_np_arr = (PyArrayObject *) PyArray_ZEROS(1, arr_dims, NPY_INT, 0);
    if (hamming_distance_np_arr == NULL) {
        PyErr_Format(PyExc_ValueError, "rollingHammingDistance: failed to create Numpy array");
        return NULL;
    }
    hamming_distance_arr = (int *) PyArray_DATA(hamming_distance_np_arr);
    ss_rollingHamming((const char *) seq1, (const char *) seq2, (int)len1,
                      (int)len2, hamming_distance_arr);
    return Py_BuildValue("N", hamming_distance_np_arr);
}

PyDoc_STRVAR(can3pMisprime__doc__,
    "Check for 3p mispriming (homology) between seq1 and seq2 at each idx of seq2\n\n"
    "seq: DNA, RNA, or amino acid sequence\n"
);

PyObject*
can3pMisprime(PyObject* self, PyObject* args) {

    char            *seq1, *seq2;
    int             len_thresholds, must_3p_mismatch, out, *thresholds_arr;
    Py_ssize_t      len1, len2;

    PyObject        *thresholds_obj=NULL;
    PyArrayObject   *thresholds_arr_obj=NULL;



    if (!PyArg_ParseTuple(args, "s#s#Oi", &seq1, &len1, &seq2, &len2, &thresholds_obj, &must_3p_mismatch)) {
        return NULL;
    }

    if (len1 > len2) {
        PyErr_Format(PyExc_ValueError, "can3pMisprime: length of seq 2 must be >= length of seq 1");
        return NULL;
    }

    thresholds_arr_obj = (PyArrayObject *) PyArray_FROM_OTF(thresholds_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (thresholds_arr_obj == NULL) {
        return NULL;
    }

    /* Make sure lengths match */
    len_thresholds = (int)PyArray_DIM(thresholds_arr_obj, 0);
    if (len_thresholds != len1) {
        Py_DECREF(thresholds_arr_obj);
        return NULL;
    }

    thresholds_arr = (int *) PyArray_DATA(thresholds_arr_obj);

    out = ss_can3pMisprime((const char *) seq1, (const char *) seq2, (int)len1,
                           (int)len2, thresholds_arr, must_3p_mismatch);

    Py_DECREF(thresholds_arr_obj);
    return Py_BuildValue("i", out);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static PyMethodDef seqstr_methods[] = {
    { "reverse", reverse, METH_VARARGS, reverse__doc__ },
    { "complement", complement, METH_VARARGS, complement__doc__ },
    { "reverseComplement", reverseComplement, METH_VARARGS,
                           reverseComplement__doc__ },
    { "hammingDistance", hammingDistance, METH_VARARGS,
                            hammingDistance__doc__ },
    { "rollingHammingDistance", rollingHammingDistance, METH_VARARGS,
                            rollingHammingDistance__doc__ },
    { "can3pMisprime", can3pMisprime, METH_VARARGS,
                            can3pMisprime__doc__ },
    { NULL, NULL}
};

MOD_INIT(seqstr) {
    #if PY_MAJOR_VERSION >= 3
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "seqstr",           /* m_name */
            seqstr__doc__,      /* m_doc */
            -1,                 /* m_size */
            seqstr_methods,     /* m_methods */
            NULL,               /* m_reload */
            NULL,               /* m_traverse */
            NULL,               /* m_clear */
            NULL,               /* m_free */
        };
        PyObject* m = PyModule_Create(&moduledef);
        import_array();
        return m;
    #else
        Py_InitModule3("seqstr", seqstr_methods, seqstr__doc__);
        import_array();
    #endif
};

