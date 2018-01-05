#include <Python.h>
#include <inttypes.h>
#include "si_seqint.h"

#if PY_MAJOR_VERSION >= 3
/* see http://python3porting.com/cextensions.html */
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ C API functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

PyDoc_STRVAR(seqint__doc__,
              "Compact integer representations of DNA sequences\n");

PyDoc_STRVAR(seq2Int__doc__,
    "Convert a DNA sequence to an integer representation\n\n"
    "seq: DNA sequence comprised of the bases A/a, T/t, G/g, and C/c\n"
);

static PyObject*
seq2Int(PyObject *self, PyObject *args){

    char*           seq;
    int             length;
    uint64_t        seqint;

    if (!PyArg_ParseTuple(args, "s#", &seq, &length)) {
        return NULL;
    }

    if (length > 30) {
        PyErr_SetString(PyExc_ValueError, "Sequences over 30 bp cannot be converted");
        return NULL;
    }

    seqint = si_seq2Int(seq, length);

    return Py_BuildValue("K", seqint);
}


PyDoc_STRVAR(reverseComplement__doc__,
    "Return the reverse complement seqint representation of `seqint`\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
reverseComplement(PyObject *self, PyObject *args){

    uint64_t         seqint, seqintrc;
    int              length;

    if (!PyArg_ParseTuple(args, "Ki", &seqint, &length)) {
        return NULL;
    }

    seqintrc = si_revCompSeqInt(seqint, length);

    return Py_BuildValue("K", seqintrc);
}

PyDoc_STRVAR(complement__doc__,
    "Return the complement seqint representation of `seqint`\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
complement(PyObject *self, PyObject *args){

    uint64_t         seqint, seqintrc;
    int              length;

    if (!PyArg_ParseTuple(args, "Ki", &seqint, &length)) {
        return NULL;
    }

    seqintrc = si_compSeqInt(seqint, length);

    return Py_BuildValue("K", seqintrc);
}

PyDoc_STRVAR(reverse__doc__,
    "Return the reverse complement seqint representation of `seqint`\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
reverse(PyObject *self, PyObject *args){

    uint64_t         seqint, seqintrc;
    int              length;

    if (!PyArg_ParseTuple(args, "Ki", &seqint, &length)) {
        return NULL;
    }

    seqintrc = si_revSeqInt(seqint, length);

    return Py_BuildValue("K", seqintrc);
}

PyDoc_STRVAR(addBase__doc__,
    "Add a base to the right-hand side of a seqint\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "base: base (A/a, T/t, G/g, or C/c) to be added\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
addBase(PyObject *self, PyObject *args){

    char*       base;
    uint64_t    seqint, seqintmod;

    if (!PyArg_ParseTuple(args, "Ks", &seqint, &base)) {
        return NULL;
    }

    seqintmod = si_addBase(seqint, base[0]);

    return Py_BuildValue("K", seqintmod);
}


PyDoc_STRVAR(removeBase__doc__,
    "Remove a base from the right-hand side of a seqint\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
removeBase(PyObject *self, PyObject *args){

    uint64_t    seqint, seqintmod;
    int         length;

    if (!PyArg_ParseTuple(args, "Ksi", &seqint, &length)) {
        return NULL;
    }

    seqintmod = si_removeBase(seqint);

    return Py_BuildValue("K", seqintmod);
}


PyDoc_STRVAR(addToWindow__doc__,
    "Add a base to the right-hand side of a seqint (fixed length window)\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "base: base (A/a, T/t, G/g, or C/c) to be added\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
addToWindow(PyObject *self, PyObject *args){

    char*       base;
    uint64_t    seqint, seqintmod;
    int         length;

    if (!PyArg_ParseTuple(args, "Ksi", &seqint, &base, &length)) {
        return NULL;
    }

    seqintmod = si_addToWindow(seqint, base[0], length);

    return Py_BuildValue("K", seqintmod);
}

PyDoc_STRVAR(getSubstring__doc__,
    "Return the seqint representing the defined substring\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "start_idx: substring start index (inclusive, 0-based)\n"
    "start_idx: substring end index (exclusive)\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
getSubstring(PyObject *self, PyObject *args){

    uint64_t    seqint, seqintmod;
    int         length, start_idx, end_idx;

    if (!PyArg_ParseTuple(args, "Kiii", &seqint, &start_idx, &end_idx,
                          &length)) {
        return NULL;
    }

    seqintmod = si_seqIntSubstring(seqint, start_idx, end_idx, length);

    return Py_BuildValue("K", seqintmod);
}


PyDoc_STRVAR(int2Seq__doc__,
    "Return the DNA sequence string represented by `seqint`\n\n"
    "seqint: integer representation of a DNA sequence\n"
    "length: length of the represented sequence in bases\n"
);

static PyObject*
int2Seq(PyObject *self, PyObject *args){

    uint64_t    seqint;
    int         length;
    PyObject*   ret_obj;

    if (!PyArg_ParseTuple(args, "Ki", &seqint, &length)) {
        return NULL;
    }

    if (length == 0)
        return Py_BuildValue("s", "");

    char* out_seq = (char *) malloc((length+1)*sizeof(char));
    if (out_seq == NULL) {
        PyErr_SetString(PyExc_IOError, "Could not allocate memory for sequence.");
        return NULL;
    }

    if (si_int2Seq(seqint, out_seq, length)) {
        free(out_seq);
        PyErr_SetString(PyExc_IOError, "Could not convert integer to sequence.");
        return NULL;
    }

    ret_obj = Py_BuildValue("s", out_seq);
    free(out_seq);

    return ret_obj;
}

/* ~~~~~~~~~~~~~~~~~~~~~ Needs documentation / not used ~~~~~~~~~~~~~~~~~~~~ */

// PyDoc_STRVAR(hasNoRepeats__doc__,
//     "Check for repeats of 6 A's, 6 T's, 3 C's and 3 G's\n\n"
//     "1 means no repeats, 0 means there are repeats\n\n"
// );

// static PyObject*
// hasNoRepeats(PyObject *self, PyObject *args){
//     uint64_t  seqint;
//     uint64_t seqint_copy;
//     uint64_t  temp;
//     int   i, il, bases;
//     int res = 1;
//     if (!PyArg_ParseTuple(args, "Ki", &seqint, &bases)) {
//         return NULL;
//     }
//     seqint_copy = seqint;

//     // separate by mask length: 3 bases
//     for (i=0, il=bases-3+1; i < il; i++){
//         temp = SI_MASK3 & seqint_copy;
//         if (temp == SI_G_X_3) {
//             res = 0;
//             break;
//         }
//         if (temp == SI_C_X_3) {
//             res = 0;
//             break;
//         }
//         seqint_copy >>= 2;
//     }

//     // separate by mask length: 6 bases
//     for (i=0, il=bases-6+1; i < il; i++){
//         temp = seqint & SI_MASK6;
//         if (temp == SI_A_X_6) {
//             res = 0;
//             break;
//         }
//         if (temp == SI_T_X_6) {
//             res = 0;
//             break;
//         }
//         seqint >>= 2;
//     } 

//     return Py_BuildValue("i", res);
// }

// /*****************************************************************************/
// /* begin from _bitarray.c by Ian Schnell
// * add compatibility with bitarray objects for repeats
// * see https://github.com/ilanschnell/bitarray/blob/master/bitarray/_bitarray.c
// */
// typedef long long int idx_t;
// #define BITMASK(endian, i)  (((char) 1) << ((endian) ? (7 - (i)%8) : (i)%8))
// static void
// setbit(char* ob_item, int endian, int i, int bit) {
//     char *cp, mask;

//     mask = BITMASK(endian, i);

//     cp = ob_item + i / 8;
//     if (bit)
//         *cp |= mask;
//     else
//         *cp &= ~mask;
// }

// typedef struct {
//     PyObject_VAR_HEAD
// #ifdef WITH_BUFFER
//     int ob_exports;             /* how many buffer exports */
// #endif
//     char *ob_item;
//     Py_ssize_t allocated;        how many bytes allocated 
//     idx_t nbits;                /* length og bitarray */
//     int endian;                 /* bit endianness of bitarray */
//     PyObject *weakreflist;      /* list of weak references */
// } bitarrayobject;

// /* end from _bitarray.c */ 
// /*****************************************************************************/
// PyDoc_STRVAR(setRepeats__doc__,
//     "set repeats\n\n"
// );

// static PyObject*
// setRepeats(PyObject *self, PyObject *args){
//     PyObject *bitarr_obj;
//     int  num_bases;

//     bitarrayobject * bitarr;
//     char * obitem;
    
//     uint64_t filter1;
//     uint64_t filter2;

//     uint64_t mask, nmask;

//     int   i, il, j, jl;
//     int res = 1;
//     if (!PyArg_ParseTuple(args, "O", &bitarr_obj)) {
//         return NULL;
//     }
//     bitarr = (bitarrayobject *) bitarr_obj;
//     obitem = bitarr->ob_item;
//     int nbits = bitarr->nbits;
//     num_bases = (int) log2(nbits)/2;
//     // printf("this many %d\n", num_bases); 
//     int endian = bitarr->endian;

//     // group by length: 3 bases
//     filter1 = SI_G_X_3;
//     filter2 = SI_C_X_3;
//     mask = SI_MASK3;

//     for (i=0, il=num_bases-3+1; i < il; i++) {
//         nmask = ~mask;
//         for (j=0; j < nbits; j++) {
//             jl = j & nmask;
//             setbit(obitem, endian, jl | filter1, 1);
//             setbit(obitem, endian, jl | filter2, 1);
//         }
//         filter1 <<= 2;
//         filter2 <<= 2;
//         mask <<= 2;
//     }

//     // group by length: 6 bases
//     filter1 = SI_A_X_6;
//     filter2 = SI_T_X_6;
//     mask = SI_MASK6;

//     for (i=0, il=num_bases-6+1; i < il; i++) {
//         nmask = ~mask;
//         for (j=0; j < nbits; j++) {
//             jl = j & nmask;
//             setbit(obitem, endian, jl | filter1, 1);
//             setbit(obitem, endian, jl | filter2, 1);
//         }
//         filter1 <<= 2;
//         filter2 <<= 2;
//         mask <<= 2;
//     }

//     return Py_BuildValue("i", res);
// };

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

static PyMethodDef seqint_methods[] = {
    { "seq2Int", seq2Int, METH_VARARGS, seq2Int__doc__ },
    { "reverseComplement", reverseComplement, METH_VARARGS, reverseComplement__doc__ },
    { "complement", complement, METH_VARARGS, complement__doc__ },
    { "reverse", reverse, METH_VARARGS, reverse__doc__ },
    { "addBase", addBase, METH_VARARGS, addBase__doc__ },
    { "removeBase", removeBase, METH_VARARGS, removeBase__doc__ },
    { "addToWindow", addToWindow, METH_VARARGS, addToWindow__doc__ },
    { "getSubstring", getSubstring, METH_VARARGS, getSubstring__doc__ },
    { "int2Seq", int2Seq, METH_VARARGS, int2Seq__doc__ },
    /* ~~~~~~~~~~~~~~~~~~~ Needs documentation / not used ~~~~~~~~~~~~~~~~~~ */
    // { "hasNoRepeats", hasNoRepeats, METH_VARARGS, hasNoRepeats__doc__ },
    // { "setRepeats", setRepeats, METH_VARARGS, setRepeats__doc__ },
    { NULL, NULL}
};

MOD_INIT(seqint) {
    #if PY_MAJOR_VERSION >= 3
        static struct PyModuleDef moduledef = {
            PyModuleDef_HEAD_INIT,
            "seqint",           /* m_name */
            seqint__doc__,      /* m_doc */
            -1,                 /* m_size */
            seqint_methods,     /* m_methods */
            NULL,               /* m_reload */
            NULL,               /* m_traverse */
            NULL,               /* m_clear */
            NULL,               /* m_free */
        };
        PyObject* m = PyModule_Create(&moduledef);
        return m;
    #else
        Py_InitModule3("seqint", seqint_methods, seqint__doc__);
    #endif
};
