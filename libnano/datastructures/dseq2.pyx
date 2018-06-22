cimport c_util
cimport libnano.seqstr

from libnano.cymem.cymem cimport Pool

cdef class DSeq2:
    cdef int* fwd_gaps
    cdef int* rev_gaps
    cdef char* fwd_cstr
    cdef Py_ssize_t fwd_cstr_len
    cdef char* rev_cstr
    cdef Py_ssize_t rev_cstr_len

    cdef Pool mem   # Memory manager

    def __cinit__(self, fwd: str, rev: str):
        self.mem = None
        self.fwd_gaps = NULL
        self.rev_gaps = NULL
        self.fwd_cstr = NULL
        self.fwd_cstr_len = 0
        self.rev_cstr = NULL
        self.rev_cstr_len = 0


    def __init__(self, fwd: str, rev: str):
        self.mem = mem = Pool()
        self.fwd_cstr = c_util.obj_to_cstr_len(fwd, &self.fwd_cstr_len)
        self.rev_cstr = c_util.obj_to_cstr_len(rev, &self.rev_cstr_len)

    def awesome(self):
        print("AWESOME")
# end class