
from libnano.helpers cimport c_util
from cpython.ref cimport Py_INCREF
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libnano.cymem.cymem cimport Pool
from cython.operator cimport postincrement as inc
from libc.string cimport memcpy
from libc.stdio cimport printf

cdef extern from "shl_seedhashlist.h":
    # leave out a few members
    ctypedef struct tuple_t:
        int idx
        int offset

    ctypedef struct seed_search_obj_t:
        pass

    int shl_seed_weight(char*, int)
    int shl_parse_seed(char*,int,int*)
    int shl_free_hash_obj(seed_search_obj_t*)
    int shl_buildSeedTable(int*, int, int, char*, int, int*, int*, int, seed_search_obj_t*)
    int shl_findSingleMatches( int, char*, int, char*, int, int*, int*, seed_search_obj_t*, tuple_t**, int)
    int checkht(seed_search_obj_t*)

cdef class SeedMatcher:
    """ Generate the hash table from a list of
    sequences for submers of lenght m with k mismatches

    uses seed.pyx to generate seeds or seed pairs
    """
    cdef list reference_list
    cdef Py_ssize_t reference_list_length

    cdef char* ref_cstr_buffer  # buffer to put all c-strings of
    cdef Py_ssize_t ref_cstr_buf_length

    cdef int* reference_lengths     # length of each string in library
    cdef int* reference_start_idxs    # start indices of each string in the library

    cdef int* seed_idxs
    cdef int seed_idxs_len

    cdef object seed
    cdef Py_ssize_t seed_len

    cdef char* seed_cstr
    cdef seed_search_obj_t* search_obj

    cdef Pool mem   # Memory manager

    def __cinit__(self, seed, list reference_list, debug=False):
        self.mem = None
        self.reference_list = None
        self.reference_list_length = 0
        self.seed = None
        self.seed_len = 0

        # the concatenated c string version of the reference, null terminated
        self.ref_cstr_buffer = NULL

        # this is length of each string in reference_list
        self.reference_lengths = NULL

        # this is the start index of each string in reference_list in
        #       self.ref_cstr_buffer
        self.reference_start_idxs = NULL


        self.seed_idxs = NULL
        self.seed_cstr = NULL
        self.search_obj = NULL

    def __init__(self, seed, list reference_list, debug=False):
        cdef Py_ssize_t i, j, k
        cdef Py_ssize_t reference_list_length, total_seq_length
        cdef Py_ssize_t seed_len
        cdef int offset_last, check

        cdef char* temp_ref_cstr_buff = NULL
        cdef int* temp_ref_lengths = NULL
        cdef int* temp_ref_start_idxs = NULL
        cdef char* seq_cstr = NULL
        cdef Py_ssize_t seq_cstr_len

        self.mem = mem = Pool()
        self.seed = seed

        # 1 Setup the seed
        self.seed_cstr = c_util.obj_to_cstr_len(seed, &seed_len)
        self.seed_idxs_len = shl_seed_weight(self.seed_cstr, seed_len)
        self.seed_idxs = <int*> mem.malloc(self.seed_idxs_len, sizeof(int))

        shl_parse_seed(self.seed_cstr, <int> seed_len, self.seed_idxs)

        self.seed_len = seed_len

        # 2. Get constants from the list of sequences
        self.reference_list = reference_list
        reference_list_length = len(reference_list)
        total_seq_length = sum(len(x) for x in reference_list)


        # 3. Memory allocations
        self.ref_cstr_buf_length = total_seq_length + reference_list_length
        temp_ref_cstr_buff = <char*> mem.malloc(self.ref_cstr_buf_length, sizeof(char))
        temp_ref_lengths = <int*> mem.malloc(reference_list_length, sizeof(int))
        temp_ref_start_idxs = <int*> mem.malloc(reference_list_length, sizeof(int))

        self.search_obj = <seed_search_obj_t *> mem.calloc(1, sizeof(seed_search_obj_t));


        # 4. copy strings into one buffer
        i = 0
        k = 0
        offset_last = 0
        for seq in reference_list:
            seq_cstr = c_util.obj_to_cstr_len(seq, &seq_cstr_len)
            temp_ref_lengths[k] = <int> inc(seq_cstr_len)  # notice the inc
            temp_ref_start_idxs[inc(k)] = offset_last
            offset_last += <int> seq_cstr_len
            # copy everything including the NULL
            memcpy(&temp_ref_cstr_buff[i], seq_cstr, seq_cstr_len)
            i += seq_cstr_len
        # end for
        #print(i, self.ref_cstr_buf_length)

        # assign pointers
        self.ref_cstr_buffer = temp_ref_cstr_buff
        self.reference_lengths = temp_ref_lengths
        self.reference_start_idxs = temp_ref_start_idxs

        if debug:
            print("building table", self.seed_idxs[0])
        check = shl_buildSeedTable(
            self.seed_idxs,
            self.seed_idxs_len,
            self.seed_len,
            self.ref_cstr_buffer,
            self.ref_cstr_buf_length,
            self.reference_lengths,
            self.reference_start_idxs,
            reference_list_length,
            self.search_obj)
        if debug:
            print("table built")
        if check < 0:
            raise OSError("something happened")
    #  end def

    def __dealloc__(self):
        #print("freeing")
        shl_free_hash_obj(self.search_obj)
    # end def

    cdef list toList(self, tuple_t* matches, int repeat_count):
        cdef Py_ssize_t i
        cdef list outlist = [None]*repeat_count
        for i in range(repeat_count):
            outlist[i] = (matches[i].idx, matches[i].offset)
        return outlist
    # end def

    def _checkht(self):
        checkht(self.search_obj)

    def match(self, target, int mismatches):
        cdef int repeat_count
        cdef char* target_cstr = NULL
        cdef Py_ssize_t target_len
        cdef tuple_t* matches = NULL


        cdef int matches_buffer_size = self.ref_cstr_buf_length;

        target_cstr = c_util.obj_to_cstr_len(target, &target_len)

        repeat_count = shl_findSingleMatches(
                            mismatches,
                            target_cstr,
                            target_len,
                            self.ref_cstr_buffer,
                            self.ref_cstr_buf_length,
                            self.reference_lengths,
                            self.reference_start_idxs,
                            self.search_obj,
                            &matches,
                            matches_buffer_size
                            )
        if repeat_count > 0:
            #print("Found", repeat_count)
            self.mem.own(matches, repeat_count, sizeof(tuple_t))
            outlist = self.toList(matches, repeat_count)
            self.mem.free(matches)
            return outlist
        else:
            return []
    # end def
# end class
