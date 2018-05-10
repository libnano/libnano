
cimport c_util

__doc__ = "Simple filtering of DNA sequences"

cdef extern from "sf_seqscreen.h":
    int sf_containsRun(char* seq, int seq_length,
                    int maxA, int maxT, int maxG,
                    int maxC, int maxAT, int maxGC)

    int sf_gcWindow(char* seq, int length,
                int gc_min_percent, int gc_max_percent,
                int window_size, int* window_start, int* gc_count)

def containsRun(object seq_obj, int maxA, int maxT,
                            int maxG, int maxC,
                            int maxAT, int maxGC):
    """
    Check a DNA sequence for A/T/G/C/AT/GC runs of a max length
    seq: DNA sequence
    maxA: Maximum allowable run of A
    maxT: Maximum allowable run of T
    maxG: Maximum allowable run of G
    maxC: Maximum allowable run of C
    maxAT: Maximum allowable run of AT
    maxGC: Maximum allowable run of GC
    """
    cdef int     rc
    cdef char*   seq
    cdef Py_ssize_t  length

    seq = c_util.obj_to_cstr_len(seq_obj, &length)

    rc = sf_containsRun(seq, length, maxA, maxT, maxG, maxC, maxAT, maxGC);

    return True if rc else False
# end def

def gcWindow(object seq_obj, int gc_min_percent, int gc_max_percent,
                            int window_size):
    """
    Check a DNA sequence for GC content in a sliding window of fixed size
    Returns a tuple of (Passed Filter [True/False], start of window,
                        gc count of window
    seq: DNA sequence
    min_gc_percent: minimum GC percent within window (integer from 0-99)
    max_gc_percent: maximum GC percent within window (integer from 1-100)
    window_size: size of the sliding window in bases
    """

    cdef int         rc
    cdef int         window_start = 0
    cdef int         gc_count = 0
    cdef char*       seq
    cdef Py_ssize_t  length

    seq = c_util.obj_to_cstr_len(seq_obj, &length)

    if gc_min_percent < 0 or gc_min_percent > 100:
        raise ValueError("gc_min_percent must be within 0-100")

    if gc_max_percent < 0 or gc_max_percent > 100:
        ValueError("gc_max_percent must be within 0-100")

    if gc_min_percent > gc_max_percent:
        ValueError("gc_min_percent cannot be greater than gc_max_percent")

    rc = sf_gcWindow(seq, length, gc_min_percent, gc_max_percent, window_size,
                   &window_start, &gc_count);

    if rc == -1:
        raise ValueError("window_size is larger than sequence length")

    ret = True if rc else False
    return ret, window_start, gc_count
# end def
