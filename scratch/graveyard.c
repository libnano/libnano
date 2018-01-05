/*
DO NOT DELETE THIS FILE

These are C and Cython functions that we played with in the past but do not
intend on maintaining / using in the near term. Each function / group of
related functions should be preceeded with a brief comment explaining the
purpose and why it was relegated to this file.

*/

/*
can3pMisprime

This family of functions was intended to be a loose heuristic for identifying
potentially problematic 3 prime mispriming. The current implementation is
hard to grasp conceptually, and its reliability as a heuristic is dubious.
*/

// C function
int ss_can3pMisprime(const char *s1, const char *s2, int len1, int len2,
                     const int *thresholds, int must_3p_mismatch) {
    /*

    s1 is the sequence to compare, len1 is it's length
    s2 is the target sequence/genome to compare against, len2 is it's length

    thresholds is an array of length len1 where the values correspond to the
    mininum number of mismatches that must exist within the BASES_FROM_3P length
    from the 3 prime end of the sequence.  this array is indexed by the hamming
    distance.  for instance

        thresholds = {-1, -1, 2, 2, 2, 1, 1};

    Would mean that for a hamming distance of 0 no amount of mismatches are
    allowed and a hamming distance of 3 2 mismatches are allowed in the
    BASES_FROM_3P last bases.

    thresholds also set the minimum hamming distance requirement:

        base_threshold = sum(thresholds != -1)

    must_3p_mismatch is a flag indicating whether to enforce if the 3p end
    must mismatch in addition to the hamming distance being greater than
    the base_threshold

    return 1 means a pass
    return 0 means a fail
    */

    const int num_positions = len2 - len1 + 1;
    int hamming_distance;
    int hamming_distance_3p;
    int is_3p_mismatch;
    int num_3p_mismatches;

    char *s1_ptr;
    char *s2_ptr;
    char *s2_start_ptr = (char *) s2;
    #define BASES_FROM_3P    5
    int i;
    int base_threshold = len1;  // default to minimum

    for (i = 0; i < len1; i++) {
        if (thresholds[i] != -1) {
            base_threshold = i + 1;
            break;
        }
    }
    if (base_threshold == len1) {
        // fail if no threshold
        return 0;
    }

    for (i = 0; i < num_positions; i++) {
        s1_ptr = (char *) s1;
        s2_ptr = s2_start_ptr++;
        hamming_distance = 0;
        hamming_distance_3p = 0;
        is_3p_mismatch = 0; // reset this if set
        const char *s2_ptr_lim = s2_ptr + len1 - 1;
        const char *s2_ptr_3plim = s2_ptr + len1 - BASES_FROM_3P;
        // 1. get hamming distance from start of s1[0:(len1 - BASES_FROM_3P)]
        while (s2_ptr < s2_ptr_3plim) {
            if (*s1_ptr++ != *s2_ptr++) {
                hamming_distance++;
            }
        }
        // 2. now cover s1[(len1 - BASES_FROM_3P):len1-1] and increment separate
        // hamming_distance_3p counter on mismatch
        while (s2_ptr < s2_ptr_lim) {
            if (*s1_ptr++ != *s2_ptr++) {
                hamming_distance++;
                hamming_distance_3p++;
            }
        }
        // 3. now cover s1[len1 - 1] and increment separate
        // hamming_distance_3p counter on mismatch
        // and flag a is_3p_mismatch
        if (*s1_ptr != *s2_ptr) {
            hamming_distance++;
            hamming_distance_3p++;
            is_3p_mismatch = 1;
        }
        // printf("%s, %s\n", s1_ptr, s2_ptr);

        // 4. if hamming distance is too small
        if (hamming_distance < base_threshold) {
            // printf("hamming %d/%d\n", hamming_distance, base_threshold);
            return 0;
        } else if (hamming_distance == len1) {
            continue;
        }

        num_3p_mismatches = thresholds[hamming_distance];

        if (hamming_distance_3p >= num_3p_mismatches) {
            if (must_3p_mismatch) {
                if ((hamming_distance_3p == num_3p_mismatches) && (is_3p_mismatch == 0)) {
                    // printf("hamming %d/%d %d\n", hamming_distance_3p,
                    //     num_3p_mismatches, is_3p_mismatch);
                    return 0;
                } else {
                    continue;
                }
            } else {
                continue;
            }
        } else {
            // printf("hamming %d/%d\n", hamming_distance, base_threshold);
            return 0;
        }
    } // end for loop
    return 1;
}

// C header entry
int ss_can3pMisprime(const char *s1, const char *s2, int len1, int len2,
                     const int *thresholds, int must_3p_mismatch);

// Cython header entry
cdef inline int can3pMisprime_c(char* seq1, char* seq2, int len1, int len2,
                    int* thresholds, bint must_3p_mismatch)

// Cython functions
def can3pMisprime(object seq1_obj, object seq2_obj,
                  object thresholds_obj, bint must_3p_mismatch):
    """Check for mispriming
    Check for 3p mispriming (homology) between seq1 and seq2 at each
        idx of seq2

    Args:
        seq1_obj (str): DNA, RNA, or amino acid sequence
        seq2_obj (str):
        thresholds: an sequence of allowable mismatches indexed
            by Hamming distance of the length of sequence 1
        ``must_3p_mismatch``:

    Returns:
        ``True`` if seq1 can misprime against seq2
        ``False`` otherwise

    Raises:
        ``ValueError``
        ``OSError``
    """

    cdef char*  seq1_cstr
    cdef char*  seq2_cstr
    cdef int            len_thresholds, i, out
    cdef int*           thresholds_arr=NULL
    cdef Py_ssize_t     len1, len2

    seq1_cstr = c_util.obj_to_cstr_len(seq1_obj, &len1)

    seq2_cstr = c_util.obj_to_cstr_len(seq2_obj, &len2)

    if len1 > len2:
        raise ValueError("can3pMisprime: length of seq 2 must be >= length of seq 1")

    len_thresholds = len(thresholds_obj)
    if len_thresholds != len1:
        raise ValueError("can3pMisprime: length mismatch %d, %d" % (len1, len_thresholds))

    thresholds_arr = <int *> PyMem_Malloc(len_thresholds * sizeof(int))
    if thresholds_arr == NULL:
        raise OSError("can3pMisprime: out of memory")

    for i in range(len_thresholds):
        thresholds_arr[i] = thresholds_obj[i]

    out = can3pMisprime_c(seq1_cstr, seq2_cstr,
                          <int> len1, <int> len2,
                          thresholds_arr, must_3p_mismatch)

    PyMem_Free(thresholds_arr)
    return out
# end def


cdef inline int can3pMisprime_c(char* seq1, char* seq2, int len1, int len2,
                    int* thresholds, bint must_3p_mismatch):
    return ss_can3pMisprime(<const char *> seq1, <const char *> seq2,
                            <int> len1, <int> len2,
                            <const int*> thresholds, must_3p_mismatch)
# end def


