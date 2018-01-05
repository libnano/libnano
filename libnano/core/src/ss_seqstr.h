/******************************************************************************
** ss_seqstr.h
**
** C functions and respective Python C API
** functions for general manipulation of nucleic acid sequences (i.e.,
** computing complement / reverse / reverse complement etc).
**
******************************************************************************/

#ifndef SEQ_STR_H
#define SEQ_STR_H


// Compute the complement of `seq` in place.
int ss_revSeq(char* seq, int length);

// Compute the complement of `seq` in place.
int ss_compSeq(char* seq, int length);

// Compute the reverse complement of `seq` in place.
int ss_revCompSeq(char* seq, int length);
int ss_revCompSeqCopy(char *dest, const char *from, int length);

// Compute the hamming distance between s1 and s2 of
// equal lengths. The hamming distance is undefined
// for sequences of different lengths.
int ss_hamming(const char *s1, const char *s2, int len);

// Compute the hamming distance between s1 and s2 at each index of
// s2.
int ss_rollingHamming(const char *s1, const char *s2,
                        int len1, int len2,
                        int overlap,
                        int* hamming_arr_ptr);

/* check s1 versus a longer s2 in a sliding window and return the
minimum hamming distance found by scanning all of s2 versus s1
*/
int ss_minHamming(const char *s1, const char *s2, int len1, int len2);


/* check s1 versus a longer s2 in a sliding window and if any hamming
distance is less than a threshold return immediately.  Otherwise return the
minimum hamming distance found by scanning all of s2 versus s1
*/
int ss_minHammingThreshold(const char *s1, const char *s2, int len1, int len2, int threshold);


/* do hamming checks of s1 versus s2 both of length `len` in the following
* four (4) ways
* 1. full s1, full s2
* 2. MISMATCH_LENGTH 3p (right most) bases looking for a mismatch of s1,
*    and s2
* 3. full s1, full reverse s2 (not the complement)
* 4. MISMATCH_LENGTH 3p (right most) bases looking for a mismatch of s1,
*    and rev s2
* 5. MISMATCH_LENGTH 3p (right most) bases looking for a mismatch of rev s1,
*    and s2.  This is needed for symmetry
* return 1 if greater hamming distances are greater than `min_dist` and at
least one mismatch occurs at the 3p end
* return 0 if any hamming distances are too small
*/
unsigned char ss_hammingAnd3pMismatchCheck(const char *s1, const char *s2,
                                const int len, const int min_dist);


/*
* Compute the hamming distance forward and reverse of 2 strings of length `len`
* and return 1 if both hamming distances are greater than min_dist
* otherwise return 0
*/
int ss_hamming2XCheck(const char *s1, const char *s2, int len, const int min_dist);


// Lookup table for DNA complementarity
extern const char COMP_BASE_LUT[128];

#endif