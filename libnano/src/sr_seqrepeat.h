/******************************************************************************
** seqrepeat.h
**
******************************************************************************/

#ifndef SEQ_REPEAT_H
#define SEQ_REPEAT_H
#include <inttypes.h>
#include "nano.h"
#include "khash.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Data structure representing the root of one or more repeated sequences
 * sharing a common base word (of the initial word_size used to
 * buildRepeatData)
 */
typedef struct repeatroot_t {
    int count;                  // the total number of occurances of this seqint
    int start;                  // the first index of an occurance of this seqint
    int last;                   // the next index of an occurance of this seqint
                                // (last may no longer be needed)
} repeatroot_t;
KHASH_MAP_INIT_INT64(map_repeatroot, repeatroot_t*)

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The main data structure for repeat analysis. Created and populated in
 * sr_buildRepeatData.
 */
typedef struct repeatcheck_t {
    repeatroot_t *repeatroot_arr;    // the root array of length 4^word_size
    int word_size;              // the root word size
    int seq_length;             // the length of the input sequence
    int count_max;              // the maximum number of occurences of any seqint
    // repeat_idx_arrs - array of index arrays where greater than zero entries
    // indicate a start of a repeat at that index and the value is the index
    // of the next occurence of that repeat.
    // repeat_idx_arrs[k][0] would be the beginning of the array of repeat
    // starts of seeds of length word_size+k
    // an entry of -1 in an index array means that this is the last repeat in
    // the sequence
    // [0, 4, 0, 0, 9, 0, 0, 0, 0, -1, 0]
    // is an example of an index array with one repeat occuring 3 times
    // since the next start of a repeat can never be index 0, there is no problem
    // reserving zero for a non-start
    int *repeat_idx_arr;
    khash_t(map_repeatroot) *ht_root;
    // uint64_t* seqint_arr;       // the converted seqint array of a input sequence
} repeatcheck_t;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Builds the base repeat data structures -- the `repeatroot_arr` (comprised
 * of repeatroot_t objects for each possible seqint), and the `repeat_idx_arrs`
 * which denote indices for each word size at which a repeated sequence is
 * found (the value of the element at that index is described above).
 */
int sr_buildRepeatData(char* seq,                // char array of DNA sequence
                       int seq_length,           // length of char array
                       int word_size,           // repeat root word size
                       // int max_word_size,        // max repeat size to calc
                       repeatcheck_t** rcheck
                    );
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Populates arrays denoting extended repeats for word sizes greater than the
 * root word size (root_word_size), up to and including the max_word_size
 */
// int sr_buildExtendedRepeats(repeatcheck_t* rcheck);

// Helper function to free repeatcheck_t structs
int sr_freeRepeatCheck(repeatcheck_t* rcheck);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Builds an array of seqints representing the `word_size` sequence beginning
 * at each index in `seq`. Allocates memory for `seqint_array` internally.
 */
// int sr_buildSeqintArray(char* seq, int seq_length, int word_size,
//                         uint64_t** seqint_array);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Generates an array representing a "pileup" of repeat data from
 * root_word_size to max_word_size. The value at each index in the resulting
 * array will be the sum of the number of different size repeats that the
 * underlying base participates in, from root_word_size to max_word_size. This
 * is a simple heuristic for identifying potential repeat "hotspots".
 */
// int sr_repeatPileup(repeatcheck_t* rcheck, int min_word_size, int max_word_size,
//                  int *pileup_array);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Performs a rolling window analysis to look for local regions with high
 * repeated sequence content.
 */
int sr_repeatWindow(repeatcheck_t* rcheck,
    int window_size,
    int max_repeat_count,
    int** repeat_violation_idxs,
    int *violation_count,
    int** repeat_counts);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Find the index of a query string if it exists in seq
 * uses the seed as a shortcut
 */
int sr_queryRepeatIdx(repeatcheck_t* rcheck, char* seq, char* query, int query_length);
#endif
