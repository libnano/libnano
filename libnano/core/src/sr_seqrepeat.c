/******************************************************************************
** seqrepeat.c
**
******************************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> /* for NULL pointers */
#include <assert.h>
#include "sr_seqrepeat.h"
#include "si_seqint.h"

// #include "nano.h"
// #include "khash.h"

// typedef struct {
//     size_t n; int* arr;
// } vecint_t;

// KHASH_MAP_INIT_INT64(map_vecint_t, (vecint_t* ))

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Free a repeatcheck_t struct
 */
int sr_freeRepeatCheck(repeatcheck_t* rcheck) {
    // does not free the rcheck struct
    if (rcheck->repeat_idx_arr != NULL) {
        kfree(rcheck->repeat_idx_arr);
        rcheck->repeat_idx_arr = NULL;
    }
    if (rcheck->repeatroot_arr != NULL) {
        kfree(rcheck->repeatroot_arr);
        rcheck->repeatroot_arr = NULL;
    }
    if (rcheck->ht_root != NULL) {
        kh_destroy(map_repeatroot, rcheck->ht_root);
        rcheck->ht_root = NULL;
    }
    return 0;
};


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Initialize a repeatcheck_t struct for `seq` of length `seq_length`. The
 * initial repeat data structures will be generated using the provided
 * `word_size` (i.e., the minimal measurable repeat unit in bp).
 */
int sr_buildRepeatData( char* seq, int seq_length,
                        int word_size, repeatcheck_t** rcheck) {

    uint64_t        seq_int;
    char            *seq_ptr = seq + word_size;
    khint_t k;
    int                 i, il, ret;
    repeatroot_t        *cur_root; // pointer to root structure representing
                                   // the current repeat
    int                 count_max = 0;  // the max occurences of any one seq.

    // Array of repeatroot_t structs
    repeatroot_t* repeatroot_arr = NULL;
    repeatroot_t* repeatroot_arr_ptr = NULL;

    // Hash table for storing counts of root words and start and end
    // indices
    khash_t(map_repeatroot) *ht_root = NULL;
    ht_root = kh_init(map_repeatroot);
    assert(ht_root != NULL);


    /* Array indexed by position in the original sequence. Element values:
     *       0:  the word starting at the current index is not repeated
     *       >0: the word starting at the current position is repeated and
     *           the next occurance is at the index indicated by the value
     *       -1: the word starting at the current position is repeated but
     *           the current instance is the last occurance in the sequence
     */
    int* repeat_idx_array = NULL;

    // pointers get effectively NULLed by calloc
    // http://stackoverflow.com/questions/5857588/calloc-with-structure-with-pointers-in-c
    repeatroot_arr = (repeatroot_t*) kcalloc(seq_length - word_size + 1,
                                            sizeof(repeatroot_t));
    if (repeatroot_arr == NULL) {
        goto brd_fail;
    }
    repeatroot_arr_ptr = repeatroot_arr;

    repeat_idx_array = (int *) kcalloc(seq_length, sizeof(int));
    if (repeat_idx_array == NULL) {
        goto brd_fail;
    }

    seq_int = si_seq2Int(seq, word_size);
    k = kh_put(map_repeatroot, ht_root, seq_int, &ret);
    if (ret == 1) {
        cur_root = repeatroot_arr_ptr++;
        kh_value(ht_root, k) = cur_root;  // assign a pointer
        cur_root->start = 0;
        cur_root->count++;
        count_max = 1;
        cur_root->last = 0;
    } else {
        printf("failed building ht: %d\n", ret);
        goto brd_fail;
    }
    for (i = 1, il = seq_length - word_size + 1; i < il; i++) {
        // Shift base into seqint while shifting the first base out
        seq_int = si_addToWindow(seq_int, *seq_ptr++, word_size);
        // search for hash pointer
        k = kh_get(map_repeatroot, ht_root, seq_int);   // query the hash table
        if (k == kh_end(ht_root)) {
            // it's missing so install it in the table and set start
            k = kh_put(map_repeatroot, ht_root, seq_int, &ret);
            if (ret == 1) {
                cur_root = repeatroot_arr_ptr++;
                kh_value(ht_root, k) = cur_root;  // assign a pointer
                cur_root->start = i;
            } else {
                printf("failed building ht: %d\n", ret);
                goto brd_fail;
            }
        } else {
            cur_root = kh_value(ht_root, k);
            // set the element at the index of the previous instance of this
            // repeat to the index of the current instance
            repeat_idx_array[cur_root->last] = i;
            // -1 is a terminator of a repeats sequence
            repeat_idx_array[i] = -1;
        }
        cur_root->count++;
        if (cur_root->count > count_max) {
            count_max++;
        }
        cur_root->last = i;
    }

    // create the repeatcheck_t structure and populate
    (*rcheck)->word_size = word_size;
    (*rcheck)->seq_length = seq_length;
    (*rcheck)->repeatroot_arr = repeatroot_arr;
    (*rcheck)->repeat_idx_arr = repeat_idx_array;
    (*rcheck)->count_max = count_max;
    (*rcheck)->ht_root = ht_root;
    return 1;
brd_fail:
    if (repeatroot_arr != NULL) {
        kfree(repeatroot_arr);
    }
    if (ht_root != NULL) {
        kh_destroy(map_repeatroot, ht_root);
    }
    if (repeat_idx_array != NULL) {
        kfree(repeat_idx_array);
    }
    sr_freeRepeatCheck(*rcheck);
    return -1;
};

int sr_queryRepeatIdx(repeatcheck_t* rcheck, char* seq,
                        char* query, int query_length) {
    int idx;
    int word_size  = rcheck->word_size;
    if (query_length < word_size){
        return -1;
    }
    // capture just the root word size
    khash_t(map_repeatroot) *ht_root = rcheck->ht_root;
    uint64_t seq_int = si_seq2Int(query, word_size);
    khint_t k = kh_get(map_repeatroot, ht_root, seq_int);   // query the hash table
    if (k == kh_end(ht_root)) {
        return -1;
    }
    repeatroot_t *rroot = kh_value(ht_root, k);
    if (rroot->count > 0) {
        idx = rroot->start;
    } else {
        return -1;
    }

    if (query_length == word_size) {
        return idx;
    } else {
        int* arr_base = rcheck->repeat_idx_arr;
        int seq_length = rcheck->seq_length;
        while (idx < seq_length) {
            if (strncmp(seq + idx, query, query_length) == 0) {
                return idx;
            }
            idx = arr_base[idx];
            if (idx < 0) {
                return idx;
            }
        }
    }
    return -1;
}

int sr_repeatWindow(repeatcheck_t* rcheck,
                    int window_size,
                    int max_repeat_count,
                    int** repeat_violation_idxs,
                    int *violation_count,
                    int** repeat_counts) {

    /*
    window_size - size of the rolling window
    max_repeat_count - the count threshold not to exceed
    repeat_violation_idxs - the limit violating indices
    violation_count - the length of the violating indices
    repeat_counts - the counts of repeats all told
    */

    int     i, j, test_val, jl;
    int seq_length = rcheck->seq_length;

    int il = seq_length - window_size + 1;
    int word_size = rcheck->word_size;

    int did_alloc = 1;

    if (window_size > seq_length) {
       return -1;
    }

    // initialize things for pointer return values

    int* rpt_counts = NULL;
    // repeat_counts can be a NUMPY array pointer passed in
    if (*repeat_counts == NULL) {
        rpt_counts = (int *) kcalloc(il, sizeof(int));
        if (rpt_counts == NULL) {
            return -1;
        }
    } else {
        did_alloc = 0;
        rpt_counts = *repeat_counts;
    }

    // printf("did_alloc: %d, %d\n", did_alloc, il);

    int vcount = 0;
    // over allocate the buffer
    int* violation_idxs = (int *) kcalloc(il, sizeof(int));
    if (violation_idxs == NULL) {
        if (did_alloc) {
            kfree(rpt_counts);
        }
        return -1;
    }

    int *repeat_idx_array = rcheck->repeat_idx_arr;
    if (repeat_idx_array == NULL) {
        printf("Error repeat index array is NULL\n");
    }

    for (i=0; i < il; i++) {

        for (j=i, jl=i + window_size - word_size + 1; j < jl; j++)  {
            test_val = repeat_idx_array[j];
            if ( (test_val > 0) && (test_val < jl ) ) {
                rpt_counts[i] += 1;
            }
        }
    }
    for (i=0; i < il; i++) {
        if (rpt_counts[i] > max_repeat_count) {
            violation_idxs[vcount++] = i;
        }
    }
    *repeat_violation_idxs = violation_idxs;
    if (did_alloc) {
        *repeat_counts = rpt_counts;
    }
    *violation_count = vcount;
    return 1;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


