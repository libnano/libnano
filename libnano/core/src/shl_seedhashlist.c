#include <string.h> // for null pointer
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "nano.h"
#include "ss_seqstr.h"
#include "shl_seedhashlist.h"

#include "ksort.h"
#include "pairedsort.h"

// initialize klib style sort stuff
#define const_lt_str(a, b) (strcmp(((const char *) a), ((const char *) b)) < 0)
PAIREDSORT_INIT(keyslave, char*, tuple_t, const_lt_str);
#define lt_tup(a, b) ( ((a).idx < (b).idx) ? 1 : ( ((a).idx != (b).idx) ? 0 : ((a).offset < (b).offset) ) )

KSORT_INIT(tuple_t, tuple_t, lt_tup);

// #define SHL_DEBUG

int shl_seed_weight(char* seed,
                    int seed_len) {
    int i;
    int w = 0;
    char* seed_ptr = seed;
    for (i=0; i < seed_len; i++, seed_ptr++) {
        if (*seed_ptr == '#') {
            ++w;
        }
    }
    return w;
}

int shl_parse_seed(char* seed,
                    int seed_len,
                    int* seed_idxs) {
    int i;
    char* seed_ptr = seed;
    int* seed_idx_ptr = seed_idxs;
    for (i=0; i < seed_len; i++, seed_ptr++) {
        if (*seed_ptr == '#') {
            *seed_idx_ptr++ = i;
        }
    }
    return 0;
}


int shl_free_hash_obj(seed_search_obj_t* seed_search_obj) {
    if (seed_search_obj != NULL) {
        if (seed_search_obj->idx_buffer != NULL) {
            kfree(seed_search_obj->idx_buffer);
            seed_search_obj->idx_buffer = NULL;
        }
        if (seed_search_obj->all_hashes != NULL) {
            kfree(seed_search_obj->all_hashes);
            seed_search_obj->all_hashes = NULL;
        }
        if (seed_search_obj->vecbuffer != NULL) {
            kfree(seed_search_obj->vecbuffer);
            seed_search_obj->vecbuffer = NULL;
        }
        // COMMENTED OUT. seed_idxs not allocated here
        // if (seed_search_obj->seed_idxs != NULL) {
        //     kfree(seed_search_obj->seed_idxs);
        //     seed_search_obj->seed_idxs = NULL;
        // }
        if (seed_search_obj->ht != NULL) {
            kh_destroy(str2vectuple, seed_search_obj->ht);
            seed_search_obj->ht = NULL;
        }
        return 0;
    } else {
        return -1;
    }
}

int shl_buildSeedTable(
    int* seed_idxs,         // Array of non-wildcard indices in the seed
    int seed_idxs_len,      // Length of the seed_idxs array
    int seed_len,           // Full length of the seed
    char* reference,        // Library of sequences concatenated for which to build hash table
    int reference_count,    // Length of reference seq char array
    int* reference_lengths, // length of each element in the library of sequences
    int* reference_idxs,    // start offset, and length for next sequence in reference
    int reference_idxs_len, // len of the reference_idxs array equal to number of sequences
    seed_search_obj_t* seed_search_obj) // Dict pointer to be populated w/
                            // lookups based on seed hashes
{
    int ret;
    int i, j, q;
    int local_num_m_mers;
    int idx_offset;

    int* idx_ptr = NULL;

    // 1. Setup constants
    const int* idx_ptr_lim = seed_idxs + seed_idxs_len;

    // const int num_m_mers = reference_len - seed_len + 1; // for one
    // subtract off length of null terminators, then the shifted seed lengths
    const int num_m_mers = (reference_count - reference_idxs_len) + \
                             (1 - seed_len)*reference_idxs_len;

    const size_t hash_key_size = seed_idxs_len + 1; // add 1 for 0 terminator required by khash

    // 2. Allocate temporaries
    char *hash_key = (char*) kmalloc((hash_key_size)*sizeof(char));
    hash_key[seed_idxs_len] = '\0';    // 0 terminate string

    char* hash_key_ptr = NULL;

    vectuple_t* vb_temp = NULL;
    char* ah_temp = NULL;

    // 3. declare HEAP vars
    // buffer containing one copy of each key concatenated
    char* all_hashes = NULL;
    // buffer containing all of the indexes to be sorted by key
    tuple_t* idx_buffer = NULL;
    // sortable pointer array to the buffer of keys
    char** all_hashes_sortable = NULL;

    // buffer contains the values in the hash table
    vectuple_t * vecbuffer = NULL;
    // hash table herself
    khash_t(str2vectuple) *ht = NULL;

    khint_t k; // local key


    // 4. allocate memory
    // num_m_mers could be big so put on HEAP
    all_hashes_sortable = (char **) kmalloc(num_m_mers*sizeof(char *));
    assert(all_hashes_sortable != NULL);

    char** ahs_ptr =  all_hashes_sortable;

    all_hashes = (char *) kmalloc((seed_idxs_len+1)*num_m_mers*sizeof(int));
    assert(all_hashes != NULL);

    char* all_hashes_ptr = all_hashes;

    idx_buffer = (tuple_t *) kmalloc(num_m_mers*sizeof(tuple_t));
    assert(idx_buffer != NULL);

    tuple_t* idx_buffer_ptr = idx_buffer;


    vecbuffer = (vectuple_t *) kcalloc(num_m_mers, sizeof(vectuple_t));
    assert(vecbuffer != NULL);

    vectuple_t* vecbuffer_ptr = vecbuffer;

    ht = kh_init(str2vectuple);
    assert(ht != NULL);

    // 5.  Hash all reference relative to string.
    q = 0;
    for (i=0; i < reference_idxs_len; i++) {
        idx_offset = reference_idxs[i];
        local_num_m_mers = reference_lengths[i] - seed_len + 1;
        for (j=0; j < local_num_m_mers; j++) {
            for (idx_ptr = seed_idxs, hash_key_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
                *hash_key_ptr++ = reference[idx_offset + j + *idx_ptr++];
            };
            // search for hash pointer
            k = kh_get(str2vectuple, ht, hash_key);   // query the hash table
            if (k == kh_end(ht)) {
                // it's missing so install it in the table
                memcpy(all_hashes_ptr, hash_key, hash_key_size);

                k = kh_put(str2vectuple, ht, all_hashes_ptr, &ret);
                if (ret == 1) {
                    kh_value(ht, k) = vecbuffer_ptr++;  // assign a pointer
                } else {
                    printf("failed building ht: %d\n", ret);
                    goto bt_fail;
                }
                all_hashes_ptr += hash_key_size;
            }
            kh_value(ht, k)->n++; // increment the number of indices found for this key

            idx_buffer_ptr->idx = i;
            idx_buffer_ptr->offset = j;
            idx_buffer_ptr++;

            all_hashes_sortable[q++] = (char *) kh_key(ht, k);  // retrieve pointer to key
        }
    }
    assert(q == num_m_mers);

    // 6. sort both keys and idx_buffer by key.  The sort needs to be stable
    // i.e. preserves relative order of idx_buffer elements with equivalent key strings
    idx_buffer_ptr = idx_buffer;
    // specialMergesort(ahs_ptr, idx_buffer_ptr, num_m_mers);
    mergesort_paired_keyslave(num_m_mers, ahs_ptr, idx_buffer_ptr);

    // 7. put in the hashtable
    // 7a. do iteration zero to set ahs_ptr to non-NULL
    ahs_ptr = all_hashes_sortable;
    idx_buffer_ptr = idx_buffer;
    k = kh_get(str2vectuple, ht, *ahs_ptr);
    if (k == kh_end(ht)) {
        printf("DOUBLE CRAAPPPP1\n");
        goto bt_fail;
    }
    vecbuffer_ptr = kh_value(ht, k);
    vecbuffer_ptr->arr = idx_buffer_ptr;    // assign pointer
    idx_buffer_ptr += vecbuffer_ptr->n;

    // 7b. now do the rest of the iterations
    for (i=1; i < num_m_mers; i++) {
        // test to see if we have encountered a new key in the sortable list
        if (*ahs_ptr != all_hashes_sortable[i]) {
            ahs_ptr = &(all_hashes_sortable[i]);
            k = kh_get(str2vectuple, ht, *ahs_ptr);
            if (k == kh_end(ht)) {
                printf("DOUBLE CRAAPPPP2%s\n", *ahs_ptr);
                goto bt_fail;
            } else {
                vecbuffer_ptr = kh_value(ht, k);
                vecbuffer_ptr->arr = idx_buffer_ptr;    // assign pointer to index buff
                idx_buffer_ptr += vecbuffer_ptr->n;
            }
        }
    }

    // 8. realloc memory of buffer for actual number of keys and assign to output struct
    // We can't trust realloc to leave memory locations put so we must
    // add the memory offset delta wherever we are storing pointers
    seed_search_obj->ht = ht;
    khint_t num_keys = kh_size(ht);

    vb_temp = (vectuple_t *) krealloc(vecbuffer, num_keys*sizeof(vectuple_t));
    assert(vb_temp != NULL);
    // you can subtract two pointers of the same type.
    Py_ssize_t vb_delta = vb_temp - vecbuffer;
    seed_search_obj->vecbuffer = vb_temp;

    ah_temp = (char *) krealloc(all_hashes, (seed_idxs_len+1)*num_keys*sizeof(char));
    assert(ah_temp != NULL);
    // you can subtract two pointers of the same type.
    Py_ssize_t ah_delta =  ah_temp -  all_hashes;
    seed_search_obj->all_hashes = ah_temp;

    for (k = kh_begin(ht); k != kh_end(ht); ++k) {
        if (kh_exist(ht, k)) {
            kh_value(ht, k)  = kh_value(ht, k) + vb_delta;
            kh_key(ht, k) = kh_key(ht, k) + ah_delta;
        }
    }


    // 9. copy pointers to output data structure
    seed_search_obj->idx_buffer = idx_buffer;
    seed_search_obj->seed_idxs = seed_idxs;
    seed_search_obj->seed_idxs_len = seed_idxs_len;
    seed_search_obj->seed_len = seed_len;
    // assign to output pointer
    kfree(hash_key);
    kfree(all_hashes_sortable);
    // printf("built\n");

    return 0;

// failure label
bt_fail:
    kfree(hash_key);
    kfree(all_hashes_sortable);
    kfree(vecbuffer);
    kfree(all_hashes);
    kfree(idx_buffer);
    kh_destroy(str2vectuple, ht);
    return -1;
}

int checkht(seed_search_obj_t* seed_search_obj) {
    khash_t(str2vectuple) *ht = seed_search_obj->ht;
    if (kh_get(str2vectuple, ht, "TCTGAC") != kh_end(ht)) {
        printf("::::::::: in the table\n");
    }
    printf("HT point C: %p, %p, %d\n", ht, ht->keys, (int) kh_size(ht));
    return 0;
}

int shl_findSingleMatches(  int mismatches,
                            const char* target,                 // Target seq for which to build hash table
                            int target_len,                     // the length of the m-mer AKA m
                            char* reference,                    // the reference string (genome)
                            int reference_count,
                            int* reference_lengths,             // length of each element in the library of sequences
                            int* reference_idxs,                // start offset, and length for next sequence in referenc
                            seed_search_obj_t* seed_search_obj, // Dict pointer
                            tuple_t** matches,                  // assumed to be NULL
                            int buffer_size)                    // size of the match buffe to be allocated

{
    // assumes a non-bogus input (seed_len != 1, not a reference of all A's etc.)
    int i, j;
    khiter_t k;

#ifdef SHL_DEBUG
    char debugbuff[target_len+1];
    debugbuff[target_len] = '\0';
#endif

    const int* seed_idxs = seed_search_obj->seed_idxs;           // Array of non-wildcard indices in the seed
    const int seed_idxs_len = seed_search_obj->seed_idxs_len;    // Length of the seed_idxs array
    const int seed_len = seed_search_obj->seed_len;              // Full length of the seed

    const char* target_local = target;

    int offset;
    int offset_last = -1;   // used for pre deduping filtering
    int ref_idx;
    int ref_idx_last = -1;

    const int* idx_ptr;
    const int *idx_ptr_lim = seed_idxs + seed_idxs_len;

    char *hash_key = (char*) kmalloc((seed_idxs_len + 1)*sizeof(char));

    hash_key[seed_idxs_len] = '\0';    // 0 terminate string
    char* hash_ptr = NULL;

    vectuple_t *hash_vec = NULL;
    khash_t(str2vectuple) *ht = seed_search_obj->ht;

    // heuristic still, no obvious way to predict how big this buffer needs to be other than
    // the full reference_len
    int do_free_matches_for_fail = 1;
    assert(*matches == NULL);
    tuple_t* match_buffer = NULL;
    match_buffer = (tuple_t *) kmalloc(buffer_size*sizeof(tuple_t));
    assert(match_buffer != NULL);
    if (match_buffer == NULL) { goto fsm_fail; }

    tuple_t* match_ptr = match_buffer;
    int match_count = 0;

    // const int max_offset = reference_len - target_len;
    // collect all potential matches into match_buffer
    const int hashes_per_m_mer = target_len - seed_len + 1;
    for (i=0; i < hashes_per_m_mer; i++) {
        for (idx_ptr = seed_idxs, hash_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
            *hash_ptr++ = target_local[i + *idx_ptr++];
        }
        k = kh_get(str2vectuple, ht, hash_key);   // query the hash table

        if (k == kh_end(ht)) {
            // goto fsm_fail;
            continue;
            // not there
        } else {
            hash_vec = (vectuple_t *) kh_value(ht, k);
            const int jl = (int) (hash_vec->n);
            for (j=0; j < jl;) {
                ref_idx = hash_vec->arr[j].idx;
                offset =  hash_vec->arr[j++].offset - i;

                if ( ( offset != offset_last ) || (ref_idx != ref_idx_last) ) {
                    // Add integer index of current hash to the end of list in reference_table
                    // printf("$ref_idx of match: %d\n", ref_idx);
                    match_ptr->idx = ref_idx_last = ref_idx;
                    match_ptr->offset = offset_last = offset;
                    match_ptr++;
                    match_count++;
                }
            }
        }
    }
    // sort match_buffer
    ks_introsort(tuple_t, match_count, match_buffer);

    // dedupe the match buffer and do the Hamming distance as necessary
    offset_last = -1;
    ref_idx_last = -1;
    int repeat_count = 0; // count of the uniques of each index comprising a total count
    int dist;
    int start_idx, byte_offset, match_length;
    const int must_match = target_len - mismatches;

    const tuple_t *match_ptr_lim = match_buffer + match_count;
    for (match_ptr=match_buffer; match_ptr < match_ptr_lim;) {
        ref_idx = match_ptr->idx;
        offset = match_ptr->offset;
        match_ptr++;
        if ((offset != offset_last) || (ref_idx != ref_idx_last)) {
            start_idx = reference_idxs[ref_idx];
            match_length = reference_lengths[ref_idx];
            byte_offset = start_idx + offset;
            // printf("* check %d, %d, %d, %d\n", ref_idx, offset, target_len, match_length);
            if (offset > 0) {
                // printf("check %d, %d, %d\n", offset, target_len, match_length);
                if ((offset + target_len) < match_length) {;
                    dist = ss_hamming(  target_local,
                                    (const char*) reference + byte_offset,
                                    target_len);
                } else { // case where at the right end of a string
                    if ((match_length - offset) < must_match) {
                        dist = mismatches + 1;  // can't match
                        // printf("doop %d %d %d\n", must_match, match_length, offset);
                    } else {  // limit the hamming distance
                        dist = ss_hamming(  target_local,
                                        (const char*) reference + byte_offset,
                                        match_length - offset);
                    }
                }
            } else { // left edge of string
                dist = ss_hamming(  target_local - offset,
                                    (const char*) reference + start_idx,
                                    target_len + offset);
            }

            if (dist <= mismatches) {
                ref_idx_last = match_buffer[repeat_count].idx = ref_idx;
                offset_last = match_buffer[repeat_count].offset = offset;
                repeat_count++;
                #ifdef SHL_DEBUG
                printf("hit %d, %d at ref_idx: %d, start_idx: %d, offset: %d\n", dist, target_len, ref_idx, start_idx, offset);
                // printf("%s\n", reference+start_idx+offset);
                #endif
            }
#ifdef SHL_DEBUG
            else {
                printf("false+ %d, %d at ref_idx: %d, start_idx: %d, offset: %d\n", dist, target_len, ref_idx, start_idx, offset);
                strncpy(debugbuff, reference+ref_idx+offset, target_len);
                printf("%s %s\n", target_local, debugbuff);
            }
#endif
        }
    }

    // prepare for return
    if (do_free_matches_for_fail) {
        if (repeat_count > 0) {
            *matches = (tuple_t *) krealloc(match_buffer, sizeof(tuple_t) * repeat_count);
            if (*matches == NULL) {
                kfree(match_buffer);
                return -1;
            }
        } else {
            kfree(match_buffer);
        }
    }
    kfree(hash_key);
#ifdef SHL_DEBUG
    // kfree(debugbuff);
#endif
    return repeat_count;

fsm_fail:
    kfree(hash_key);
#ifdef SHL_DEBUG
    // kfree(debugbuff);
#endif
    if (do_free_matches_for_fail) {
        kfree(match_buffer);
    }
    return -1;
}
