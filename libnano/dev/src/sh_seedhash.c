#include <string.h> // for null pointer
#include <stdlib.h>
#include <stdio.h>
#include "ss_seqstr.h"
#include "sh_seedhash.h"

#include "ksort.h"
#include "pairedsort.h"
// #include <x86intrin.h>

// initialize klib style sort stuff
#define const_lt_str(a, b) (strcmp(((const char *) a), ((const char *) b)) < 0)
PAIREDSORT_INIT(keyval, char*, int, const_lt_str);
KSORT_INIT_GENERIC(int);

int sh_free_hash_obj(seed_table_obj_t* seed_table_obj) {
    free(seed_table_obj->idx_buffer);
    free(seed_table_obj->all_hashes);
    free(seed_table_obj->vecbuffer);
    free(seed_table_obj->seed_idxs);
    kh_destroy(str2vec, seed_table_obj->ht);
    free(seed_table_obj);
    return 0;
}

int sh_buildSeedTable(int* seed_idxs,  // Array of non-wildcard indices in the seed
                int seed_idxs_len,     // Length of the seed_idxs array
                int seed_len,          // Full length of the seed
                char* reference,       // Target seq for which to build hash table
                int reference_len,     // Length of reference seq char array
                seed_table_obj_t** seed_table_obj) // Dict pointer to be populated w/  
                                        // lookups based on seed hashes
{
    int ret;
    int i;
    int *idx_ptr;
    const int *idx_ptr_lim = seed_idxs + seed_idxs_len;

    const int num_m_mers = reference_len - seed_len + 1;   

    const size_t hash_key_size = seed_idxs_len + 1; // add 1 for 0 terminator required by khash 
#ifdef _WIN32
    char *hash_key = (char*) malloc((hash_key_size)*sizeof(char));
#else
    char hash_key[hash_key_size];
#endif
    hash_key[seed_idxs_len] = '\0';    // 0 terminate string
    
    char* hash_key_ptr;

    // declare HEAP vars
    // buffer containing one copy of each key concatenated
    char* all_hashes = NULL;
    // buffer containing all of the indexes to be sorted by key
    int* idx_buffer = NULL;
    // sortable pointer array to the buffer of keys
    char** all_hashes_sortable = NULL;

    // buffer contains the values in the hash table
    vecint_t * vecbuffer = NULL;
    // hash table herself
    khash_t(str2vec) *ht = NULL;
    khint_t k; // local key

    all_hashes_sortable = (char **) malloc(num_m_mers*sizeof(char *)); // num_m_mers could be big so put on HEAP
    if (all_hashes_sortable == NULL) {
        goto bt_fail;
    }
    char** ahs_ptr =  all_hashes_sortable;

    all_hashes = (char *) malloc((seed_idxs_len+1)*num_m_mers*sizeof(int));
    if (all_hashes == NULL) {
        goto bt_fail;
    }
    char * all_hashes_ptr = all_hashes;

    idx_buffer = (int *) malloc(num_m_mers*sizeof(int));
    if (idx_buffer == NULL) {
        goto bt_fail;
    }
    int* idx_buffer_ptr = idx_buffer;

    
    vecbuffer = (vecint_t *) calloc(num_m_mers, sizeof(vecint_t));
    if (vecbuffer == NULL) {
        goto bt_fail;
    }

    vecint_t* vecbuffer_ptr = vecbuffer;
    
    ht = kh_init(str2vec);
    if (ht == NULL) {
        goto bt_fail;
    }
    // COMMENTED OUT, resizing up front is slow since we don't know the exact size
    // if (kh_resize(str2vec, ht, n) < 0) {
    //     goto bt_fail;
    // }

    // 1.  Hash all reference relative to string.
    for (i=0; i < num_m_mers; i++) {
        for (idx_ptr = seed_idxs, hash_key_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
            *hash_key_ptr++ = reference[i + *idx_ptr++];
        }
        // search for hash pointer
        k = kh_get(str2vec, ht, hash_key);   // query the hash table
        if (k == kh_end(ht)) {
            // it's missing so install it in the table   
            strcpy(all_hashes_ptr, (const char *) hash_key);
            // printf("%s %s %s %s\n", all_hashes_ptr, hash_key, &reference[reference_len-1], &reference[n-1]); fflush(stdout);
            k = kh_put(str2vec, ht, all_hashes_ptr, &ret);
            if (ret == 1) {
                kh_value(ht, k) = vecbuffer_ptr++;  // assign a pointer
            } else {
                printf("failed building ht: %d\n", ret);
                goto bt_fail;
            }
            all_hashes_ptr += hash_key_size;
        }
        kh_value(ht, k)->n++; // increment the number of indices found for this key
        *idx_buffer_ptr++ = i;
        all_hashes_sortable[i] = (char *) kh_key(ht, k);  // retrieve pointer to key
    }

    // 2. sort both keys and idx_buffer by key.  The sort needs to be stable
    // i.e. preserves relative order of idx_buffer elements with equivalent key strings
    idx_buffer_ptr = idx_buffer;
    // specialMergesort(ahs_ptr, idx_buffer_ptr, num_m_mers);
    mergesort_paired_keyval(num_m_mers, ahs_ptr, idx_buffer_ptr);

    // 3. put in the hashtable
    // 3a. do iteration zero to set ahs_ptr to non-NULL
    ahs_ptr = all_hashes_sortable;
    idx_buffer_ptr = idx_buffer;
    k = kh_get(str2vec, ht, *ahs_ptr);
    if (k == kh_end(ht)) {
        goto bt_fail;
    }
    vecbuffer_ptr = kh_value(ht, k);
    vecbuffer_ptr->arr = idx_buffer_ptr;    // assign pointer
    idx_buffer_ptr += vecbuffer_ptr->n;
    // 3b. now do the rest of the iterations
    for (i=1; i < num_m_mers; i++) {
        // test to see if we have encountered a new key in the sortable list
        if (*ahs_ptr != all_hashes_sortable[i]) {
            ahs_ptr = &(all_hashes_sortable[i]);
            k = kh_get(str2vec, ht, *ahs_ptr);
            if (k == kh_end(ht)) {
                goto bt_fail;
            } else {
                vecbuffer_ptr = kh_value(ht, k);
                vecbuffer_ptr->arr = idx_buffer_ptr;    // assign pointer to index buff
                idx_buffer_ptr += vecbuffer_ptr->n;
            }
        }
    }

    //4. realloc memory of buffer for actual number of keys
    seed_table_obj_t* out = NULL;
    out = (seed_table_obj_t *) malloc(sizeof(seed_table_obj_t));
    if (out == NULL) {
        goto bt_fail;
    }
    out->ht = ht;
    khint_t num_keys = kh_size(ht);
    // COMMENTED OUT, resizing at end is slow since we don't know the exact size
    // may make sense to try another permutation of counting keys but that's TBD
    // if (kh_resize(str2vec, ht, num_keys) < 0) {
    //     goto bt_fail;
    // }
    if ( (out->vecbuffer = (vecint_t *) realloc(vecbuffer, num_keys*sizeof(vecint_t))) == NULL) {
        goto bt_fail;
    }
    if ( (out->all_hashes = (char *) realloc(all_hashes, (seed_idxs_len+1)*num_keys*sizeof(char))) == NULL) {
        goto bt_fail;
    }
    // 5. copy pointers to output data structure
    out->idx_buffer = idx_buffer;
    out->seed_idxs = seed_idxs;
    out->seed_idxs_len = seed_idxs_len;
    out->seed_len = seed_len;
    // assign to output pointer
    *seed_table_obj = out;
#ifdef _WIN32
    free(hash_key);
#endif
    free(all_hashes_sortable);
    // printf("built\n");
    return 0;

    bt_fail:
#ifdef _WIN32
        free(hash_key);
#endif
        free(all_hashes_sortable);
        free(vecbuffer);
        free(all_hashes);
        free(idx_buffer);
        free(seed_idxs);
        kh_destroy(str2vec, ht);
        return -1;
}

int sh_findSingleMatches(   int mismatches,
                            char* target,       // Target seq for which to build hash table
                            int target_len,     // the length of the m-mer AKA m
                            char* reference,    // the reference string (genome)
                            int reference_len,
                            int buffer_size, 
                            seed_table_obj_t* seed_table_obj,     // Dict pointer
                            int** matches)  
{
    // assumes a non-bogus input (seed_len != 1, not a reference of all A's etc.)
    int i;
    int j;
    khiter_t k;
    
#ifdef SH_DEBUG
    char debugbuff[target_len+1];
    debugbuff[target_len] = '\0';
#endif
    
    const int* seed_idxs = seed_table_obj->seed_idxs;           // Array of non-wildcard indices in the seed
    const int seed_idxs_len = seed_table_obj->seed_idxs_len;    // Length of the seed_idxs array
    const int seed_len = seed_table_obj->seed_len;              // Full length of the seed

    int offset;
    int offset_last = -1;   // used for pre deduping filtering
    const int* idx_ptr;
    const int *idx_ptr_lim = seed_idxs + seed_idxs_len;

    
#ifdef _WIN32
    char *hash_key = (char*) malloc((seed_idxs_len + 1)*sizeof(char));
#else
    char hash_key[seed_idxs_len + 1];
#endif
    hash_key[seed_idxs_len] = '\0';    // 0 terminate string
    char* hash_ptr;
    
    vecint_t *hash_vec;
    khash_t(str2vec) *ht = seed_table_obj->ht;

    // heuristic still, no obvious way to predict how big this buffer needs to be other than 
    // the full reference_len
    int do_free_matches_for_fail = 1;
    int* match_buffer = NULL;
    if (*matches == NULL) {
        int *match_buffer = (int *) malloc(buffer_size*sizeof(int));
        if (match_buffer == NULL) {
            return -1;
        }
    } else {
        do_free_matches_for_fail = 0;
        match_buffer = *matches;    // matches needs to be buffer_size
    }


    int* match_ptr = match_buffer;
    int match_count = 0;

    // const int max_offset = reference_len - target_len;
    // collect all potential matches into match_buffer
    const int hashes_per_m_mer = target_len - seed_len + 1;
    for (i=0; i < hashes_per_m_mer; i++) {
        for (idx_ptr = seed_idxs, hash_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
            *hash_ptr++ = target[i + *idx_ptr++];
        }

        k = kh_get(str2vec, ht, hash_key);   // query the hash table

        if (k == kh_end(ht)) {
            goto fsm_fail;
            // not there
        } else {
            hash_vec = (vecint_t *) kh_value(ht, k);
            const int jl = (int)hash_vec->n; 
            for (j=0; j < jl;) {
                offset = (int) (hash_vec->arr[j++] - i);
                // Following if else if is for edge cases at beginning and end of 
                // the reference sequence
                if (offset < 0) {
                    offset = 0;     // clip at zero is an edge case
                } 
                if ( offset != offset_last ) {
                    // Add integer index of current hash to the end of list in reference_table
                    *match_ptr++ = offset_last = offset;
                    match_count++;
                }
            }
        }
    }
    // sort match_buffer
    ks_introsort(int, match_count, match_buffer);

    // dedupe the match buffer and do the Hamming distance as necessary
    offset_last = -1;
    int repeat_count = 0; // count of the uniques of each index comprising a total count
    int dist;

    const int *match_ptr_lim = match_buffer + match_count;
    for (match_ptr=match_buffer; match_ptr < match_ptr_lim;) {
        if ((offset = *match_ptr++) != offset_last) {
            dist = ss_hamming((const char *) target, (const char *) reference+offset, target_len);
            if (dist <= mismatches) {
                offset_last = match_buffer[repeat_count++] = offset;
                #ifdef SH_DEBUG
                // printf("hit %d, %d at %d\n", dist, target_len, offset);
                #endif
            } 
#ifdef SH_DEBUG
            else {
                // printf("false+ %d, %d at index %d\n", dist, target_len, offset);
                strncpy(debugbuff, reference+offset, target_len);
                printf("%s %s\n", target, debugbuff);
            }
#endif
        }
    }

    // prepare for return 
    if (do_free_matches_for_fail) {
        if (repeat_count > 0) {
            *matches = (int *) realloc(match_buffer, sizeof(int) * repeat_count);
            if (*matches == NULL) {
                free(match_buffer);
                return -1;
            }
        } else {
            free(match_buffer);
        }
    } 
    // printf("repeat count is %d\n", repeat_count);
#ifdef _WIN32
    free(hash_key);
#ifdef SH_DEBUG
    free(debugbuff);
#endif
#endif
    return repeat_count;

fsm_fail:
#ifdef _WIN32
    free(hash_key);
#ifdef SH_DEBUG
    free(debugbuff);
#endif
#endif
    if (do_free_matches_for_fail) {
        free(match_buffer);
    }
    return -1;
}

int sh_findAllMatches(  int mismatches,
                        int target_len,     // the length of the m-mer AKA m
                        char* reference,    // the reference string (genome)
                        int reference_len,
                        int buffer_size, 
                        seed_table_obj_t* seed_table_obj,     // Dict pointer
                        int** repeat_counts                   // pointer to array which will contain the repeat count at each index
                        )  {
    // assumes a non-bogus input (seed_len != 1, not a reference of all A's etc.)
    int i,j, q;
    khiter_t k;
    int repeat_count_max = 0;
#ifdef SH_DEBUG
#ifdef _WIN32
    char* debugbuff = (char*) malloc((target_len+1)*sizeof(char));
#else
    char debugbuff[target_len+1];
#endif
    debugbuff[target_len] = '\0';
#endif
    
    const int* seed_idxs = seed_table_obj->seed_idxs;           // Array of non-wildcard indices in the seed
    const int seed_idxs_len = seed_table_obj->seed_idxs_len;    // Length of the seed_idxs array
    const int seed_len = seed_table_obj->seed_len;              // Full length of the seed
    const int num_m_mers = reference_len - seed_len + 1; 

    const int* idx_ptr;
    const int *idx_ptr_lim = seed_idxs + seed_idxs_len;

    
#ifdef _WIN32
    char* hash_key = (char*) malloc((seed_idxs_len + 1)*sizeof(char));
#else
    char hash_key[seed_idxs_len + 1];
#endif
    hash_key[seed_idxs_len] = '\0';    // 0 terminate string
    char* hash_ptr;
    
    vecint_t *hash_vec;
    khash_t(str2vec) *ht = seed_table_obj->ht;

    // heuristic still, no obvious way to predict how big this buffer needs to be other than 
    // the full reference_len

    const int hashes_per_m_mer = target_len - seed_len + 1;

    // reverse complement index delta for a given seed to get the position of a rc hit
    const int rc_idx_delta = target_len - seed_len;

    // collect pointers to potential matches a buffer since the i+1 m-mer uses
    // the hashes_per_m_mer - 1 of the hits of the i-th m-mer
#ifdef _WIN32
    vecint_t **active_hash_vals_buffer = (vecint_t**) malloc(hashes_per_m_mer*sizeof(vecint_t*));   // the array of indices of the end of a match per hash
    vecint_t **active_hash_vals_rc_buffer = (vecint_t**) malloc(hashes_per_m_mer*sizeof(vecint_t*));
#else
    vecint_t *active_hash_vals_buffer[hashes_per_m_mer];   // the array of indices of the end of a match per hash
    vecint_t *active_hash_vals_rc_buffer[hashes_per_m_mer];
#endif

    int ahv_idx = 0;                       // keep track of where we are at in the buffer of potential indice hits

    int *match_buffer = NULL;
    int *match_rc_buffer = NULL;
    int *out = NULL;
    int do_free_out_for_fail = 1;

    // tracking variables
    int match_count = 0;
    int match_rc_count = 0;
    int offset;
    int offset_last = -1;   // used for pre deduping filtering
    int offset_rc;
    int offset_rc_last = -1;   // used for pre deduping filtering

#ifdef _WIN32
    char* target_rc_buffer = NULL;
    char* rc_buffer = NULL;
    rc_buffer = (char*) malloc((target_len+1)*sizeof(char));
    target_rc_buffer = (char*) malloc((target_len+1)*sizeof(char));
#else
    char rc_buffer[target_len+1];
    char target_rc_buffer[target_len+1];
#endif
    rc_buffer[target_len] = '\0';
    target_rc_buffer[target_len] = '\0';

    // 1. preamble fill up ahv pointer buffer
    // difference between the start of an m-mer and the last index of the last hash of it
    const int index_diff = hashes_per_m_mer - 1;
    
    if (*repeat_counts == NULL) {
        printf("mallocing\n"); fflush(stdout);
        int* out = (int *) malloc(num_m_mers*sizeof(int)); 
        if (out == NULL) {
            repeat_count_max = -1;
            goto fam_cleanup;
        }
    } else {
        do_free_out_for_fail = 0;
        out = *repeat_counts;
    }

    match_buffer = (int *) malloc(buffer_size*sizeof(int));
    if (match_buffer == NULL) {
        repeat_count_max = -1;
        goto fam_cleanup;
    }
    int *match_ptr = match_buffer;

    match_rc_buffer = (int *) malloc(buffer_size*sizeof(int));
    if (match_rc_buffer == NULL) {
        repeat_count_max = -1;
        goto fam_cleanup;
    }
    int *match_rc_ptr = match_rc_buffer;

    for (i = 0; i < index_diff;) {
        // A. get the reverse complement first
        const char * rev_ref_root = reference+i+seed_len-1;
        for (idx_ptr = seed_idxs, hash_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
                *hash_ptr++ = COMP_BASE_LUT[(unsigned char) *(rev_ref_root - *idx_ptr++)];
        }
        k = kh_get(str2vec, ht, hash_key);   // query the hash table
        if (k != kh_end(ht)) {      // the reverse comp does not have to be there
            // printf("non-null\n"); fflush(stdout);
            active_hash_vals_rc_buffer[i] = (vecint_t *) kh_value(ht, k);
        } else {
            active_hash_vals_rc_buffer[i] = NULL;
        }

        // B. get the forward now and increment i as required
        for (idx_ptr = seed_idxs, hash_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
                *hash_ptr++ = reference[i + *idx_ptr++];
        }
        k = kh_get(str2vec, ht, hash_key);   // query the hash table
        if (k == kh_end(ht)) {      // the forward has to be there
            repeat_count_max = -1;
            goto fam_cleanup;
            // not there
        } else {
            active_hash_vals_buffer[i++] = (vecint_t *) kh_value(ht, k);
        }
    }
    // 2. now collect the ith hash value and generate the match set for
    // the i - hashes_per_m_mer;
    ahv_idx = i;
    // ahv_idx should be equal to index_diff
    // const int limit_idx = reference_len - seed_len + 1;
    const int limit_idx = reference_len - target_len + 1;
    for (; i < limit_idx; i++) {
        // A. get the reverse complement first
        match_rc_count = 0;
        match_rc_ptr = match_rc_buffer; // reset the match buffer
        const char * rev_ref_root = reference+i+seed_len-1;
        for (idx_ptr = seed_idxs, hash_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
                *hash_ptr++ = COMP_BASE_LUT[(unsigned char) *(rev_ref_root - *idx_ptr++)];
        }
        k = kh_get(str2vec, ht, hash_key);   // query the hash table
        if (k != kh_end(ht)) {      // the reverse comp does not have to be there
            // printf("non-null\n"); fflush(stdout);
            active_hash_vals_rc_buffer[ahv_idx] = (vecint_t *) kh_value(ht, k);
        } else {
            active_hash_vals_rc_buffer[ahv_idx] = NULL;
        }

        // B. get the forward now
        for (idx_ptr = seed_idxs, hash_ptr=hash_key; idx_ptr < idx_ptr_lim;) {
            *hash_ptr++ = reference[i + *idx_ptr++];
        }

        k = kh_get(str2vec, ht, hash_key);   // query the hash table
        if (k == kh_end(ht)) {
            repeat_count_max = -1;
            goto fam_cleanup;
            // not there
        } else {
            int iter_ahv_idx;
            // ahv_idx is currently pointing to the last hash of a m-mer 
            active_hash_vals_buffer[ahv_idx++] = (vecint_t *) kh_value(ht, k);

            // A. now copy all reverse comp potential hits to the match_buffer for sorting
            iter_ahv_idx = ahv_idx %= hashes_per_m_mer;
            // ahv_idx is now pointing to the first hash of a m-mer
            for (q = 0; q < hashes_per_m_mer; q++) {
                hash_vec = active_hash_vals_rc_buffer[iter_ahv_idx++ % hashes_per_m_mer];
                if (hash_vec != NULL) {
                    const int jl = (int)hash_vec->n;
                    for (j = 0; j < jl;) {
                        // same offset calculation for rc as we'll be calculating the actual m-mer RC in step 3
                        offset_rc = (int) (hash_vec->arr[j++] - rc_idx_delta + q);
                        // Following if else if is for edge cases at beginning and end of 
                        // the reference sequence
                        if (offset_rc < 0) {
                            offset_rc = 0;     // clip at zero is an edge case
                        } 
                        if ( offset_rc != offset_rc_last ) {
                            // Add integer index of current hash to the end of list in reference_table
                            *match_rc_ptr++ = offset_rc_last = offset_rc;
                            match_rc_count++;
                        }
                    } // end for j
                } // end NULL check
            } // end for q

            // B. now copy all forward potential hits to the match_buffer for sorting  
            iter_ahv_idx = ahv_idx;
            // ahv_idx is now pointing to the first hash of a m-mer 
            match_count = 0;
            match_ptr = match_buffer; // reset the match buffer
            for (q = 0; q < hashes_per_m_mer; q++) {
                hash_vec = active_hash_vals_buffer[iter_ahv_idx++ % hashes_per_m_mer];
                const int jl = (int)hash_vec->n;
                for (j = 0; j < jl;) {
                    offset = (int) (hash_vec->arr[j++] - q);
                    // Following if else if is for edge cases at beginning and end of 
                    // the reference sequence
                    if (offset < 0) {
                        offset = 0;     // clip at zero is an edge case
                    } 
                    if ( offset != offset_last ) {
                        // Add integer index of current hash to the end of list in reference_table
                        *match_ptr++ = offset_last = offset;
                        match_count++;
                    }
                } // end for j
            } // end for q
        } // end else
        
        if ((match_count + match_rc_count) > num_m_mers) {
            repeat_count_max = -1;
            goto fam_cleanup;
        }

        // sort match_buffer
        ks_introsort(int, match_rc_count, match_rc_buffer);
        ks_introsort(int, match_count, match_buffer);

        // dedupe the match buffer and do the Hamming distance as necessary
        offset_last = -1;
        int repeat_count = 0;       // count of the uniques of each index comprising a total count
        int dist;

        const int *match_ptr_lim = match_buffer + match_count;
        const int repeat_idx = i - index_diff;
        const char* target = reference + repeat_idx;
        for (match_ptr=match_buffer; match_ptr < match_ptr_lim;) {
            if ((offset = *match_ptr++) != offset_last) {
                dist = ss_hamming(target, (const char *) reference + offset, target_len);
                if (dist <= mismatches) {
                    offset_last = match_buffer[repeat_count++] = offset;
                    #ifdef SH_DEBUG
                    printf("hit %d, %d at %d\n", dist, target_len, offset);
                    #endif
                } 
            }
        }
        // Reverse complement
        const int *match_rc_ptr_lim = match_rc_buffer + match_rc_count;
        ss_revCompSeqCopy(target_rc_buffer, target, target_len);
        int repeat_rc_count = 0;
        offset_rc_last = -1;
        for (match_rc_ptr = match_rc_buffer; match_rc_ptr < match_rc_ptr_lim;) {
            if ((offset_rc = *match_rc_ptr++) != offset_rc_last) {
                ss_revCompSeqCopy(rc_buffer, reference + offset_rc, target_len);
                dist = ss_hamming(target_rc_buffer, (const char *) rc_buffer, target_len);
                if (dist <= mismatches) {
                    offset_rc_last = match_rc_buffer[repeat_rc_count++] = offset_rc;
                    #ifdef SH_DEBUG
                    printf("rc hit %d, %d at %d\n", dist, target_len, offset_rc);
                    #endif
                } 
            }
        }
        repeat_count += repeat_rc_count;
        out[repeat_idx] = repeat_count;
        if (repeat_count > repeat_count_max) { repeat_count_max = repeat_count; }
    } // end for i
    if (do_free_out_for_fail) {
        *repeat_counts = out;
    }
fam_cleanup:
#ifdef _WIN32
    free(rc_buffer);
    free(target_rc_buffer);
    free(hash_key);
    free(active_hash_vals_buffer);
    free(active_hash_vals_rc_buffer);
#ifdef SH_DEBUG
    free(debugbuff);
#endif
#endif
    free(match_buffer);
    free(match_rc_buffer);
    if (repeat_count_max < 0) {
        if (do_free_out_for_fail) { free(out); }
        return -1;
    } else {
        return repeat_count_max;

    }
}