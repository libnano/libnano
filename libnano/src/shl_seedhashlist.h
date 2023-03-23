#ifndef SEEDHASHLIST_H
#define SEEDHASHLIST_H

#include "khash.h"
// #include "kvec.h"

// #define SH_DEBUG

typedef struct {
    int idx;        /* index in the pool of sequences */
    int offset;     /* index into the sequence */
} tuple_t;

typedef struct {
    size_t n; tuple_t* arr;
} vectuple_t;

// initialize the khash table type required
KHASH_MAP_INIT_STR(str2vectuple, vectuple_t*);

typedef struct {
    khash_t(str2vectuple)* ht;  // the hash table
    vectuple_t* vecbuffer;      // array of subarrays pointing to the idx_list of a given hash key
    char* all_hashes;           // a buffer containing all hash keys
    tuple_t* idx_buffer;        // the buffer containing the by-key sorted indices of potential occurences

    int* seed_idxs;     // non-joker indices of seeds
    int seed_idxs_len;  // AKA seed weight
    int seed_len;       // length of seed mask
} seed_search_obj_t;

/*
return the weight of a seed (the number of #'s)
*/
int shl_seed_weight(char* seed, int seed_len);

/*
Convert a seed to an integer array of indices of the #'s
*/
int shl_parse_seed(char* seed, int seed_len, int* seed_idxs);

/*
free memory of seed_search_obj
*/
int shl_free_hash_obj(seed_search_obj_t* seed_search_obj);

/*
Build a hash table for a given seed against a list of strings in reference
*/
int shl_buildSeedTable(
    int* seed_idxs,         // Array of non-wildcard indices in the seed
    int seed_idxs_len,      // Length of the seed_idxs array
    int seed_len,           // Full length of the seed
    char* reference,        // Library of sequences concatenated for which to build hash table
    int reference_count,    // Length of reference seq char array
    int* reference_lengths, // length of each element in the library of sequences
    int* reference_idxs,    // start offset, and length for next sequence in reference
    int reference_idxs_len, // len of the reference_idxs array equal to number of sequences
    seed_search_obj_t* seed_search_obj); // Dict pointer to be populated w/
                                        // lookups based on seed hashes
int checkht(seed_search_obj_t* seed_search_obj);
/*
Find all hits with less than mismatches mismatches.
mismatches should be less the the k mismatches the seed the seed_search_obj_t
was built with
*/
int shl_findSingleMatches(   int mismatches,                    // number of mismatches allowed
                            const char* target,                 // Target seq for which to build hash table
                            int target_len,                     // the length of the m-mer AKA m
                            char* reference,                    // the reference string (genome)
                            int reference_count,
                            int* reference_lengths,             // length of each element in the library of sequences
                            int* reference_idxs,                // start offset, and length for next sequence in referenc
                            seed_search_obj_t* seed_search_obj, // Dict pointer
                            tuple_t** matches,
                            int buffer_size);                   // size of the match buffe to be allocated

// int shl_findAllMatches(  int mismatches,
//                         int target_len,     // the length of the m-mer AKA m
//                         char* reference,    // the reference string (genome)
//                         int reference_len,
//                         int buffer_size,
//                         seed_search_obj_t* seed_search_obj,     // Dict pointer
//                         int** unique_counts);
#endif
