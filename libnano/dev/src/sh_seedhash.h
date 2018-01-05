#ifndef SEEDHASH_H
#define SEEDHASH_H

#include "khash.h"
// #include "kvec.h"

// #define SH_DEBUG

typedef struct {
    size_t n; int* arr;
} vecint_t;


KHASH_MAP_INIT_STR(str2vec, vecint_t*);

// commented out kvec versionfields
// typedef kvec_t(int) vecint2;
// KHASH_MAP_INIT_STR(str2vec2, vecint2*);

typedef struct {
    // commented out kvec versionfields
    // khash_t(str2vec2) *ht;
    // vecint2 * vecbuffer; // array of subarrays pointing to the idx_list of a given hash key
 
    khash_t(str2vec) * ht; // the hash table
    vecint_t * vecbuffer; // array of subarrays pointing to the idx_list of a given hash key
    char * all_hashes;  // a buffer containing all hash keys
    int * idx_buffer;     // the buffer containing the by-key sorted indices of potential occurences

    int* seed_idxs;     // non-joker indices of seeds
    int seed_idxs_len;  // AKA seed weight
    int seed_len;       // length of seed mask
} seed_table_obj_t;

int sh_free_hash_obj(seed_table_obj_t* seed_table_obj);

int sh_buildSeedTable(int* seed_idxs,   // Array of non-wildcard indices in the seed
                int seed_idxs_len,  	// Length of the seed_idxs array
                int seed_len,       	// Full length of the seed
                char* reference,       	// Target seq for which to build hash table
                int reference_len,     	// Length of reference seq char array
                seed_table_obj_t** seed_table_obj);


int sh_findSingleMatches(  int mismatches,
                            char* target,       // Target seq for which to build hash table
                            int target_len,     // Length of target seq char array
                            char* reference,    // the reference string (genome)
                            int reference_len,
                            int buffer_size, 
                            seed_table_obj_t* seed_table_obj, 
                            int** matches);

int sh_findAllMatches(  int mismatches,
                        int target_len,     // the length of the m-mer AKA m
                        char* reference,    // the reference string (genome)
                        int reference_len,
                        int buffer_size, 
                        seed_table_obj_t* seed_table_obj,     // Dict pointer
                        int** unique_counts);
#endif