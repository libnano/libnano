#ifndef NANO_KHASH_H
#define NANO_KHASH_H

#include "nano.h"
#include "khash.h"

// type 1.
typedef struct {
    int idx;        /* index in the pool of sequences */
    int offset;     /* index into the sequence */
} tuple_t;

typedef struct {
    size_t n; tuple_t* arr;
} vectuple_t;


KHASH_MAP_INIT_STR(str2vectuple, vectuple_t*)

// type 2.

#endif
