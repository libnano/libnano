#ifndef NANO_H
#define NANO_H

#include <Python.h>

static inline void* ikcalloc(size_t N, size_t Z) {
    void* p;
    if ( (p = PyMem_Malloc((N) * (Z))) ) {
        memset(p, 0, (N)*(Z));
        return p;
    }
    return NULL;
}

#define kcalloc(N,Z) ikcalloc((N), (Z))

#define kmalloc(Z) PyMem_Malloc((Z))

#define krealloc(P,Z) PyMem_Realloc((P), (Z))

#define kfree(P) PyMem_Free((void*) P)

#endif
