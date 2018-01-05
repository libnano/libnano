cdef class Pool:
    cdef readonly size_t size
    cdef readonly dict addresses

    cdef void* malloc(self, size_t number, size_t elem_size) except NULL
    cdef void* calloc(self, size_t number, size_t size) except NULL
    cdef void* realloc(self, void* addr, size_t n) except NULL
    cdef void free(self, void* addr) except *
    cdef void own(self, void* p, size_t number, size_t elem_size) except *
    cdef void disown(self, void* p) except *


cdef class Address:
    cdef void* ptr