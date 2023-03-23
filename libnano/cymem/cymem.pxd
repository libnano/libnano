cdef class Pool:
    cdef:
        readonly size_t size
        readonly dict addresses

        void* malloc(self, size_t number, size_t elem_size) except NULL
        void* calloc(self, size_t number, size_t size) except NULL
        void* realloc(self, void* addr, size_t n) except NULL
        void free(self, void* addr) except *
        void own(self, void* p, size_t number, size_t elem_size) except *
        void disown(self, void* p) except *


cdef class Address:
    cdef void* ptr
