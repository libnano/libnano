# cython: embedsignature=True

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memset
from libc.string cimport memcpy

"""
modified to be more like standard lib malloc, calloc, and realloc

branched from https://github.com/syllog1sm/cymem but
see: https://github.com/explosion/cymem for the current version
"""

cdef class Pool:
    """Track allocated memory addresses, and free them all when the Pool is
    garbage collected.  This provides an easy way to avoid memory leaks, and
    removes the need for deallocation functions for complicated structs.

    >>> from cymem.cymem cimport Pool
    >>> cdef Pool mem = Pool()
    >>> data1 = <int*> mem.calloc(10, sizeof(int))
    >>> data2 = <float*> mem.calloc(12, sizeof(float))

    Attributes:
        size (size_t): The current size (in bytes) allocated by the pool.
        addresses (dict): The currently allocated addresses and their sizes. Read-only.
    """

    def __cinit__(self):
        self.size = 0
        self.addresses = {}

    def __dealloc__(self):
        cdef size_t addr
        for addr in self.addresses:
            PyMem_Free(<void*>addr)
        self.size = 0
        self.addresses = {}

    cdef void* malloc(self, size_t number, size_t elem_size) except NULL:
        """Allocate a 0-initialized block_size-byte block of memory, and
        remember its address. The block will be freed when the Pool is garbage
        collected.
        """
        cdef void* p = NULL
        cdef size_t block_size = number*elem_size
        p = PyMem_Malloc(block_size)
        if p == NULL:
            raise MemoryError("Couldn't allocate")
        self.addresses[<size_t>p] = block_size
        self.size += block_size
        return p

    cdef void* calloc(self, size_t number, size_t elem_size) except NULL:
        """Allocate a 0-initialized number*elem_size-byte block of memory, and
        remember its address. The block will be freed when the Pool is garbage
        collected.
        """
        cdef void* p = NULL
        cdef size_t block_size = number*elem_size
        p = PyMem_Malloc(block_size)
        if p == NULL:
            raise MemoryError("Couldn't allocate")
        memset(p, 0, block_size)
        self.addresses[<size_t>p] = block_size
        self.size += block_size
        return p

    cdef void* realloc(self, void* p, size_t new_size) except NULL:
        """Resizes the memory block pointed to by p to new_size bytes, returning
        a non-NULL pointer to the new block. new_size must be larger than the
        original.

        If p is not in the Pool or new_size is 0, a MemoryError is raised.

        This differs from normal Cymem in that it doesn't always memcpy when
        resizing
        """
        cdef void* new = NULL
        cdef size_t old_size

        if <size_t>p not in self.addresses:
            raise MemoryError("Pointer %d not found in Pool %s" % (<size_t>p, self.addresses))
        if new_size == 0:
            raise MemoryError("Realloc requires new_size != 0")

        old_size = self.addresses[<size_t>p]

        new = PyMem_Realloc(p, new_size)
        if new == NULL:
            # __dealloc__ should handle freeing the memory in p
            raise MemoryError("Realloc failed")

        self.size -= old_size
        self.size += new_size

        # only need to pop address if there as a change
        if new != p:
            self.addresses.pop(<size_t>p)

        self.addresses[<size_t>new] = new_size
        return new

    cdef void free(self, void* p) except *:
        """Frees the memory block pointed to by p, which must have been returned
        by a previous call to Pool.alloc.  You don't necessarily need to free
        memory addresses manually --- you can instead let the Pool be garbage
        collected, at which point all the memory will be freed.

        If p is not in Pool.addresses, a KeyError is raised.
        """
        self.size -= self.addresses.pop(<size_t>p)
        PyMem_Free(p)

    cdef void own(self, void* p, size_t number, size_t elem_size) except *:
        """ Take ownership of a pointer not allocated by Pool
        New to libnano's version of cymem
        """
        cdef size_t block_size = number*elem_size
        if <size_t>p in self.addresses:
            raise MemoryError("Pointer %d already owned" % (<size_t>p))
        elif p == NULL:
            raise MemoryError("Can't own NULL")
        else:
            self.addresses[<size_t>p] = block_size
            self.size += block_size
    # end def

    cdef void disown(self, void* p) except *:
        """ Like free but no actually memory freeing, allows for a another
        item to take over the memory
        """
        if <size_t>p not in self.addresses:
            raise MemoryError("Pointer %d not found in Pool %s" % (<size_t>p, self.addresses))
        self.size -= self.addresses.pop(<size_t>p)
    # end def

cdef class Address:
    """A block of number * size-bytes of 0-initialized memory, tied to a Python
    ref-counted object. When the object is garbage collected, the memory is freed.

    >>> from cymem.cymem cimport Address
    >>> cdef Address address = Address(10, sizeof(double))
    >>> d10 = <double*>address.ptr

    Args:
        number (size_t): The number of elements in the memory block.
        elem_size (size_t): The size of each element.

    Attributes:
        ptr (void*): Pointer to the memory block.
        addr (size_t): Read-only size_t cast of the pointer.
    """
    def __cinit__(self, size_t number, size_t elem_size):
        self.ptr = NULL

    def __init__(self, size_t number, size_t elem_size):
        self.ptr = PyMem_Malloc(number * elem_size)
        memset(<void*> self.ptr, 0, number * elem_size)

    property addr:
        def __get__(self):
            return <size_t>self.ptr

    def __dealloc__(self):
        PyMem_Free(self.ptr)
