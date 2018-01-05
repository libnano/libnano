from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cpython cimport Py_XINCREF, Py_XDECREF
from cpython.ref cimport PyObject

from libc.string cimport memmove

DEF _INIT_VEC_CAP = 32

# Py_ssize_t
cdef struct FeatureInstance:
    Py_ssize_t start
    Py_ssize_t end
    int feature_type
    PyObject* obj

cdef class FeatureVector:
    cdef:
        FeatureInstance* data
        size_t n, capacity 

    def __cinit__(self):
        self.n = self.capacity = 0
        self.data = NULL
        self.resize(_INIT_VEC_CAP)
    # end def

    def __dealloc__(self):
        cdef int i
        if self.data != NULL:
            for i in range(self.n):
                Py_XDECREF(self.data[i].obj)
            PyMem_Free(self.data)
    # end def

    def append(self, o):
        cdef FeatureInstance fi
        fi.start = o.location[1]
        fi.end = o.location[2]
        fi.feature_type = 0
        fi.obj = <PyObject *> o
        self.push(fi)
    # end def

    cdef inline FeatureInstance* take(self, int i):
        """ return a pointer
        """
        return &self.data[i]
    # end def

    cdef inline resize(self, size_t s):
        cdef FeatureInstance* temp = <FeatureInstance*> PyMem_Realloc(self.data,
                                sizeof(FeatureInstance)*s)
        if temp != NULL:
            self.data = temp
            self.capacity = s
        else:
            #PyMem_Free(self.data)
            raise MemoryError("Out of memory")

    cdef inline pop(self, FeatureInstance* out):
        """ pop off the end of the list
        """
        if self.n != 0:
            self.n -= 1
            out[0] = self.data[self.n]
            Py_XDECREF(out.obj)
        else:
            raise ValueError("OB")

    def __len__(self):
        return self.n

    cdef inline size_t size(self):
        return self.n

    cdef inline push(self, FeatureInstance x):
        """ assign by value
        """
        cdef FeatureInstance* temp
        cdef size_t cap = self.capacity
        if self.n == cap:
            cap = cap << 1 if cap != 0 else 2
            temp = <FeatureInstance*> PyMem_Realloc(self.data,
                                sizeof(FeatureInstance) * cap)
            if temp != NULL:
                self.capacity = cap
                self.data = temp
            else:
                #PyMem_Free(self.data)
                raise MemoryError("Out of memory")
        self.data[self.n] = x
        Py_XINCREF(x.obj)
        self.n += 1
    # end def

    cdef inline insertAt(self, FeatureInstance x, int idx):
            """ At most resizes by one
            """
            cdef size_t elements_to_move = self.n - idx
            cdef FeatureInstance* fvec = self.data
            if x.obj == fvec[idx].obj:
                raise ValueError("Feature already exists at %d" % (idx))
            if self.n == self.capacity:
                self.resize(self.capacity << 1)
            if idx < self.capacity:
                if idx < self.n:
                    memmove(<void*> &fvec[idx + 1], 
                            <void*> &fvec[idx], 
                            sizeof(FeatureInstance)*elements_to_move)
                    fvec[idx] = x
                    Py_XINCREF(x.obj)
                elif idx == self.n:
                    fvec[idx] = x
                    Py_XINCREF(x.obj)
                else:
                    raise IndexError("index %d out of element range %d",
                            idx, self.n)
            else:
                raise IndexError("index %d out of capacity range %d", 
                        idx, self.capacity)
    # end def

    cdef inline int bisect_left(self, FeatureInstance x) except -1:
        """ Locate the insertion point for x in a to maintain sorted order.
        """
        cdef int i = 0
        cdef int start = x.start
        cdef int end = x.end
        cdef FeatureInstance* fvec = self.data

        cdef int lo = 0
        cdef int hi = self.n

        while lo < hi:
            mid = (lo + hi) // 2
            if fvec[mid].start < start:
                lo = mid + 1
            else:
                hi = mid

        # now find the 3 prime in case of features with the same start
        for i in range(lo, self.n):
            if end < fvec[i].end:
                continue
            else:
                lo = i

        return lo
    # end def

    cdef index(self, FeatureInstance x):
        'Locate the leftmost value exactly equal to x'
        i = self.bisect_left(x)
        if i != self.n:
            while self.data[i].start == x.start:
                if self.data[i].obj == x.obj:
                    return i
        raise ValueError("index: FeatureInstance not in FeatureVector")

    cdef inline int remove_inst(self, FeatureInstance x) except -1:
        cdef int idx = self.index(x)
        cdef size_t elements_to_move = self.n - idx - 1
        memmove(<void*> &self.data[idx], 
                <void*> &self.data[idx + 1], 
                sizeof(FeatureInstance)*elements_to_move)
        Py_XDECREF(x.obj)
        self.n -= 1
        return 0
    # end def

    def remove(self, o):
        cdef FeatureInstance fi
        fi.start = o.location[1]
        fi.end = o.location[2]
        fi.feature_type = 0
        fi.obj = <PyObject *> o
        self.remove_inst(fi)
    # end def

    def insert(self, o):
        cdef FeatureInstance fi
        fi.start = o.location[1]
        fi.end = o.location[2]
        fi.feature_type = 0
        fi.obj = <PyObject *> o
        self.insort_left(fi)
    # end def

    cdef inline insort_left(self, FeatureInstance x):
        cdef int i
        cdef int lo = 0
        cdef int hi = self.n
        cdef int start = x.start
        cdef int end = x.end
        cdef FeatureInstance* fvec = self.data

        # find insertion of 5 prime
        while lo < hi:
            mid = (lo + hi) // 2

            if fvec[mid].start < start:
                lo = mid + 1
            else:
                hi = mid

        # now find the 3 prime in case of features with the same start
        for i in range(lo, self.n):
            if end < fvec[i].end:
                continue
            else:
                lo = i
        self.insertAt(x, lo)
    # end def