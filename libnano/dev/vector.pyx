from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from libc.string cimport memset
from libc.string cimport memcpy

cdef inline fused_t v_init(fused_vec_t v, size_t i):
    v.n = v.capacity = 0
    v.data = NULL
# end def

cdef inline fused_t v_take(fused_vec_t v, size_t i):
    return  v.data[i]
# end def

cdef inline v_resize(fused_vec_t v, size_t s, size_t tsize) except -1:
    cdef void* temp = <void*> PyMem_Realloc(v.data, tsize*(v.m))
    if temp != NULL:
        v.data = temp
        v.capacity = s
    else:
        PyMem_Free(v.data)
        raise MemoryError("Out of memory")

cdef inline fused_t v_pop(fused_vec_t v):
    v.n -= 1
    return v.data[v.n]  # v.n is now the last element

cdef inline size_t size(fused_vec_t v):
    return v.n

cdef inline v_push(fused_vec_t v, fused_t x, size_t tsize) except -1:
    cdef void* temp
    if v.n == v.capacity:
        v.capacity = v.capacity << 1 if v.capacity != 0 else 2
        temp = <void*> PyMem_Realloc(v.data, tsize * v.capacity)
        if temp != NULL:
            <void*> v.data = temp
        else:
            PyMem_Free(v.data)
            raise MemoryError("Out of memory")
    v.n += 1
    v.data[v.n] = x
# end def

cdef inline destroy(fused_vec_t v):
    PyMem_Free(v.data)

cdef inline v_insert(fused_vec_t v, fused_t x, size_t index, size_t tsize) except -1:
        """ At most resizes by one
        """
        cdef size_t elements_to_move = v.n - index
        if v.n == v.capacity:
            v_resize(v, v.capacity + 1, tsize)
        if index < v.capacity:
            if index < v.n:
                memcpy(<void*> &v.data[index+1], <void*> &v.data[index], tsize*elements_to_move)
                v.data[index] = x
            elif index == v.n:
                v.data[index] = x
            else:
                raise IndexError("index %d out of element range %d", index, v.n)
        else:
            raise IndexError("index %d out of capacity range %d", index, v.capacity)

cdef class ObjectVector:
    cdef pyobject_vec_t ovec

    def __cinit__(self, size_t m):
        v_init(self.ovec, m)

    cdef inline push_c(self, object obj):
        # stealing a reference means not calling Py_XDECREF
        v_push(self.ovec, obj, sizeof(PyObject*))

    cdef inline object pop_c(self, object obj):
        return v_pop(self.ovec)

    cdef inline insert_c(self, object obj, size_t index):
        return v_insert(self.ovec, obj, index, sizeof(PyObject*))

    def __dealloc__(self):
        cdef size_t i;
        for i in range(self.ovec.n):
            Py_XDECREF((<PyObject*>self.ovec[i])
        destroy(self.ovec)  