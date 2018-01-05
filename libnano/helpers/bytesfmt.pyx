cdef extern from "_bytesfmt.h":
    object PyBytes_Format(object fmt, object args)

def bformat(object fmt, object args):
    return PyBytes_Format(fmt, args)