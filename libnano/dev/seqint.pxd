"""
seqintx.pxd
~~~~~~~~~~~~~~~~~~

Cython header file for seqintx.pyx -- allows for cross-project Cython /
C integration of the low-level thermodynamic analysis bindings.

This header mainly serves to help wrap functions written in raw c
that are wrapped in seqintx.pyx such that those functions need not
be recompiled nor the associated headers included
"""
from libnano.helpers.inttypes cimport uint64_t
cdef inline int2Seq_c(uint64_t, char*, int)


