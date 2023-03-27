# cython: language_level=3, boundscheck=False, wraparound=False
# Copyright (C) 2014-2018. Nick Conway & Ben Pruitt; Wyss Institute
# Copyright (C) 2023 Nick Conway & Ben Pruitt;
# See LICENSE.TXT for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
libnano.seqgraph
~~~~~~~~~~~~~~~~

Graph operations on sets of DNA, RNA, and amino acid sequences

'''
from typing import (
    Generator,
    List,
)

cimport cython

import numpy as np

cimport numpy as cnp

cnp.import_array()

cimport c_util

__doc__ = "Graph operations on sets of DNA, RNA, and amino acid sequences"

cdef extern from "ss_seqstr.h":
    int ss_hamming(const char*, const char*, int)
    int ss_minHamming(const char*, const char*, int, int)
    unsigned char ss_hammingAnd3pMismatchCheck(const char *, const char*, int, int)
    int ss_minHammingThreshold(const char*, const char*, int, int, int)
    int ss_hamming2XCheck(const char*, const char*, int, const int)

# structure for hamming distance
cdef struct hamming_node_t:
    cnp.int32_t distance, index

cdef struct hamming_node_3p_t:
    cnp.int32_t distance, mismatch3p, index


@cython.boundscheck(True)
@cython.wraparound(True)
def _find_cliques(
        cnp.uint8_t[:,:] data_graph,
) -> Generator[List[int]]:
    """Adapted from networkx 1.9.1 clique.py BSD license

    Search for all maximal cliques in a graph.
    Maximal cliques are the largest complete subgraph containing
    a given node.  The largest maximal clique is sometimes called
    the maximum clique.


    NOTE: for some reason I can't type data as a buffer type due to some function
    local variable thing in cython that I don't fully get why this is different

    Uses cython memory views for args:
        http://docs.cython.org/src/userguide/memoryviews.html

    Args:
        data_graph (numpy 2D ndarray): undirected graph representing the
            connected system

    Returns:
        generator of lists: generator of member list for each maximal clique

    See Also
    --------
    find_cliques_recursive :
    A recursive version of the same algorithm
    Notes
    -----
    To obtain a list of cliques, use list(find_cliques(G)).
    Based on the algorithm published by Bron & Kerbosch (1973) [1]_
    as adapted by Tomita, Tanaka and Takahashi (2006) [2]_
    and discussed in Cazals and Karande (2008) [3]_.
    The method essentially unrolls the recursion used in
    the references to avoid issues of recursion stack depth.
    This algorithm is not suitable for directed graphs.
    This algorithm ignores self-loops and parallel edges as
    clique is not conventionally defined with such edges.
    There are often many cliques in graphs.  This algorithm can
    run out of memory for large graphs.
    References
    ----------
    .. [1] Bron, C. and Kerbosch, J. 1973.
       Algorithm 457: finding all cliques of an undirected graph.
       Commun. ACM 16, 9 (Sep. 1973), 575-577.
       http://portal.acm.org/citation.cfm?doid=362342.362367
    .. [2] Etsuji Tomita, Akira Tanaka, Haruhisa Takahashi,
       The worst-case time complexity for generating all maximal
       cliques and computational experiments,
       Theoretical Computer Science, Volume 363, Issue 1,
       Computing and Combinatorics,
       10th Annual International Conference on
       Computing and Combinatorics (COCOON 2004), 25 October 2006, Pages 28-42
       http://dx.doi.org/10.1016/j.tcs.2006.06.015
    .. [3] F. Cazals, C. Karande,
       A note on the problem of reporting maximal cliques,
       Theoretical Computer Science,
       Volume 407, Issues 1-3, 6 November 2008, Pages 564-568,
       http://dx.doi.org/10.1016/j.tcs.2008.05.010
    """
    if len(data_graph) == 0:
        return
    cdef Py_ssize_t i, j, u
    cdef list Q = [None]

    cdef Py_ssize_t n = len(data_graph)

    cdef set subg_q, adj_q, cand_q
    cdef set subg = set(range(n))
    cdef set cand = set(range(n))

    cdef set ext_u
    cdef list stack = []

    # construct adjacency
    cdef dict adj_graph = {x: set() for x in range(n)}
    for i in range(n):
        for j in range(n):
            if j != i and data_graph[i,j] == 1:
                adj_graph[i].add(j)

    u = max(
        subg,
        key=lambda x: len(cand & adj_graph[x]),
    )
    ext_u = cand - adj_graph[u]

    try:
        while True:
            if ext_u:
                q = ext_u.pop()
                cand.remove(q)
                Q[-1] = q
                adj_q = adj_graph[q]
                subg_q = subg & adj_q
                if not subg_q:
                    yield Q[:]
                else:
                    cand_q = cand & adj_q
                    if cand_q:
                        stack.append((subg, cand, ext_u))
                        Q.append(None)
                        subg = subg_q
                        cand = cand_q
                        u = max(
                            subg,
                            key=lambda u: len(cand & adj_graph[u]),
                        )
                        ext_u = cand - adj_graph[u]
            else:
                Q.pop()
                subg, cand, ext_u = stack.pop()
    except IndexError:
        pass
# end def

def find_cliques(
        data_graph,
        num_seeds_needed: int,
        cutoff: int,
) -> List[List[int]]:
    """Compute all maximal cliques of
    Keep the cutoff something reasonable and this search is very fast.  Over
    50 nodes starts to get really slow if no cutoff exists

    Args:
        data_graph (numpy 2D ndarray): Undirected graph representing the
            connected system
        num_seeds_needed (int): Minimum size of a passing clique
        cutoff (int): Number of permissible solutions to allow

    Returns:
        List[List[int]]: a list of lists of indices no longer than length cutoff
            but possibly shorter depending on the number of maximal cliques in
            the graph
    """
    cdef list clique
    cdef list out_list = []
    cdef int solution_count = 0

    for clique in _find_cliques(data_graph):
        if len(clique) > num_seeds_needed - 1:
            out_list.append(clique)
            solution_count += 1
            if solution_count > cutoff:
                break
    return out_list


def hammingGraphEqual(list sequences) -> np.ndarray:
    """Compute hamming distance between all strings in a list

    Args:
        sequences (list): list of python strings,
            (bytes or unicode, all of the same length

    Returns:
        ndarray: A numpy 2D array of the (hamming distances, index) between all
            sequences in the list.  The index is included to allow for sorting
    """
    cdef char* seq1_cstr
    cdef char* seq2_cstr
    cdef Py_ssize_t len1, len2

    cdef int* hm_arr

    cdef Py_ssize_t len_sequences = len(sequences)
    cdef Py_ssize_t i, j, row_idx
    cdef int res

    # The index field is to aid in sorting each row by hamming distance
    cdef cnp.ndarray[hamming_node_t, ndim=2] hm_np
    hm_np = np.zeros(
        (len_sequences, len_sequences),
        dtype=np.dtype(
            [('distance', np.int32), ('index', np.int32)],
            align=True,
        )
    )

    # Since hm_np is aligned, we can write directlty to each field
    hm_arr = <int*> cnp.PyArray_DATA(hm_np)

    for i in range(len_sequences - 1):
        seq1_cstr = c_util.obj_to_cstr_len(
            sequences[i],
            &len1,
        )
        row_idx = i*len_sequences*2
        for j in range(i + 1, len_sequences):
            #print(i, j)
            seq2_cstr = c_util.obj_to_cstr_len(
                sequences[j],
                &len2,
            )
            res = ss_hamming(
                <const char *> seq1_cstr,
                <const char *> seq2_cstr,
                <int> len1,
            )
            hm_arr[row_idx + 2*j] = res
            # mirror result
            hm_arr[(j * len_sequences + i) * 2] = res

    # write index fields
    for i in range(len_sequences):
        row_idx = i*len_sequences*2 + 1 # add the offset here
        for j in range(len_sequences):
            hm_arr[row_idx + 2 * j] = j

    return hm_np


def hammingGraph(list sequences) -> np.ndarray:
    """Compute hamming distance between all strings in a list. Strings
    can be different lengths but must be sorted longest to shortest

    Args:
        sequences (list): list of python strings,
            (bytes or unicode, all of the same length

    Returns:
        A numpy 2D array of the (hamming distances, index) between all sequences
            in the list.  The index is included to allow for sorting
    """
    cdef char* seq1_cstr    # the shorter sequence
    cdef char* seq2_cstr    # the longer sequence
    cdef Py_ssize_t len1, len2

    cdef int* hm_arr

    cdef Py_ssize_t len_sequences = len(sequences)
    cdef Py_ssize_t i, j, row_idx
    cdef int res

    # the index field is to aid in sorting each row by hamming distance
    cdef cnp.ndarray[hamming_node_t, ndim=2] hm_np
    hm_np = np.zeros(
        (len_sequences, len_sequences),
        dtype=np.dtype(
            [('distance', np.int32), ('index', np.int32)],
            align=True,
        )
    )

    # since hm_np is aligned, we can write directlty to each field
    hm_arr = <int*> cnp.PyArray_DATA(hm_np)

    for i in range(len_sequences - 1):
        seq1_cstr = c_util.obj_to_cstr_len(
            sequences[i],
            &len1,
        )
        row_idx = i * len_sequences*2
        for j in range(i + 1, len_sequences):
            #print(i, j)
            seq2_cstr = c_util.obj_to_cstr_len(
                sequences[j],
                &len2,
            )
            res = ss_minHamming(
                <const char *> seq2_cstr,
                <const char *> seq1_cstr,
                <int> len2,
                <int> len1,
            )
            hm_arr[row_idx + 2*j] = res
            # mirror result
            hm_arr[(j*len_sequences + i)*2] = res
        # end for
    # end for

    # write index fields
    for i in range(len_sequences):
        row_idx = i*len_sequences*2 + 1 # add the offset here
        for j in range(len_sequences):
            hm_arr[row_idx+2*j] = j

    return hm_np


def hamming2XGraph(
        list sequences,
        int min_distance=2,
) -> np.ndarray:
    """Compute the graph of having a hamming distance at least `min_distance`
    between every sequence in both the forward and reverse direction

    Args:
        sequences (list): list of python strings, (bytes or unicode, all
            of the same length.
        min_distance (int): minimum hamming distance between sequences

    Returns:
        ndarray: A 2D numpy array of np.uint8 that represent the compatibility
            of any two elements in the list
    """
    cdef char* seq1_cstr
    cdef char* seq2_cstr
    cdef Py_ssize_t len1, len2

    cdef unsigned char* check_arr

    cdef Py_ssize_t len_sequences = len(sequences)
    cdef Py_ssize_t i, j
    cdef Py_ssize_t row_idx = 0
    cdef unsigned char res

    cdef int solution_count = 0

    # the index field is to aid in sorting each row by hamming distance
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] check_np
    check_np = np.zeros(
        (len_sequences, len_sequences),
        dtype=np.uint8,
    )


    # since check_np is aligned, we can write directlty to each field
    check_arr = <unsigned char*> cnp.PyArray_DATA(check_np)

    # create the check matrix
    for i in range(len_sequences - 1):
        seq1_cstr = c_util.obj_to_cstr_len(
            sequences[i],
            &len1,
        )

        row_idx = i*len_sequences
        for j in range(i + 1, len_sequences):
            #print(i, j)
            seq2_cstr = c_util.obj_to_cstr_len(
                sequences[j],
                &len2,
            )

            res = ss_hamming2XCheck(
                <const char *> seq1_cstr,
                <const char *> seq2_cstr,
                <int> len1,
                min_distance,
            )
            check_arr[row_idx + j] = res
            # mirror result
            check_arr[j*len_sequences + i] = res

    return check_np


def hammingAnd3pMismatchGraph(
        list sequences,
        int min_distance = 2,
) -> np.ndarray:
    """Compute the graph of having a hamming distance
    greater than min_distance between elements and having at least one mismatch
    in the last 3 basepairs at the 3 prime end

    Args:
        sequences (list): list of python strings,
            (bytes or unicode, all of the same length.
        min_distance (int): minimum hamming distance between sequences

    Returns:
        A 2D numpy array of np.uint8  that represent the compatibility of
        any two elements in the list
    """
    cdef char* seq1_cstr
    cdef char* seq2_cstr
    cdef Py_ssize_t len1, len2

    cdef unsigned char* check_arr

    cdef Py_ssize_t len_sequences = len(sequences)
    cdef Py_ssize_t i, j
    cdef Py_ssize_t row_idx = 0
    cdef unsigned char res

    cdef int solution_count = 0

    # the index field is to aid in sorting each row by hamming distance
    cdef cnp.ndarray[cnp.uint8_t, ndim=2] check_np
    check_np = np.zeros(
        (len_sequences, len_sequences),
        dtype=np.uint8,
    )


    # since check_np is aligned, we can write directlty to each field
    check_arr = <unsigned char*> cnp.PyArray_DATA(check_np)

    # create the check matrix
    for i in range(len_sequences - 1):
        seq1_cstr = c_util.obj_to_cstr_len(
            sequences[i],
            &len1,
        )

        row_idx = i*len_sequences
        for j in range(i + 1, len_sequences):
            #print(i, j)
            seq2_cstr = c_util.obj_to_cstr_len(
                sequences[j],
                &len2,
            )

            res = ss_hammingAnd3pMismatchCheck(
                <const char *> seq1_cstr,
                <const char *> seq2_cstr,
                <int> len1,
                min_distance,
            )
            check_arr[row_idx + j] = res
            # mirror result
            check_arr[j*len_sequences + i] = res

    return check_np
