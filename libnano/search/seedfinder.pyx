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
libnano.search.seedfinder
~~~~~~~~~~~~~~~~~~~~~~~~~

Generate seed strings for use in homology searches

refs:
Darling, Aaron E., Bob Mau, and Nicole T. Perna. "progressiveMauve: multiple
genome alignment with gene gain, loss and rearrangement." PloS one 5.6 (2010):
e11147.

Xu, Jinbo, et al. "Optimizing multiple spaced seeds for homology search."
Journal of Computational Biology 13.7 (2006): 1355-1368.

Kucherov G, Noe L, Roytberg M (2005) Multiseed lossless filtration.
IEEE/ACM Trans Comput Biol Bioinformatics 2: 51–61.
'''
from typing import (
    Iterator,
    List,
    Tuple,
    Union,
)


def checkSeed(
        seed: str,
        m: int,
        k_min: int = 1,
        k_max: int = 3,
) -> int:
    '''Check that a seed is valid

    Args:
        seed (str): candidate seed string composed of `#` and `-` characters
        m (int): length of the mer to evaluate
        k_min (Optional[int]): minimum mismatch count to satisfy
        k_max (Optional[int]): maximum mismatch count to satisfy

    Returns:
        int: -1 if no k mismatch found, or a value in the range `k_min` to
            `k_max` representing the most mismatches the seed satisfies
    '''
    seed_length = len(seed)
    # get indices of non-jokers
    seed_idxs = [idx for idx, pos in enumerate(seed) if pos == "#"]
    seed_weight = len(seed_idxs)
    num_tests = m - seed_length + 1
    best_k = -1
    # print("seed idxs", seed_idxs)
    for k in range(k_min, k_max+1):
        # print("number of perms for k of %d: %d" % (k, len(perms)))
        for perm in combo(m, k):
            did_mismatches = seedhashseq(seed_idxs, perm, num_tests)
            if sum(did_mismatches) == num_tests: # no one passed
                return best_k
        best_k = k
    return best_k


def checkSeedAtK(
        seed: str,
        m: int,
        perm_list: List[int],
        perm_limit: int,
) -> Tuple[bool, Union[List[str], None]]:
    '''Given a list of permutations for a given mismatch k (calculated by
    caller) determine if a seed satisfies finding those mismatches

    Args:
        seed (str): candidate seed string composed of `#` and `-` characters
        m (int): length of the mer to evaluate
        perm_list (List[int]): A list of permutations of k mismatches in a m-mer
        perm_limit (int): the number of allowable failed permutations before the
            evaluation of this seed breaks out of the loop

    Returns:
        Tuple[bool, List]: True if no failure in the given permutations
            otherwise False
            list of failed_perms in string form or None on pass

    '''
    seed_length = len(seed)
    # get indices of non-jokers
    seed_idxs = [
        idx for idx, pos in enumerate(seed) if pos == "#"
    ]
    seed_weight = len(seed_idxs)
    num_tests = m - seed_length + 1
    # print("number of perms for k of %d: %d" % (k, len(perms)))
    perm_failure = False
    failed_perms = []
    for perm in perm_list:
        did_mismatches = seedhashseq(
            seed_idxs,
            perm,
            num_tests,
        )
        if sum(did_mismatches) == num_tests: # no one passed
            perm_failure = True
            failed_perms.append(str(perm))
            if len(failed_perms) > perm_limit:
                break
    if perm_failure:
        if len(failed_perms) < perm_limit:
            return False, failed_perms
        else:
            return False, None
    else:
        return True, None


def checkSeeds(
        seeds: List[str],
        m: int,
        k_min: int = 1,
        k_max: int = 3,
) -> int:
    '''Check a list of seeds for their combined mismatch screening power

    Args:
        seeds (List[str]): candidate seed strings composed of `#` and `-`
            characters
        m (int): length of the mer to evaluate
        k_min (Optional[int]): minimum mismatch count to satisfy
        k_max (Optional[int]): maximum mismatch count to satisfy

    Returns:
        int: -1 if no k mismatch found, or a value in the range `k_min` to
            `k_max` representing the most mismatches the list of seeds satisfy
    '''
    # num_seeds = len(seeds)
    # get indices of non-jokers
    seed_idxs_list = []
    num_tests_list = []
    for seed in seeds:
        seed_length = len(seed)
        seed_idxs_list.append(
            [idx for idx, pos in enumerate(seed) if pos == "#"]
        )
        num_tests_list.append(
            m - seed_length + 1
        )
    # total_tests = sum(num_tests_list)
    best_k = -1
    # print("seed idxs", seed_idxs)
    for k in range(k_min, k_max+1):
        # print("number of perms for k of %d: %d" % (k, len(perms)))
        for perm in combo(m, k):
            did_mismatches = []
            at_least_one_catches = False
            for seed_idxs, num_tests in zip(seed_idxs_list, num_tests_list):
                #did_mismatches += seedhashseq(seed_idxs, perm, num_tests)
                res = seedhashseq(seed_idxs, perm, num_tests)
                did_mismatches.append((res, perm))
                if sum(res) < num_tests:
                    at_least_one_catches = True

            if not at_least_one_catches:
                return best_k
            #if sum(did_mismatches) == total_tests: # no one passed
            #    return best_k
        best_k = k
    return best_k


cdef list seedhashseq(
        seed_idxs: List[int],
        seq: List[int],
        int num_tests,
):
    '''Use a seed to hash a sequence object since the length of the seed
    will be shorter than the length of the seq

    Args:
        seed_idxs [List(int)]: List of indices of a seed string that equal "#"
        seq (List[int]): the list of 1's and 0's representing the sequence with
            mismatches where a 1 represents a correct base and a 0 represents a
            mismatch
        num_tests (int): the number of tests a sequence can be hashed by a seed
            computed by caller.  Equal to `len(seq) - len(seed) + 1`

    Returns:
        List[int]: Whether or not a mismatch is caught (1) or not (0) for each
            test in `num_tests`
    '''
    did_mismatches = [0 for i in range(num_tests)]
    for i in range(num_tests):
        for j in seed_idxs:
            if seq[i + j] != 1:   # a 1 is a correct base
                # caught a mismatch
                did_mismatches[i] = 1
                break

    return did_mismatches


def combo(
        m: int,
        k: int,
) -> Iterator[List[int]]:
    '''Iterator of all lists of length m
    with k substitions/mismatches represented with a 0
    and correct bases represented with a 1

    Args:
        m (int): length of the test mer to create
        k (int): mismatch count

    Returns:
        Iterator[List[int]]: Iterator of 1, 0 represention of a correct bases
            (1's) and incorrect bases (0's)
    '''
    k0 = k
    joker = [0]*k

    def nest(m, kt, q, start):
        # print("kt is", kt)
        for i in range(start, m-k+1+q):
            joker[k - kt] = i
            if kt > 1:
                for val in nest(m, kt - 1, q + 1, i + 1):
                    yield val
            else:
                out = [1]*m    # 1 means a desired base
                for z in joker:
                    out[z] = 0 # mismatch is a zero
                yield out
    # end def
    return nest(m, k0, 0, 0)
# end def

def seed_combo(n: int, weight: int) -> Iterator[str]:
    '''Given a seed of length n get all of a certain weight
    uses combo to do the comination despite creating something different
    than other calls to combo in this mode

    Args:
        n (int): length of the sedd to create
        weight (int): weight of the seed (number of `#`'s)

    Returns:
        Iterator[str]: Iterator Seed strings with `#` representing
            must matches and `-` representing don't care /skip for the hashing
    '''
    for val in combo(n, weight):
        yield "".join("-" if x == 1 else '#' for x in val)

'''
from pg 12 of

Kucherov G, Noe L, Roytberg M (2005) Multiseed lossless filtration.
IEEE/ACM Trans Comput Biol Bioinformatics 2: 51–61.
ref: http://hal.inria.fr/docs/00/35/48/10/PDF/KucherovNoeRoytberg.pdf

Example 1:
50,5 is solved by
#-#-#---#-----#-#-#---#-----#-#-#---# and ###-#--###-#--###-#
the first seed is an expansion of the second which is a solution
to the 25,2 problem.  This is because:

Lemma 2: If a family Fi = i ⊗ F solves an (im, k)-problem,
then F solves both the (im, k)-problem and the (m, ⌊k/i⌋)-problem

Lemma 1: If a family F solves an (m, k)-problem,
then the (im,(i + 1)k − 1)-problem is solved both by family
F and by its i-regular expansion Fi = i ⊗ F.

'''
def findSeed(
        m: int,
        k: int,
        wmin: float = 0.4,
        wmax: float = 0.8,
        nmin: float = 0.74,
        nmax: float = 0.76,
) -> List[Tuple[int, int, str, int, int]]:
    '''Find a seed that satisfies m, k givven weight restrictions and length of
    seed restrictions

    Args:
        m (int): length of the test mer to create
        k (int): mismatch count
        nmin, nmax (float): multipliers of `len(m)` to consider seed lengths
        wmin, wmax (float): multipliers of `len(n)` to consider seed weights
            (number of `#`'s)

    Returns:
        List[Tuple[int, int, str, int, int]]: List of tuples of::

            (m, k, seed, n, w)

            where `n` is the seed length and `w` is seed weight
    '''
    n_min = int(nmin*m)   # seed length
    n_max = int(nmax*m)
    seed_list = []
    for n in range(n_min, n_max + 1):
        w_min = int(wmin*n)
        w_max = int(wmax*n)
        # print(n, w_min, w_max)
        for w in range(w_min, w_max):
            for seed in seed_combo(n, w):
                best_k = checkSeed(seed, m, k_min=k, k_max=k+1)
                if best_k == k:
                    seed_list.append((m, k, seed, n, w))
                    break # onto the next weight
    return seed_list


SP_T = Tuple[str, str, int, float, float, int, float, float, int]


def findSeedPairs(
        m: int,
        k: int,
        w_min_frac: float = 0.6,
        w_max_frac: float = 0.8,
        w_exact: Tuple[float, float] = None,
        n_min_frac: float = 0.7,
        n_max_frac: float = 0.9,
        n_exact: Tuple[float, float] = None,
        is_greedy: bool = True,
        debug: bool = False,
) -> SP_T:
    '''Find a pair of seeds A and B for the `m, k` problem
    keep the highest of weight only.  Greedy so continues once it gets
    one answer for a given weight.  This needs a little cleaning up.

    Args:
        m (int): length of the test mer to create
        k (int): mismatch count

    Returns:
        List[Tuple]: List of tuples of::

                (seed_A, seed_B, best_k, n_A, w_A, length_A, n_B, w_B, length_B)

            where `n` is a seed length and `w` a seed weight, and `length`
            is a length of a seeds indidual failure permutations at solving the
            `m, k` problem
    '''
    if n_exact is None:
        n_min = int(n_min_frac*m)   # seed length
        n_max = int(n_max_frac*m)
    else:
        n_min, n_max = n_exact

    if w_exact is None:
        w_min = int(w_min_frac*n_min)
        w_max = int(w_max_frac*n_max)
    else:
        w_min, w_max = w_exact

    seed_list = []
    perm_list = [p for p in combo(m, k)]

    '''
    set the allowable failures of the permutation of the seed along sequence m
    since we can use more the one seed to find all hits.
    '''
    temp = len(perm_list) // 2
    perm_limit = temp if temp < 50 else 50
    if debug:
        print(f'Perm Limit: {perm_limit}')

    got_seed = False
    for w in range(w_min, w_max):
        if debug:
            print(f'Now searching weight: {w}')
        # look for seed families over constant weight
        failures_list = []
        for n in range(n_min, n_max + 1):
            if w > w_max_frac*n:
                continue
            for seed in seed_combo(n, w):
                did_pass, failed_perms = checkSeedAtK(
                    seed,
                    m,
                    perm_list,
                    perm_limit,
                )
                if did_pass and n > 0.6*m:
                    # see if the seed solves k+1
                    if checkSeed(
                        seed,
                        m,
                        k_min=k+1,
                        k_max=k+1,
                    ) == -1:
                        seed_list.append((seed, n, w))
                        if is_greedy:
                            break # onto the next weight
                elif failed_perms is not None:  # failures
                    failures_list.append(
                        (seed, set(failed_perms), n, w),
                    )
        # Now search through list for non-overlapping perm failures.
        failures_list_len = len(failures_list)
        if debug:
            print(
                f'Finished k pass with {failures_list_len} to check'
            )
        if failures_list_len > 1:
            for i in range(failures_list_len - 1):
                if got_seed:    # got a seed for this weight
                    got_seed = False
                    break
                a = failures_list[i]
                for j in range(i+1, failures_list_len):
                    b = failures_list[j]
                    # if the two seeds solve it, then the intersection
                    # should be empty
                    if len(a[1].intersection(b[1])) == 0:
                        # see if the seeds solve k+1
                        best_k = checkSeeds(
                            [a[0], b[0]],
                            m,
                            k_min=k+1,
                            k_max=k+1,
                        )
                        # best_k should be -1
                        if best_k != k:
                            s_a, fail_set_a, n_a, w_a = a
                            l_a = len(fail_set_a)

                            s_b, fail_set_b, n_b, w_b = b
                            l_b = len(fail_set_b)

                            sol = (
                                s_a,
                                s_b,
                                best_k,
                                n_a,
                                w_a,
                                l_a,
                                n_b,
                                w_b,
                                l_b,
                            )
                            seed_list.append(sol)
                            if is_greedy:
                                got_seed = True
                                break
                        else:
                            print(
                                f'Issue: {m}, {k}, {a[0]}, {b[0]}',
                            )
    return seed_list


# def test25_2():
#     perm_list = [p for p in combo(25, 2)]
#     best_k1, failed_perms1 = checkSeedAtK("#####-##---#####-##", 25, 2, perm_list)
#     best_k2, failed_perms2 = checkSeedAtK("#-##---#####-##---####", 25, 2, perm_list)
#     if failed_perms1 is None:
#         res1 = best_k1, "k"
#     else:
#         res1 = len(failed_perms1)
#     if failed_perms2 is None:
#         res2  = best_k2, "k"
#     else:
#         res2 = len(failed_perms2)
#     print("The lengths of failures", res1, res2)
#


# if __name__ == '__main__':

#     def test_checkSeed(seed, m):
#         res = checkSeed(seed, m)
#         print("number of mismatches allowed: %d" %(res))
#         return res

#     assert(test_checkSeed("#-##--#-##", 15) == 2)
#     assert(test_checkSeed("#-##-##-##", 15) == 1)
#     assert(test_checkSeed("#######-##", 15) != 2)
#     # Example 1         #-#-#---#--###-#--###-#
#     assert(checkSeeds(["#-#-#---#-----#-#-#---#-----#-#-#---#",
#         "###-#--###-#--###-#"], 50, k_min=4, k_max=6) == 5)
# '###-###-##--'
# '#-#-#---#-#-#---#-#----'
# '#---#---#-------#---#---#-------#---#---------'

# '''
# AANNN
# ANANN
# ANNAN
# ANNNA
# NAANN
# NANAN
# NANNA
# NNAAN
# NNANA
# NNNAA

# 30,2
# [('####-#---###----------', 8),
#  ('####-#---####---------', 9),
#  ('####-#---####-#-------', 10),
#  ('####-#---####-#---#---', 11),
#  ('####-#---####-#---##--', 12),
#  ('####-#---####-#---###-', 13),
#  ('####-#---####-#---####', 14)]
# '''
