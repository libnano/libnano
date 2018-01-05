"""
 SLS
implemented from

DC Tulpan and HH Hoos, Hybrid Randomised Neighbourhoods Improve
Stochastic Local Search for DNA Code Design. 2003

"""
from os.path import abspath, dirname, join
import sys
LOCAL_DIR = abspath(dirname(__file__))
# For package imports
sys.path.insert(0, abspath(join(LOCAL_DIR, '..')))

from libnano import seqstr as ss
from libnano import seqint as si
import copy
import itertools as it
import random
from collections import deque

# ATGATTC

# TTGATTC
# AAGATTC
# ATCATTC
# ATGTTTC
# ATGAATC
# ATGATAC
# ATGATTG


def generate50GC8mers():
    out = []
    for s in it.product('ACGT', repeat=8):
        if s.count('G') + s.count('C') == 4:
            out.append(''.join(s))
    # end for
    return out
# end def

def generateInitialSet(choices, k):
    return random.sample(choices, k)
# end def

def reverseComp(seq):
    out1 = seq[-1::-1]
    out2 = ''
    for c in out1:
        if c == 'A':
            out2 += 'T'
        elif c == 'T':
            out2 += 'A'
        elif c == 'G':
            out2 += 'C'
        else:
            out2 += 'G'
    return out2
# end def

def hamming(seq1, seq2):
    dist = 0
    for i in range(len(seq1)):
        dist += 1
    return dist
# end def


# check_constraints
def checkConstraints(candidate, seq_set, hd=5):
    # check GC content
    # if candidate.count('G') + candidate.count('C') != 4:
    #     return False
    # check hamming distance
    for seq in seq_set:
        if ss.hammingDistance(candidate, seq) < hd:
            return False
        if ss.hammingDistance(candidate, ss.reverseComplement(seq)) < hd:
            return False
    return True
# end def

def checkConstraintFailedCount(candidate, idxs, seq_set):
    failed = 0
    for idx in idxs:
        seq = seq_set[idx]
        if ss.hammingDistance(candidate, seq) < 5:
        # if hamming(candidate, seq) < 5:
            failed += 1
        if ss.hammingDistance(candidate, ss.reverseComplement(seq)) < 5:
        # if hamming(candidate, reverseComp(seq)) < 5:
            failed += 1
    return failed

# # check_constraints
# def checkConstraintFailedCount(candidate, seq_set):
#     # check GC content
#     failed = [0, 0, 0]
#     # if candidate.count('G') + candidate.count('C') != 4:
#     #     satisfied[0] = 1
#     # check hamming distance
#     for seq in seq_set:
#         if ss.hammingDistance(candidate, seq) < 5:
#             failed[1] = 1
#         if ss.hammingDistance(candidate, ss.reverseComplement(seq)) < 5:
#             failed[2] = 1
#     return sum(failed)
# # end def

def countAllViolations(seq_set):
    return sum([checkConstraintFailedCount(s, seq_set) for s in seq_set])

def neighborhood(seq):
    out = [None]*8
    for i in range(8):
        c_in = seq[i]
        if c_in == 'A':
            c = 'T'
        elif c_in == 'T':
            c = 'A'
        elif c_in == 'G':
            c = 'C'
        else:
            c = 'G'
        if i < 7:
            out[i] = seq[:i] + c + seq[i+1:]
        else:
            out[7] = seq[:7] + c
    return out
# end def

def runSLS(max_tries=5, max_steps=120000):
    s_all = generate50GC8mers()
    set_max_idx = 6
    fail_idx_list = [0, 0]
    theta_list = [0]*8+[1]*2
    for i in range(max_tries):
        print("try", i, '...')
        seq_set = generateInitialSet(s_all, 7)
        set_last = copy.copy(seq_set)
        last_count = countAllViolations(set_last)
        for j in range(max_steps):
            fail = False
            for seq in seq_set:
                if not checkConstraints(seq, seq_set):
                    fail = True
                    break
            if fail == False:
                return seq_set, True
            # fail_words_count = 0
            # fail_words_last_idx = -1
            fail_idx_list = [-1, -1]
            for q in range(2):
                rand_idx = random.randint(0, set_max_idx)
                if not checkConstraints(seq_set[rand_idx], seq_set):
                    fail_idx_list[q] = rand_idx
            # now generate M
            if fail_idx_list[0] == -1:
                N1 = set()
            else:
                N1 = set(neighborhood(seq_set[fail_idx_list[0]]))
            if fail_idx_list[1] == -1:
                N2 = set()
            else:
                N2 = set(neighborhood(seq_set[fail_idx_list[1]]))
            M = list(N1.union(N2))
            if random.choice(theta_list):
                w_prime = random.choice(M)
            else:
                min_count = 8
                min_idx = 0
                for z, candidate in enumerate(M):
                    num_failed = checkConstraintFailedCount(candidate, seq_set)
                    if num_failed == 0:
                        min_idx = z
                        break
                    elif num_failed < min_count:
                        min_count = num_failed
                        min_idx = z
                w_prime = M[min_idx]
                if w_prime in N1:
                    seq_set[fail_idx_list[0]] = w_prime
                elif fail_idx_list[1] != -1:
                    seq_set[fail_idx_list[1]] = w_prime
                new_count = countAllViolations(seq_set)
                if last_count > new_count:
                    set_last = copy.copy(seq_set)
                    last_count = new_count
        # end for
    # end for
    return set_last, False, last_count
# end def

def writeOut(items, fn):
    with open(fn, 'w') as fd:
        fd.write("%d\n\n" % (len(items)))
        for item in items:
            fd.write("%d\n%s\n%s\n\n" % (item[0],
                                        str(item[1]),
                                        str(item[2])) )
# end def

def brute(start_index, pool, k=27, keep_threshold=26, index_search_size=None, write_file=True):
    seq_set = []
    seq_set_idxs = []
    max_set = []
    max_set_idxs = []
    passing_sets = []
    pool_len = len(pool)
    if index_search_size == None:
        index_search_size = pool_len
    elif index_search_size > pool_len:
        index_search_size = pool_len

    set_size_per_idx = [0]*pool_len
    idx0 = start_index
    lim = pool_len - 1
    print("starting at index %d, for %d" % (idx0, index_search_size))
    for i in range(index_search_size):
        # clear the set
        seq_set = []
        seq_set_idxs = []
        for j in range(pool_len):
            j_idx = j + idx0
            if j_idx > lim:
                j_idx -= pool_len
            candidate = pool[j_idx]
            if checkConstraints(candidate, seq_set):
                seq_set.append(candidate)
                seq_set_idxs.append(j_idx)
            if len(seq_set) > k:
                break
        # end for j
        len_ss = len(seq_set)
        set_size_per_idx[idx0] = len_ss
        # if len(max_set) < len_ss:
        #     print("new max found %d starting at index %d" % (len_ss, idx0))
        #     max_set = seq_set
        #     max_set_idxs = seq_set_idxs
        #     # if len_ss > k:
        #     #     break
        if len_ss > keep_threshold:
            print("found %d at index %d" % (len_ss, idx0))
            passing_sets.append((len_ss, seq_set_idxs, seq_set))
        idx0 += 1
        if idx0 > lim:
            idx0 = 0
    # end for i
    print("ended at index %d" % (idx0))
    if write_file:
        writeOut(passing_sets, "hamming5merGC50.txt")
    # return max_set_idxs, max_set
    return passing_sets
# end def

def sharpen(passing_sets, pool):
    pool_idxs = list(range(len(pool)))
    pool_idxs_set = set(pool_idxs)
    for len_ss, seq_set_idxs, seq_set in passing_sets:
        temp = set(seq_set_idxs)
        not_in_seq_set_idxs = list(pool_idxs_set.difference(temp))

        # check fail counts
        failed_counts = [0]*len_ss
        for i in range(len_ss):
            candidate = seq_set[i]
            failed_counts[i] = checkConstraintFailedCount(candidate, not_in_seq_set_idxs, pool)
        print(failed_counts)
        # # try sharpening
        # max_fail_counts = max(fail_counts)
        # max_fail_idxs = [i for i, j in enumerate(fail_counts) if j == m]
        # for idx in max_fail_counts:
        #     new_set = copy.copy(seq_set)
        #     new_set.pop(idx)

def edits1(word):
    alphabet = 'ACGT'
    splits     = [(word[:i], word[i:]) for i in range(len(word) + 1)]
    deletes    = [a + b[1:] for a, b in splits if b]
    replaces   = [a + c + b[1:] for a, b in splits for c in alphabet if b]
    return set(deletes + replaces)
# end def

def edits2(word):
    return set(e2 for e1 in edits1(word) for e2 in edits1(e1))

def correct(word, seq_set):
    """ correct up to two errors in 8 mer
    """
    candidates = set([word]).union(edits1(word)).union(edits2(word))
    return max(candidates, key=seq_set.count)

if __name__ == '__main__':
    import json
    # s_all = generate50GC8mers()
    # print(len(s_all))
    # print(brute(0, s_all, 27))

    # psets = brute(0, s_all, write_file=False)
    # with open('passed.json', 'w') as fd:
    #     json.dump({'passed': psets}, fd)

    with open('passed.json') as fd:
        psets = json.load(fd)['passed']
    print(len(psets))

    # print(psets[0])
    # sharpen(psets, s_all)

    # error correction test
    ss = psets[0][2]
    word = ss[0]
    # new_word = 'T' + word[1:]
    new_word = 'T' + word[1] + 'G' + word[3:]
    print("correct: ", word, "variant:", new_word)
    print("found:   ", correct(new_word, ss))

