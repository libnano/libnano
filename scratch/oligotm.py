from __future__ import print_function

import sys
import random
import math


# ~~~~~~~~~~~~~~~~~~~~~~~~ DNA Sequence Manipulation ~~~~~~~~~~~~~~~~~~~~~~~~ #

_DNAcomp = str.maketrans('ACGTacgt','TGCATGCA')

def reverseComplement(seq):
    return seq.translate(_DNAcomp)[::-1]
# end def

rc = reverseComplement

def complement(seq):
    return seq.translate(_DNAcomp)
# end def

def randomSeq(length):
    return ''.join([random.choice('ATGC') for x in range(length)])
# end def

# ~~~~~~~~~~~~~~~~~~~~ GENERAL CONSTANTS AND CALCULATIONS ~~~~~~~~~~~~~~~~~~~ #

R = 1.9872e-3  # Gas constant (1 / kcal * mol)
KELVIN = 273.15


def cToK(deg_c):
    ''' Convert degrees Celsius to degrees Kelvin '''
    return deg_c + KELVIN


def kToC(deg_k):
    ''' Convert degrees Kelvin to degrees Celsius '''
    return deg_k - KELVIN


def calcKa(dg, deg_c):
    ''' Return the association constant at a given temperature '''
    return math.e**(-dg/(R * cToK(deg_c)))


def calcDg(ds, dh, deg_c):
    ''' Return the dg at a given temp using the provided ds and dh '''
    return dh - ds * cToK(deg_c)


def calcRandCoil(dg, deg_c):
    ''' Return the percent of randomly coiled oligo with dg at deg_c degrees '''
    return 1/(calcKa(dg, deg_c) + 1)


# ~~~~~~~~~~~~~~~~~ PYTHON IMPLEMENTATION OF PRIMER3 OLIGOTM ~~~~~~~~~~~~~~~~ #

def divalentToMonovalent(divalent, dntp):
    if divalent == 0:
        dntp = 0
    if divalent < dntp:
        divalent = dntp
    print("d2m", divalent-dntp)
    return 120 * math.sqrt(divalent-dntp)

def calcThermo(seq, conc_nm=50, monovalent=50, divalent=0, dntp=0.8):
    ''' Return the thermo parameters for DNA under specified salt cond.

    '''

    enthalpies = {
        'AA': 79, 'AT': 72, 'AG': 78, 'AC': 84,
        'TA': 72, 'TT': 79, 'TG': 85, 'TC': 82,
        'GA': 82, 'GT': 84, 'GG': 80, 'GC': 98,
        'CA': 85, 'CT': 78, 'CG': 106, 'CC': 80
    }
    entropies = {
        'AA': 222, 'AT': 204, 'AG': 210, 'AC': 224,
        'TA': 213, 'TT': 222, 'TG': 227, 'TC': 222,
        'GA': 222, 'GT': 224, 'GG': 199, 'GC': 244,
        'CA': 227, 'CT': 210, 'CG': 272, 'CC': 199
    }
    dH = dS = 0
    # Calculate oligo symmetry
    sym = seq == reverseComplement(seq)
    # Calculate NN uncorrected dS and dH for oligo
    for idx in range(len(seq)-1):
        dH += enthalpies[seq[idx:idx+2]]
        dS += entropies[seq[idx:idx+2]]
    # Terminal AT penalty and initiation parameters (combined)
    if seq[0] in 'AT':
        dH += -23
        dS += -41
    else:
        dH += -1
        dS += 28
    if seq[-1] in 'AT':
        dH += -23
        dS += -41
    else:
        dH += -1
        dS += 28
    if sym:
        dS += 14
    dH *= -100.0
    dS *= -0.1
    # Convert divalent salt and dntp conc. to monovalent equivalencies
    monovalent += divalentToMonovalent(divalent, dntp)
    dS = dS + 0.368 * (len(seq) - 1) * math.log(monovalent / 1000.0)
    # Account for oligo symmetry and calculate tm
    if sym:
        tm = dH / (dS + 1.987 * math.log(conc_nm/1.0e9)) - KELVIN
    else:
        print("conc", conc_nm, dS)
        tm = dH / (dS + 1.987 * math.log(conc_nm/4.0e9)) - KELVIN
    return dH, dS, tm

def calcTm(seq, conc_nm=50, monovalent=50, divalent=0.0, dntp=0.8):
    _, _, tm = calcThermo(seq, conc_nm=conc_nm, monovalent=monovalent,
            divalent=divalent, dntp=dntp)
    return tm


if __name__ == '__main__':
    # Basic test against primer3 (assert tm equality in 10000 random seqs)
    import primer3

    for _ in range(10000):
        seq = randomSeq(random.randint(5, 59))
        assert(primer3.calcTm(seq) == calcTm(seq))
    print('Passed test against primer3, 10000 trials')
