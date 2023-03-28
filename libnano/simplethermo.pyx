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
libnano.simplethermo
~~~~~~~~~~~~~~~~~~~~

Basic thermodynamic calculations

'''

import math as _math
from typing import Tuple

from libnano.seqstr import reverse_complement as _rc  # type: ignore

# ~~~~~~~~~~~~~~~~~~~~ GENERAL CONSTANTS AND CALCULATIONS ~~~~~~~~~~~~~~~~~~~ #

cdef:
    double R = 1.9872e-3  # Gas constant (1 / kcal * mol)
    double KELVIN = 273.15


def c_to_k(deg_c: float) -> float:
    '''Convert degrees Celsius to degrees Kelvin

    Args:
        deg_c (float): temperature in [Celsius]

    Returns:
        float: temperature in [Kelvin]

    '''
    return deg_c + KELVIN


def k_to_c(deg_k: float) -> float:
    ''' Convert degrees Kelvin to degrees Celsius

    Args:
        deg_k (float): temperature in [Kelvin]

    Returns:
        float: temperature in [Celsius]

    '''
    return deg_k - KELVIN


cdef double calc_ka(
        double dg,
        double deg_c,
):
    '''Return the association constant at a given temperature

    Args:
        dg (float): the free energy in [kcal/mol]
        deg_c (float): temperature in [Celsius]

    Returns:
        float: the association constant at a given temperature

    '''
    return _math.exp**(-dg/(R * c_to_k(deg_c)))


def calc_dg(
        ds: float,
        dh: float,
        deg_c: float,
) -> float:
    '''Return the dg at a given temp using the provided ds and dh

    Args:
        ds (float): the entropy in [kcal/mol]
        dh (float): enthalpy in in [kcal/mol]
        deg_c (float): temperature in [Celsius]

    Returns:
        float: the dg at a given temp using the provided ds and dh

    '''
    return dh - ds * c_to_k(deg_c)


def calc_rand_coil(
        dg: float,
        deg_c: float,
) -> float:
    ''' Return the percent of randomly coiled oligo with dg at deg_c degrees

    Args:
        dg (float): the free energy in [kcal/mol]
        deg_c (float): temperature in [Celsius]

    Returns:
        float: the percent of randomly coiled oligo with dg at deg_c degrees
    '''
    return 1 / (calc_ka(dg, deg_c) + 1)


# ~~~~~~~~~~~~~~~~~ PYTHON EQUIVALENT OF PRIMER3 OLIGOTM  ~~~~~~~~~~~~~~~~ #

cdef float divalent_to_monovalent(
        float divalent,
        float dntp,
):
    '''Calculate equivalent monovalent effect for a given
    divalent parameters

    Args:
        divalent (float): Divalent concentration
        dntp (float): DNTP concetration

    Returns:
        float: equivalent monovalent concentration percent
    '''
    if divalent == 0.:
        dntp = 0.
    if divalent < dntp:
        divalent = dntp
    return 120. * _math.sqrt(divalent - dntp)

def calc_thermo(
        seq: str,
        conc_nm: float = 50,
        monovalent: float = 50,
        divalent: float = 0.01,
        dntp: float = 0.0,
) -> Tuple[float, float, float]:
    ''' Return the thermo parameters for DNA under specified salt cond.
        Data is from referenced from PRIMER3

    Args:
        seq (str): the sequence to analyze
        conc_nm (Optional[int]): percent concentration
        monovalent (Optional[int]): percent concentration
        divalent (Optional[float]): fractional concentration of Mg2+
        dntp (Optional[float]): fractional concentration dntp

    Returns:
        Tuple[float]: (dH, dS, Tm)
    '''
    cdef:
        float dH = 0.
        float dS = 0.
        float monovalent_use, tm
        Py_ssize_t idx
        # Calculate oligo symmetry
        bint sym = seq == _rc(seq)

    monovalent_use = monovalent
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
    # Calculate NN uncorrected dS and dH for oligo
    for idx in range(len(seq) - 1):
        dH += enthalpies[seq[idx:idx + 2]]
        dS += entropies[seq[idx:idx + 2]]
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
    monovalent_use += divalent_to_monovalent(divalent, dntp)
    dS = dS + 0.368 * (len(seq) - 1) * _math.log(monovalent / 1000.0)
    # Account for oligo symmetry and calculate tm
    if sym:
        tm = dH / (dS + 1.987 * _math.log(conc_nm/1.0e9)) - KELVIN
    else:
        tm = dH / (dS + 1.987 * _math.log(conc_nm/4.0e9)) - KELVIN
    return dH, dS, tm

def calc_tm(
        seq: str,
        conc_nm: float = 50,
        monovalent: float = 50,
        divalent: float = 0.01,
        dntp: float = 0.0,
):
    _, _, tm = calc_thermo(
        seq,
        conc_nm=conc_nm,
        monovalent=monovalent,
        divalent=divalent,
        dntp=dntp,
    )
    '''Same as `calcThermo` but returns just the melting temperature

    Args:
        seq (str): the sequence to analyze
        conc_nm (Optional[int]): percent concentration
        monovalent (Optional[int]): percent concentration
        divalent (Optional[float]): fractional concentration of Mg2+
        dntp (Optional[float]): fractional concentration dntp

    Returns:
        float: Melting temperature in [Celcius]
    '''
    return tm
