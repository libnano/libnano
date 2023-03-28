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
libnano.barcode_tools
~~~~~~~~~~~~~~~~~~~~~

Barcode generation and assignment methods

'''

import itertools
import json
import os.path as op
import random
from typing import (
    Dict,
    List,
    Optional,
)

import numpy as np  # type: ignore

from libnano import DATASET_DIR
from libnano.seqgraph import (  # type: ignore
    find_cliques,
    hamming_graph,
)

BARCODE_SETS_DIR = op.join(DATASET_DIR, 'barcode_sets')


def _get_cliques(
        seq_list: List[str],
        num_needed: int,
        min_hd: int = 2,
        cutoff: int = 1,
) -> List[List[int]]:
    '''Helper function for finding sequence groups w/ min inter-seq hd

    Args:
        seq_list: List of Python strings
        num_needed: Number of cliques needed to satisfy set
        min_hd: Minimum hamming distance between any two set members
        cutoff: Number of permissible solutions to allow

    Returns:
        Clique list of lists of integers
    '''
    hg = hamming_graph(seq_list)
    f = np.vectorize(lambda x: x[0])
    hg = f(hg)
    hd_thresh = np.zeros_like(hg)

    np.greater(
        hg,
        np.full_like(hg, min_hd - 1),
        hd_thresh,
    )
    return find_cliques(
        hd_thresh.astype(np.uint8),
        num_needed,
        cutoff,
    )


def get_barcode_set(
        set_name: str,
        dirpath: str = '',
) -> List[str]:
    '''Load a pre-computed barcode set from the `barcode_sets` folder

    Args:
        set_name: barcode set filename w/o .json extension
        dirpath: Optional directory path to load the JSON data from

    Returns:
        Barcode set list of string for the ``set_name``
    '''
    if not dirpath:
        dirpath = BARCODE_SETS_DIR
    try:
        with open(op.join(dirpath, f'{set_name}.json')) as fh:
            bc_set = json.load(fh)
    except OSError:
        raise ValueError(f'{set_name=} is not a valid set name')
    return bc_set


def _filter_barcode_list(
        bcs: List[str],
        bc_len: int,
) -> List[str]:
    '''Remove homopolymers of 4 or more and high GC content bcs

    Args:
        bcs: List of barcodes to filter
        bc_len: Desired length of a barcode sequence

    Returns:
        Filtered list of barcodes
    '''
    filtered_bcs = []
    homopols = ['AAAA', 'TTTT', 'GGGG', 'CCCC']
    for bc in bcs:
        homopol = False
        for h in homopols:
            if h in bc:
                homopol = True
                break
        if homopol:
            continue
        if (bc.count('G') + bc.count('C')) >= bc_len - 1:
            continue
        filtered_bcs.append(bc)
    return filtered_bcs


def _generate_initial_barcode_list(
        bc_len: int,
) -> List[str]:
    '''
    Args:
        bc_len: Desired length of a barcode sequence

    Returns:
        List of all {A, C, G, T} sequences of length ``bc_len`` post a
        homopolymers and GC content exclusive filter
    '''
    # self.min_hd = 0
    raw_bcs = [
        ''.join(bc) for bc in itertools.product(
            'ATGC',
            repeat=bc_len,
        )
    ]
    filtered_bcs = _filter_barcode_list(
        raw_bcs,
        bc_len,
    )
    random.shuffle(filtered_bcs)
    return filtered_bcs


class BarcodeGen(object):

    def __init__(
            self,
            bc_len: int = 5,
            initial_bcs: Optional[List[str]] = None,
    ):
        '''
        Args:
            bc_len: barcode length (nt)
            initial_bcs: initial barcode pool (if not provided, will be
                internally generated)
        '''
        self._bc_len = bc_len
        self._min_hd = 0
        self.bc_list: List[str] = (
            initial_bcs or _generate_initial_barcode_list(bc_len)
        )

    def __repr__(self):
        return (
            f'Set size of {self._bc_len} at {self._min_hd} hd '
            f'is: {len(self.bc_list)}'
        )

    def save(
            self,
            set_name: str = '',
            dirpath: str = '',
    ):
        '''Save data as a JSON file

        Args:
            set_name: Name of the JSON file. Defaults to name based on computed
                parameters
            dirpath: Directory to save to. Default is ``BARCODE_SETS_DIR``
        '''
        if not set_name:
            set_name = (
                f'{len(self.bc_list)}_{self._min_hd}hd_{self._bc_len}mer_v00'
            )
        if not dirpath:
            dirpath = BARCODE_SETS_DIR
        with open(op.join(dirpath, f'{set_name}.json')) as fh:
            json.dump(self.bc_list, fh)

    def get_random_barcode(self, do_remove: bool = False) -> str:
        '''Gets a barcode but does not remove it from the internal set

        Args:
            remove: If True, remove fetched value from the set.

        Returns:
            Random fetched barcode from :attr:`BarcodeGen.bc_list`.
        '''
        bc_seq = random.choice(self.bc_list)
        if do_remove:
            self.remove_barcode(bc_seq)
        return bc_seq

    def remove_barcode(self, bc_seq: str) -> None:
        '''Removes a barcode from the internal set

        Args:
            bc_seq: Barcode sequence to remove from internal set
        '''
        try:
            self.bc_list.remove(bc_seq)
        except ValueError:  # No longer in list
            pass

    def find_hamming_set(
            self,
            min_hd: int = 2,
            set_size: int = 50,
            cutoff: int = 1,
            raise_exc: bool = True,
    ) -> None:
        '''Find a subset of the current bc meeting the `min_hd` and `set_size`
        Updates :attr:`BarcodeGen.bc_list` with new set

        Args:
            min_hd: Minimum hamming distance between any two set members
            set_size: Minimum total set size
            cutoff: Number of permissible solutions to allow
            raise_exc: If True, whether or not to raise an exception on failure
        '''
        self._min_hd = min_hd
        cliques = _get_cliques(
            self.bc_list,
            set_size,
            min_hd,
            cutoff,
        )
        if len(cliques) < 1:
            if raise_exc:
                raise RuntimeError(
                    f'Could not find set of {set_size} barcodes of length '
                    f'{len(self.bc_list[0])} with min hamming distance '
                    f'{min_hd}',
                )
            return None
        else:
            hset = [self.bc_list[i] for i in cliques[0]]
            self.bc_list = hset
