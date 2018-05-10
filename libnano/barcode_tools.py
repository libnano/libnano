# -*- coding: utf-8 -*-
'''libnano.barcode_tools

Barcode generation and assignment methods

'''

import copy
import json
import random
import itertools
import os

import numpy as np
from libnano.seqgraph import hammingGraph, find_cliques
from libnano.seqstr import minHammingDistance

from libnano import DATASET_DIR

BARCODE_SETS_DIR = os.path.join(DATASET_DIR, 'barcode_sets')


def _getCliques(seq_list, num_needed, min_hd=2, cutoff=1):
    '''Helper function for finding sequence groups w/ min inter-seq hd
    '''
    hg = hammingGraph(seq_list)
    f = np.vectorize(lambda x: x[0])
    hg = f(hg)
    hd_thresh = np.zeros_like(hg)
    np.greater(hg, np.full_like(hg, min_hd-1), hd_thresh)
    return find_cliques(hd_thresh.astype(np.uint8), num_needed, cutoff)


def getBarcodeSet(set_name, path=None):
    '''Load a pre-computed barcode set from the `barcode_sets` folder

    Args:
        set_name (str): barcode set filename w/o .json extension
    '''
    if path is None:
        path = BARCODE_SETS_DIR
    try:
        with open(os.path.join(path, set_name + '.json')) as fh:
            bc_set = json.load(fh)
    except OSError:
        raise ValueError('"%s" is not a valid set name')
    return bc_set


class BarcodeGen(object):

    def __init__(self, bc_len=5, initial_bcs=None):
        '''
        Args:
            bc_len (int): barcode length (nt)
            initial_bcs (list of str): initial barcode pool (if not provided,
                                       will be internally generated)
        '''
        self._bc_len = bc_len
        self._min_hd = 0
        self.bcs = initial_bcs or self._genInitialBcs(bc_len)

    def __repr__(self):
        return "Set size of %s at %s hd is: %s" % (self._bc_len, self._min_hd, len(self.bcs))

    def save(self, set_name=None, path=None):
        if set_name is None:
            set_name = "%d_%dhd_%dmer_v00" % (len(self.bcs), self._min_hd, self._bc_len)
        if path is None:
            path = BARCODE_SETS_DIR
        with open(os.path.join(path, set_name + '.json')) as fh:
            json.dump(self.bcs, fh)

    def getBc(self, remove=False):
        '''Gets a barcode but does not remove it from the internal set
        '''
        bc_seq = random.choice(self.bcs)
        if remove:
            self.removeBc(bc_seq)
        return bc_seq

    def removeBc(self, bc_seq):
        '''Removes a barcode from the internal set
        '''
        try:
            self.bcs.remove(bc_seq)
        except ValueError:  # No longer in list
            pass

    def findHammingSet(self, min_hd=2, set_size=50, cutoff=1, raise_exc=True):
        '''Find a subset of the current bc meeting the `min_hd` and `set_size`

        Args:
            min_hd (int): min hamming distance between any two set members
            set_size (int): min total set size
            raise_exc (bool): whether or not to raise an exception on failure
        '''
        self._min_hd = min_hd
        cliques = _getCliques(self.bcs, set_size, min_hd, cutoff)
        if len(cliques) < 1:
            if raise_exc:
                raise RuntimeError('Could not find set of %d barcodes of '
                                   'length %d with min hamming distance %d'
                                   % (set_size, len(self.bcs[0]), min_hd))
            return None
        else:
            hset = [self.bcs[i] for i in cliques[0]]
            self.bcs = hset

    def _genInitialBcs(self, bc_len):
        self.min_hd = 0
        raw_bcs = [''.join(bc) for bc in itertools.product('ATGC',
                                                           repeat=bc_len)]
        filtered_bcs = self._filterBcs(raw_bcs, bc_len)
        random.shuffle(filtered_bcs)
        return filtered_bcs

    def _filterBcs(self, bcs, bc_len):
        '''Remove homopolymers of 4 or more and high GC content bcs
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
