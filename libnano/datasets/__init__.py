# Copyright (C) 2023. Nick Conway;
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
libnano.datasets.__init__
~~~~~~~~~~~~~~~~~~~~~~~~~

For now this only contains lazy-loaded access to the restriction enzyme
dataset, but it may act as the general access to libnano datasets should the
need arise in the future.

All datasets are provided in bytes (3.x) and str verions.
The unicode versions have a "_U" suffix.
'''

import json as _json
import os as _os
from typing import Dict

from libnano.datasets.build_codon_freq_dataset import (
    AA_TO_DNA,
    AA_TO_RNA,
    DNA_TO_AA,
    RNA_TO_AA,
    getOrganismFrequencies,
)
from libnano.datasets.build_enzyme_dataset import qcEnzymeDataset
from libnano.helpers import jsonbytes as _jsonbytes

_LOCAL_DIR: str = _os.path.dirname(_os.path.realpath(__file__))


class DatasetContainer(object):
    ''' Supports lazy loading of datasets to minimize memory footprint
    '''
    RNA_TO_AA = RNA_TO_AA
    DNA_TO_AA = DNA_TO_AA
    AA_TO_DNA = AA_TO_DNA
    AA_TO_RNA = AA_TO_RNA

    # ~~~~~~~~~~~~~~~~~~~~~~~ Codon frequency dataset ~~~~~~~~~~~~~~~~~~~~~~~ #

    def _loadCodonFreqDataset(
            self,
            as_unicode: bool = False,
    ) -> None:
        jsonlib = _json if as_unicode else _jsonbytes
        dataset_fp = _os.path.join(
            _LOCAL_DIR,
            'codon_freq_dataset.json',
        )
        with open(dataset_fp) as fd:
            raw_data = jsonlib.load(fd)
            if as_unicode:
                self._CODON_FREQ_DATASET_U = raw_data
            else:
                self._CODON_FREQ_DATASET = raw_data

    @property
    def CODON_FREQ_DATASET(self) -> dict:
        try:
            return self._CODON_FREQ_DATASET
        except AttributeError:
            self._loadCodonFreqDataset()
            return self._CODON_FREQ_DATASET

    @property
    def CODON_FREQ_DATASET_U(self) -> dict:
        try:
            return self._CODON_FREQ_DATASET_U
        except AttributeError:
            self._loadCodonFreqDataset(as_unicode=True)
            return self._CODON_FREQ_DATASET_U

    def organismFrequencies(
            self,
            organism_id: int,
            organism_class: str,
    ) -> Dict[str, Dict[str, float]]:
        return getOrganismFrequencies(
            organism_id,
            organism_class,
        )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Enzyme dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    def _loadEnzymeDataset(
            self,
            as_unicode: bool = False,
    ):
        jsonlib = _json if as_unicode else _jsonbytes
        dataset_fp = _os.path.join(_LOCAL_DIR, 'enzyme_dataset.json')
        with open(dataset_fp) as fd:
            raw_data = jsonlib.load(fd)
            if as_unicode:
                self._ENZYME_DATASET_U = raw_data[u'enzyme_data']
                self._REBASE_VERSION = raw_data[u'rebase_version']
            else:
                self._ENZYME_DATASET = raw_data[b'enzyme_data']
                self._REBASE_VERSION = raw_data[b'rebase_version']

    @property
    def ENZYME_DATASET(self) -> Dict:
        try:
            return self._ENZYME_DATASET
        except AttributeError:
            self._loadEnzymeDataset()
            return self._ENZYME_DATASET

    @property
    def REBASE_VERSION(self) -> Dict:
        try:
            return self._REBASE_VERSION
        except AttributeError:
            self._loadEnzymeDataset()
            return self._REBASE_VERSION

    @property
    def ENZYME_DATASET_U(self) -> Dict:
        try:
            return self._ENZYME_DATASET_U
        except AttributeError:
            self._loadEnzymeDataset(as_unicode=True)
            return self._ENZYME_DATASET_U


dataset_container = DatasetContainer()


__all__ = ['qcEnzymeDataset', 'dataset_container']
