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
libnano.search.restriction
~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
# import itertools
import re
# import sys
from enum import IntEnum
from typing import (
    Iterable,
    List,
    NamedTuple,
    Tuple,
    Union,
)

from libnano.datasets import dataset_container  # DatasetContainer,
from libnano.seqstr import complement  # type: ignore

# import pprint


STR_T = Union[str, bytes]
CUT_IDXS_T = Tuple[Tuple[int, int], ...]
REGEX_PAIR_T = Tuple[STR_T, STR_T]
REGEX_PAIR_COMP_T = Tuple['SRE_Pattern', 'SRE_Pattern']

class StrandDirEnum(IntEnum):
    Forward = 1
    Reverse = -1
cdef int FORWARD = 1
cdef int REVERSE = -1


class RestrictionMatch(NamedTuple):
    strand_dir: int
    start_idx: int
    end_idx: int
    enzyme: str
    cutside_number: int  # Index into the list of cutsites for the enzyme
    cut_idxs: CUT_IDXS_T
    pair_regex: str


# import pprint
# pprint.pprint(dataset_container.ENZYME_DATASET_U['M.BceSV'])

# ~~~~~~~~~~~~~~~ Restriction site search functions / classes ~~~~~~~~~~~~~~~ #

cdef class RestrictionSearcher:
    '''Base class for restriction enzyme searches. Instantiated with a list
    of restriction enzyme names and can then be used to search for presence,
    count or indicies of the respective cutsites in the provided sequence.
    '''

    cdef int _num_enzymes
    cdef object _str_type
    cdef tuple _enzyme_names
    cdef list _core_regexs
    cdef list _core_regexs_comp
    cdef list _full_regexs
    cdef list _full_regexs_comp
    cdef list _cutsite_idxs_list

    def __cinit__(
            self,
            *enzyme_names: Tuple[STR_T, ...],
    ):
        self._enzyme_names = enzyme_names
        self._num_enzymes = len(enzyme_names)
        self._str_type = type(enzyme_names[0])
        '''
        core_regexs == No ambiguous bases outside of core
        full_regexs == Includes ambiguous bases outside of core.  Usually we want this!
        '''
        cdef:
            list core_regexs = []          # type: List[REGEX_PAIR_T]
            list core_regexs_compiled = [] # type List[REGEX_PAIR_COMP_T]
            list full_regexs = []          # type: List[REGEX_PAIR_T]
            list full_regexs_compiled = [] # List[REGEX_PAIR_COMP_T]
            list cutsite_idxs_list = []    # type: List[CUT_IDXS_T]

        if isinstance(enzyme_names[0], str):
            coerce_type = lambda s: s
            enzyme_dataset: dict = dataset_container.ENZYME_DATASET_U
        else:
            coerce_type = lambda s: s.encode('utf-8')
            enzyme_dataset: dict = dataset_container.ENZYME_DATASET

        key_cutsites: STR_T = coerce_type('cutsites')
        key_core_regex: STR_T = coerce_type('core_regex')
        key_full_regex: STR_T = coerce_type('full_regex')
        key_cutsite_idxs: STR_T = coerce_type('cutsite_idxs')

        for en in enzyme_names:
            try:
                cutsites: List[dict] = enzyme_dataset[en][key_cutsites]
                for cs in cutsites:
                    # convert lists to tuples
                    str_temp = tuple(x for x in cs[key_core_regex])
                    core_regexs.append(str_temp)
                    core_regexs_compiled.append(tuple(re.compile(x) for x in str_temp))
                    str_temp = tuple(x for x in cs[key_full_regex])
                    full_regexs.append(str_temp)
                    full_regexs_compiled.append(tuple(re.compile(x) for x in str_temp))
                    cutsite_idxs_list.append(tuple(tuple(x) for x in cs[key_cutsite_idxs]))
            except KeyError:
                raise ValueError('%s is not a valid enzyme or is not present '
                                 'in the NEB Rebase database v.%d' % (en,
                                 dataset_container.REBASE_VERSION))
        self._core_regexs = core_regexs
        self._core_regexs_comp = core_regexs_compiled
        self._full_regexs = full_regexs
        self._full_regexs_comp = full_regexs_compiled
        self._cutsite_idxs_list: List[List[List[int]]] = cutsite_idxs_list

    @property
    def str_type(self):
        '''Check str_type of searcher, for testing purposes
        '''
        return self._str_type

    @property
    def enzyme_names(self) -> Tuple[STR_T, ...]:
        '''Tuple of enzyme names provided at instantiation

        '''
        return self._enzyme_names

    @enzyme_names.setter
    def enzyme_names(self, name_iterable: Iterable):
        self._num_enzymes = len(name_iterable)
        self._enzyme_names = tuple(name_iterable)

    @property
    def num_enzymes(self) -> int:
        '''Number of enzyme names provided at instantiation

        '''
        return self._num_enzymes

    def getCutSite(
            self,
            match,
    ) -> CUT_IDXS_T:
        return self._cutsite_idxs_list[match.cutsite_number]

    def sitesPresent(
            self,
            seq: STR_T,
            full_sites: bool = True,
    ) -> bool:
        '''
        Args:
            seq: sequence to check
            full_sites: If True, use the full regexes

        Returns:
            A boolean ``True`` or ``False`` if any of the restriction sites
            are present within `seq` or its reverse complement.

        Raises:
            TypeError: Cannot perform search on a `seq` of type provided

        '''
        regexs = self._full_regexs if full_sites else self._core_regexs
        try:
            for regex_fwd, regex_rev in regexs:
                if re.search(regex_fwd, seq) or re.search(regex_rev, seq):
                    return True
            return False
        except TypeError:
            raise TypeError(
                f'Cannot perform search on a `seq` of type {type(seq)}'
                ' using a RestrictionSearcher instantiated'
                f' with {self._str_type} type enzyme names'
            ) from None

    def countSites(
            self,
            seq: STR_T,
            full_sites: bool = True,
    ) -> List[int]:
        '''Return a list of counts indexed by enzyme (order as provided at
        instantiation). Palindromic sites will be counted twice (once per
        strand)

        Args:
            seq: sequence to check
            full_sites: If True, use the full regexes

        Returns:
            Return a list of counts indexed by enzyme (order as provided at
            instantiation). Palindromic sites will be counted twice (once per
            strand)

        Raises:
            TypeError: Cannot perform search on a `seq` of type provided

        '''
        regexs: List[REGEX_PAIR_T] = self._full_regexs if full_sites else self._core_regexs
        counts: List[int] = [0] * self._num_enzymes
        try:
            for idx, regex_pair in enumerate(regexs):
                counts[idx] +=  (   len(re.findall(regex_pair[0], seq)) +
                                    len(re.findall(regex_pair[1], seq))
                                )
            return counts
        except TypeError:
            raise TypeError(
                f'Cannot perform search on a `seq` of type {type(seq)}'
                ' using a RestrictionSearcher instantiated'
                f' with {self._str_type} type enzyme names'
            )

    def findSites(
            self,
            seq: STR_T,
            full_sites: bool = True,
    ) -> List[List[RestrictionMatch]]:
        '''
        Forward strand: 1
        Reverse strand: -1

        Indexing is from the 5' end of the forward strand.

        Args:
            seq: sequence to check
            full_sites: If True, use the full regexes

        Returns:
            A a list of lists, with the outer list indexed by enzyme
            (order as provided at instantiation) and will inner lists comprised
            of tuples indicating the (<strand>, <start index>, <end index>) of
            cutsites.

        Raises:
            TypeError: Cannot perform search on a `seq` of type provided

        '''
        cdef int idx = 0

        regexs = self._full_regexs if full_sites else self._core_regexs
        regexs_compiled = (
            self._full_regexs_comp if full_sites else self._core_regexs_comp
        )
        idx_list: List[List[RestrictionMatch]] = [
            [] for _ in range(self._num_enzymes)
        ]

        enzyme_names = self._enzyme_names

        cutsite_idxs_list = self._cutsite_idxs_list

        try:
            for regex_pair, regex_pair_compiled, in zip(
                regexs,
                regexs_compiled,
            ):
                regex_fwd, regex_rev = regex_pair
                regex_fwd_comp, regex_rev_comp = regex_pair_compiled
                idx_list[idx] += [
                    RestrictionMatch(
                            StrandDirEnum.Forward,
                            m.start(),
                            m.end(),
                            enzyme_names[idx],
                            idx,
                            cutsite_idxs_list[idx],
                            regex_rev
                    )
                    for m in regex_fwd_comp.finditer(seq)
                ]
                idx_list[idx] += [
                    RestrictionMatch(
                        StrandDirEnum.Reverse,
                        m.start(),
                        m.end(),
                        enzyme_names[idx],
                        idx,
                        cutsite_idxs_list[idx],
                        regex_fwd
                    )
                    for m in regex_rev_comp.finditer(seq)
                ]
                idx += 1
            return idx_list
        except TypeError:
            raise TypeError(
                f'Cannot perform search on a `seq` of type {type(seq)}'
                ' using a RestrictionSearcher instantiated'
                f' with {self._str_type} type enzyme names'
            )

    def __repr__(self):
        return '<libnano.restriction.RestrictionSearcher> for enzymes: %s' %  \
                ', '.join(self._enzyme_names)

    def __str__(self):
        return self.__repr__()

    def displaySites(
            self,
            seq: STR_T,
            sites: List[List[RestrictionMatch]],
    ):
        '''Print sites to the terminal

        Args:
            seq: Sequence to search
            sites: Sites to display
        '''
        print(seq)
        complement_seq: str = complement(seq)
        print(complement_seq)
        for en, en_sites in zip(self._enzyme_names, sites):
            print(en)
            print(f'Found {len(en_sites)} sites')
            # print(en_sites)
            for rs_match in en_sites:
                if rs_match.strand > 0:
                    print("forward")
                    print(
                        seq[
                            rs_match.start_idx:
                            rs_match.end_idx+1
                        ]
                    )
                else:
                    print("reverse")
                    print(
                        complement_seq[
                            rs_match.start_idx:
                            rs_match.end_idx+1
                        ]
                    )
            print('')


def getEnzymeRegexs(
        enzyme_names: List[STR_T],
        full_sites: bool = True,
) -> Tuple[List[STR_T], List[CUT_IDXS_T]]:
    '''Get the enzyme regexs for the provided enzyme names. If `full_sites` is
    True, the regexs will be for the full restriction sites, otherwise the
    regexs will be for the core restriction sites.

    Args:
        enzyme_names: List of enzyme names
        full_site: If True, search using full regex

    Returns:
        Tuple of the form::

            <regexes matched>, <cutsite_idxs_list>

    Raises:
        ValueError: Not a valid enzyme or is not present in the NEB Rebase DB
    '''
    cdef list regexs = []
    cdef list cutsite_idxs_list = []

    if isinstance(enzyme_names[0], str):
        coerce_type = lambda s: s
        enzyme_dataset: dict = dataset_container.ENZYME_DATASET_U
    else:
        coerce_type = lambda s: s.encode('utf-8')
        enzyme_dataset: dict = dataset_container.ENZYME_DATASE

    key_core_regex = coerce_type('core_regex')
    key_full_regex = coerce_type('full_regex')
    key_cutsites = coerce_type('cutsites')
    key_cutsite_idxs = coerce_type('cutsite_idxs')

    key_use_regex = key_full_regex if full_sites else key_core_regex

    for en in enzyme_names:
        cutsites: List[dict] = enzyme_dataset[en][key_cutsites]
        try:
            for cs in cutsites:
                # convert lists to tuples
                regexs.append(
                    tuple(x for x in cs[key_use_regex]),
                )
                cutsite_idxs_list.append(
                    tuple(tuple(x) for x in cs[key_cutsite_idxs]),
                )
        except KeyError:
            raise ValueError(
                f'{en} is not a valid enzyme or is not present '
                'in the NEB Rebase database v.'
                f'{dataset_container.REBASE_VERSION}'
            ) from None
    return regexs, cutsite_idxs_list


def sitesPresent(
        seq: STR_T,
        enzyme_names: List[STR_T],
        full_sites: bool = True,
) -> bool:
    '''
    Args:
        enzyme_names: List of enzyme names
        full_site: If True, search using full regex

    Returns:
        A boolean True or False if any of the restriction sites
        are present within `seq` or its reverse complement.

    Raises:
        TypeError: Cannot perform search on a `seq` of type provided
    '''
    regexs, _ = getEnzymeRegexs(enzyme_names, full_sites=full_sites)
    try:
        for regex in regexs:
            if re.search(regex[0], seq) or re.search(regex[1], seq):
                return True
        return False
    except TypeError:
        raise TypeError(
            f'Cannot perform search on a `seq` of type {type(seq)}'
            f' using enzyme names of type {type(enzyme_names[0])}'
        )


def countSites(
        seq: STR_T,
        enzyme_names: List[STR_T],
        full_sites: bool = True,
) -> List[int]:
    '''
    Args:
        seq: Sequence to inspect
        enzyme_names: List of enzyme names
        full_site: If True, search using full regex

    Returns:
        A list of counts indexed by enzyme (order as provided at
        instantiation). Palindromic sites will be counted twice (once per
        strand)

    Raises:
        TypeError: Cannot perform search on a `seq` of type provided
    '''
    regexs: List[STR_T]
    regexs, _ = getEnzymeRegexs(
        enzyme_names,
        full_sites=full_sites,
    )
    counts: List[int] = [0] * len(enzyme_names)
    try:
        for idx, regex_pair in enumerate(regexs):
            counts[idx] += (
                len(re.findall(regex_pair[0], seq)) +
                len(re.findall(regex_pair[1], seq))
            )
        return counts
    except TypeError:
        raise TypeError(
            f'Cannot perform search on a `seq` of type {type(seq)}'
            f' using enzyme names of type {type(enzyme_names[0])}'
        ) from None


def findSites(
        seq: STR_T,
        enzyme_names: List[STR_T],
        full_sites: bool = True,
) -> List[List[RestrictionMatch]]:
    '''Return a list of lists, with the outer list indexed by enzyme
    (order as provided at instantiation) and will inner lists comprised
    of tuples indicating the (<strand>, <start index>, <end index>) of
    cutsites.

    Forward strand: 1
    Reverse strand: -1

    Indexing is from the 5' end of the forward strand.

    Args:
        seq: the sequence string
        enzyme_names: List of enzyme names known in the database
        full_sites: If True, search using full regex

    Returns:
        List of Lists of :obj:`RestrictionMatch` ``NamedTuple``s

    Raises:
        TypeError: Cannot perform search on a `seq` of type provided
    '''
    regexs, cutsite_idxs_list = getEnzymeRegexs(
        enzyme_names,
        full_sites=full_sites,
    )
    idx_list: List[List[RestrictionMatch]] = [
        [] for _ in range(len(enzyme_names))
    ]
    try:
        for idx, regex_pair in enumerate(regexs):
            idx_list[idx] += [
                RestrictionMatch(
                    StrandDirEnum.Forward,
                    m.start(),
                    m.end(),
                    enzyme_names[idx],
                    idx,
                    cutsite_idxs_list[idx],
                    regex_pair[1]
                )
                for m in re.finditer(regex_pair[0], seq)
            ]
            idx_list[idx] += [
                RestrictionMatch(
                    StrandDirEnum.Reverse,
                    m.start(),
                    m.end(),
                    enzyme_names[idx],
                    idx,
                    cutsite_idxs_list[idx],
                    regex_pair[0]
                )
                for m in re.finditer(regex_pair[1], seq)
            ]
        return idx_list
    except TypeError:
        raise TypeError(
            f'Cannot perform search on a `seq` of type {type(seq)}'
            f' using enzyme names of type {type(enzyme_names[0])}'
        ) from None
