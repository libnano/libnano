
import itertools
import re
import sys
from typing import (
    Union,
    List
)

from libnano.datasets import dataset_container

STR_T = Union[str, bytes]

# ~~~~~~~~~~~~~~~ Restriction site search functions / classes ~~~~~~~~~~~~~~~ #

cdef class RestrictionSearcher:
    '''Base class for restriction enzyme searches. Instantiated with a list
    of restriction enzyme names and can then be used to search for presence,
    count or indicies of the respective cutsites in the provided sequence.
    '''

    cdef object _num_enzymes
    cdef object _str_type
    cdef object _enzyme_names
    cdef object _core_regexs
    cdef object _full_regexs

    def __cinit__(self, *enzyme_names):
        self._enzyme_names = enzyme_names
        self._num_enzymes = len(enzyme_names)
        self._str_type = type(enzyme_names[0])
        cdef object core_regexs = []
        cdef object full_regexs = []
        if isinstance(enzyme_names[0], str):
            coerce_type = lambda s: s
            enzyme_dataset = dataset_container.ENZYME_DATASET_U
        else:
            coerce_type = lambda s: s.encode('utf-8')
            enzyme_dataset = dataset_container.ENZYME_DATASET

        key_cutsites = coerce_type('cutsites')
        key_core_regex = coerce_type('core_regex')
        key_full_regex = coerce_type('full_regex')
        for en in enzyme_names:
            try:
                cutsites = enzyme_dataset[en][key_cutsites]
                core_regexs += itertools.chain.from_iterable(
                    [cs[key_core_regex] for cs in cutsites])
                full_regexs += itertools.chain.from_iterable(
                    [cs[key_full_regex] for cs in cutsites])
            except KeyError:
                raise ValueError('%s is not a valid enzyme or is not present '
                                 'in the NEB Rebase database v.%d' % (en,
                                 dataset_container.REBASE_VERSION))
        self._core_regexs = core_regexs
        self._full_regexs = full_regexs

    property enzyme_list:
        '''Public-facing enzyme list '''
        def __get__(RestrictionSearcher self):
            return self._enzyme_names

    property str_type:
        '''Check str_type of searcher, for testing purposes '''
        def __get__(RestrictionSearcher self):
            return self._str_type

    property _enzyme_names:
        '''List of enzyme names provided at instantiation '''
        def __set__(RestrictionSearcher self, object name_list):
            self._enzyme_names = name_list

        def __get__(RestrictionSearcher self):
            return self._enzyme_names

    property _num_enzymes:
        '''List of enzyme names provided at instantiation '''
        def __set__(RestrictionSearcher self, object num_enzymes):
            self._num_enzymes = num_enzymes

        def __get__(RestrictionSearcher self):
            return self._num_enzymes

    property _core_regexs:
        '''List of core enzyme recognition seq. regexs '''
        def __set__(RestrictionSearcher self, object regex_list):
            self._core_regexs = [re.compile(regex) for regex in regex_list]

        def __get__(RestrictionSearcher self):
            return self._core_regexs

    property _full_regexs:
        '''List of full enzyme recognition seq. regexs '''
        def __set__(RestrictionSearcher self, object regex_list):
            self._full_regexs = [re.compile(regex) for regex in regex_list]

        def __get__(RestrictionSearcher self):
            return self._full_regexs

    def sitesPresent(RestrictionSearcher self,
                    seq: STR_T,
                    full_sites: bool = True) -> bool:
        '''Return a boolean ``True`` or ``False`` if any of the restriction sites
        are present within `seq` or its reverse complement.
        '''
        regexs = self._full_regexs if full_sites else self._core_regexs
        try:
            for regex in regexs:
                if re.search(regex, seq):
                    return True
            return False
        except TypeError:
            raise TypeError('Cannot perform search on a `seq` of type {}'
                            ' using a RestrictionSearcher instantiated'
                            ' with {} type enzyme names'.format(type(seq),
                            self._str_type))

    def countSites(RestrictionSearcher self,
                    seq: STR_T,
                    full_sites: bool = True) -> int:
        '''Return a list of counts indexed by enzyme (order as provided at
        instantiation). Palindromic sites will be counted twice (once per
        strand)
        '''
        regexs = self._full_regexs if full_sites else self._core_regexs
        counts = [0] * self._num_enzymes
        try:
            for idx, regex in enumerate(regexs):
                counts[idx//2] += len(re.findall(regex, seq))
            return counts
        except TypeError:
            raise TypeError('Cannot perform search on a `seq` of type {}'
                            ' using a RestrictionSearcher instantiated'
                            ' with {} type enzyme names'.format(type(seq),
                            self._str_type))

    def findSites(RestrictionSearcher self,
                    seq: STR_T,
                    full_sites: bool = True) -> List[List[Tuple[int, int, int]]]:
        '''Return a list of lists, with the outer list indexed by enzyme
        (order as provided at instantiation) and will inner lists comprised
        of tuples indicating the (<strand>, <start index>, <end index>) of
        cutsites.

        Forward strand: 1
        Reverse strand: -1

        Indexing is from the 5' end of the forward strand.
        '''
        seq_len = len(seq)
        regexs = self._full_regexs if full_sites else self._core_regexs
        idx_list = [[] for _ in range(self._num_enzymes)]
        try:
            for idx, regex in enumerate(regexs):
                idx_list[idx//2] += [(idx%2 * -2 + 1, m.start(), m.end())
                                     for m in re.finditer(regex, seq)]
            return idx_list
        except TypeError:
            raise TypeError('Cannot perform search on a `seq` of type {}'
                            ' using a RestrictionSearcher instantiated'
                            ' with {} type enzyme names'.format(type(seq),
                            self._str_type))

    def __repr__(RestrictionSearcher self):
        return '<libnano.restriction.RestrictionSearcher> for enzymes: %s' %  \
                ', '.join(self._enzyme_names)

    def __str__(RestrictionSearcher self):
        return self.__repr__()


def getEnzymeRegexs(enzyme_names: STR_T, full_sites: bool = True) -> STR_T:
    '''Get the enzyme regexs for the provided enzyme names. If `full_sites` is
    True, the regexs will be for the full restriction sites, otherwise the
    regexs will be for the core restriction sites.
    '''
    cdef object regexs = []
    cdef object lookup = 'full_regex' if full_sites else 'core_regex'
    if isinstance(enzyme_names[0], str):
        enzyme_dataset = dataset_container.ENZYME_DATASET_U
    else:
        enzyme_dataset = dataset_container.ENZYME_DATASET
    for en in enzyme_names:
        try:
            regexs += itertools.chain.from_iterable(
                [cs[lookup] for cs in enzyme_dataset[en]['cutsites']])
        except KeyError:
            raise ValueError('%s is not a valid enzyme or is not present '
                             'in the NEB Rebase database v.%d' % (en,
                             dataset_container. REBASE_VERSION))
    return regexs


def sitesPresent(seq: STR_T,
                enzyme_names: List[STR_T],
                full_sites: bool = True) -> bool:
    '''Return a boolean True or False if any of the restriction sites
    are present within `seq` or its reverse complement.
    '''
    regexs = getEnzymeRegexs(enzyme_names, full_sites=full_sites)
    try:
        for regex in regexs:
            if re.search(regex, seq):
                return True
        return False
    except TypeError:
        raise TypeError('Cannot perform search on a `seq` of type {}'
                        ' using enzyme names of type {}'.format(type(seq),
                        type(enzyme_names[0])))



def countSites(seq: STR_T,
                enzyme_names: List[STR_T],
                full_sites: bool = True) -> int:
    '''Return a list of counts indexed by enzyme (order as provided at
    instantiation). Palindromic sites will be counted twice (once per
    strand)
    '''
    regexs = getEnzymeRegexs(enzyme_names, full_sites=full_sites)
    counts = [0] * len(enzyme_names)
    try:
        for idx, regex in enumerate(regexs):
            counts[idx//2] += len(re.findall(regex, seq))
        return counts
    except TypeError:
        raise TypeError('Cannot perform search on a `seq` of type {}'
                        ' using enzyme names of type {}'.format(type(seq),
                        type(enzyme_names[0])))


def findSites(seq: STR_T,
                enzyme_names: List[STR_T],
                full_sites: bool = True) -> List[List[Tuple[int, int, int]]]:
    '''Return a list of lists, with the outer list indexed by enzyme
    (order as provided at instantiation) and will inner lists comprised
    of tuples indicating the (<strand>, <start index>, <end index>) of
    cutsites.

    Forward strand: 1
    Reverse strand: -1

    Indexing is from the 5' end of the forward strand.
    '''
    seq_len = len(seq)
    regexs = getEnzymeRegexs(enzyme_names, full_sites=full_sites)
    idx_list = [[] for _ in range(len(enzyme_names))]
    try:
        for idx, regex in enumerate(regexs):
            idx_list[idx//2] += [(idx%2 * -2 + 1, m.start(), m.end()) for m in
                                 re.finditer(regex, seq)]
        return idx_list
    except TypeError:
        raise TypeError('Cannot perform search on a `seq` of type {}'
                        ' using enzyme names of type {}'.format(type(seq),
                        type(enzyme_names[0])))
