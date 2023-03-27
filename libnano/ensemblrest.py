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
libnano.ensemblrest
~~~~~~~~~~~~~~~~~~~
'''
import atexit
import copy
import io
import json
import os.path as op
import pickle
import sys
from pathlib import Path
from pprint import pprint
from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    Union,
)

import pandas as pd
import requests

from libnano.seqstr import reverseComplement  # type: ignore

USE_CACHE: bool = True
DO_PRINT_CACHE: bool = False
TIMEOUT_FAST: float = 2.0
TIMEOUT_SLOW: float = 6.0
TIMEOUT_REQU: float = 5 * TIMEOUT_SLOW

SERVER: str = 'https://rest.ensembl.org'

EXON_KEYS: List[str] = [
    'assembly_name',
    'db_type',
    'end',
    'id',
    'object_type',
    'seq_region_name',
    'species',
    'start',
    'strand',
    'version',
]
TRANSLATION_KEYS: List[str] = [
    'Parent',
    'db_type',
    'end',
    'id',
    'length',
    'object_type',
    'species',
    'start',
]
UTR_KEYS: List[str] = ['Parent'] + EXON_KEYS
TRANSCRIPT_KEYS: List[str] = [
    'Exon',
    'Parent',
    'Translation',
    'UTR',
    'assembly_name',
    'biotype',
    'db_type',
    'display_name',
    'end',
    'id',
    'is_canonical',
    'logic_name',
    'object_type',
    'seq_region_name',
    'source',
    'species',
    'start',
    'strand',
    'version',
]
LOOKUP_KEYS: List[str] = [
    'Transcript',
    'assembly_name',
    'biotype',
    'db_type',
    'description',
    'display_name',
    'end',
    'id',
    'logic_name',
    'object_type',
    'seq_region_name',
    'source',
    'species',
    'start',
    'strand',
    'version',
]


class Probe(NamedTuple):
    end: int
    feature_type: str
    microarray: str
    probe_length: int
    probe_name: str
    probe_set: str
    seq_region_name: str
    start: int
    strand: int


PROBE_KEYS: Tuple[Any, ...] = Probe._fields

SPECIES_NAMES: Dict = {'mouse': ['mus_musculus']}
ASSEMBLY_NAMES: Dict = {'mouse': ['GRCm38']}

DB_NAMES: List[str] = ['EntrezGene', 'MGI', 'Uniprot_gn', 'WikiGene']

OBJECT_TYPES: List[str] = [
    'Translation',
    'Exon',
    'five_prime_UTR',
    'three_prime_UTR',
    'Transcript',
]
# From overlap(x, feature='variation') call
CONSEQUENCE_TYPES: List[str] = [
    '3_prime_UTR_variant',
    'intron_variant',
    'synonymous_variant',
    'non_coding_transcript_exon_variant',
    'splice_region_variant',
    '5_prime_UTR_variant',
]
POST_JSON: Dict = {
    'Content-Type': 'application/json',
    'Accept': 'application/json',
}

home_path: str = str(Path.home())
CACHE_FILE: str = op.join(home_path, '.ENSEMBLCACHE.pickle')

SPECIES_LIST: List[str] = ['mouse', 'human']
_cache_dirty: bool
ensembl_cache: Dict
THE_FILE: str = op.basename(__file__)


def makeCache(
        species_list: Optional[List[str]] = None,
) -> Dict[str, Any]:
    global _cache_dirty
    global SPECIES_LIST
    if species_list is None:
        species_list = SPECIES_LIST
    d: Dict[str, Any] = {
        species: {}
        for species in species_list
    }
    d['species_list'] = species_list
    SPECIES_LIST = species_list
    _cache_dirty = False
    return d


def loadCache(filename: str) -> Dict:
    try:
        with io.open(filename, 'rb') as fd:
            the_cache: Dict = pickle.load(fd)
    except OSError:
        if DO_PRINT_CACHE:
            print(
                f'Could not load cache for {THE_FILE}: {filename}',
            )
        the_cache = makeCache()
    if DO_PRINT_CACHE:
        print(f'LOADED cache for {THE_FILE}: {CACHE_FILE}')
    return the_cache


if Path(CACHE_FILE).exists():
    ensembl_cache = loadCache(CACHE_FILE)
    SPECIES_LIST = ensembl_cache['species_list']
    _cache_dirty = False
else:
    ensembl_cache = makeCache()


def clearCache(
        species_list: Optional[List[str]] = None,
):
    global _cache_dirty
    global ensembl_cache
    _cache_dirty = False
    ensembl_cache = makeCache(
        species_list=species_list,
    )


def addSpecies(
        species: str,
) -> None:
    global SPECIES_LIST
    global _cache_dirty
    if species not in SPECIES_LIST:
        SPECIES_LIST.append(species)
        ensembl_cache[species] = {}
        _cache_dirty = True


def _closeCache(
) -> None:
    global _cache_dirty
    global ensembl_cache
    if _cache_dirty:
        with io.open(CACHE_FILE, 'wb') as fd:
            pickle.dump(ensembl_cache, fd)
        if DO_PRINT_CACHE:
            print(f'UPDATED {THE_FILE} cache')


atexit.register(_closeCache)


def getCache(
        url: str,
        arg: str,
        cache: Dict[str, Any],
        content_type: str = 'application/json',
) -> Any:
    global _cache_dirty
    global USE_CACHE

    res = cache.get(arg) if USE_CACHE else None

    if res is None:
        res = getURL(url, content_type=content_type)
        cache[arg] = res
        _cache_dirty = True
    return res


def getCacheList(
        url: str,
        data: Dict,
        arg_list: List[str],
        cache: Dict[str, Any],
) -> Dict[Any, Any]:
    global _cache_dirty
    global USE_CACHE
    do_post: bool = False
    res: Dict[Any, Any] = {}
    if USE_CACHE:
        for x in arg_list:
            try:
                res[x] = cache[x]
            except KeyError:
                do_post = True
                break
    else:
        do_post = True
    if do_post:
        res = postURL(url, data=data)
        cache.update(res)
        _cache_dirty = True
    return res


class Base:
    key_list: List[str] = []

    def __init__(self, d: Dict):
        self.d = d

    def __getattr__(self, key: str) -> Any:
        if key in self.key_list:
            return self.d[key]

    @property
    def idxs(self) -> Tuple[int, int]:
        return self.d['start'], self.d['end']

    @property
    def is_fwd(self) -> bool:
        return self.d['strand'] == 1

    @property
    def is_rev(self) -> bool:
        return self.d['strand'] == -1

    @property
    def chromosome(self) -> int:
        return int(self.d['seq_region_name'])


class Exon(Base):
    key_list: List[str] = EXON_KEYS


class Transcript(Base):
    key_list: List[str] = TRANSCRIPT_KEYS


class LookUp(Base):
    key_list: List[str] = LOOKUP_KEYS

    @property
    def transcripts(self) -> List[Transcript]:
        return self.d['Transcript']


def getURL(
        url: str,
        content_type='application/json',
) -> Any:
    '''GET request to the ``url`` and return the JSON response
    Exits if response is not OK

    Args:
        url: URL to post to
        headers: Headers to GET with

    Returns:
        JSON reposnse
    '''
    r = requests.get(
        url,
        headers={'Content-Type': content_type},
        timeout=TIMEOUT_REQU,
    )
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    if content_type == 'text/plain':
        return r.text
    else:
        return r.json()


def postURL(
        url: str, data: Dict,
        headers: Dict = POST_JSON,
) -> Dict:
    '''Post request to the ``url`` and return the JSON response
    Exits if response is not OK

    Args:
        url: URL to post to
        headers: Headers to post with

    Returns:
        JSON reposnse
    '''
    data_send = json.dumps(data)
    r = requests.post(
        url,
        headers=headers,
        data=data_send,
        timeout=TIMEOUT_REQU,
    )
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.json()


def archiveID(eid: str) -> Dict:
    '''Uses the given identifier to return the archived sequence

    Args:
        eid: Exon ID

    Returns:
        Dictionary of the archive ID given the Exon ID
    '''
    url = f'{SERVER}/archive/id/{eid}?'
    return getURL(url)  # type: ignore


def xRefName(species: str, name: str) -> List[Dict]:
    '''Performs a lookup based upon the primary accession or display label of an
     external reference and returning the information we hold about the entry

    Args:
        species: Species name
        name: name of gene

    Returns:
        List of dictionaries from the ``SERVER`` url

    '''
    url = f'{SERVER}/xrefs/name/{species}/{name}?'
    return getURL(url)  # type: ignore


def externalXRefs(
        species: str,
        name: str,
) -> Dict[str, List[Dict[str, Any]]]:
    '''
    Args:
        species: Species name
        name: name of gene

    Returns:
        Dictionary keyed by ``dbname`` of Lists of items with the format::

            {
                'description': ...,
                'display_id': ...,
                'info_type': ...,
                'info_text': ...,
                'primary_id': ...,
                'synonyms': ...,
                'version': ...,
            }

    '''
    res: List[Dict] = xRefName(species, name)
    out: Dict[str, List[Dict[str, Any]]] = {}
    for item in res:
        dbname = item['dbname']
        if out.get(dbname) is None:
            out[dbname] = []
        out[dbname].append({
            'description': item['description'],
            'display_id': item['display_id'],
            'info_type': item['info_type'],
            'info_text': item['info_text'],
            'primary_id': item['primary_id'],
            'synonyms': item['synonyms'],
            'version': item['version'],
        })
    return out


def xRefSymbol(species: str, name: str) -> List[Dict]:
    '''Looks up an external symbol and returns all Ensembl objects linked to
    it. This can be a display name for a gene/transcript/translation, a synonym
     or an externally linked reference. If a gene's transcript is linked to the
    supplied symbol the service will return both gene and transcript (it
    supports transient links).

    Args:
        species: Species name
        name: name of gene

    Retturns:
        List of dictionaries satisying the ``species`` and ``name`` provided

    '''
    url = f'{SERVER}/xrefs/symbol/{species}/{name}?'
    return getURL(url)  # type: ignore


def getBiotypes(species: str) -> List[Dict]:
    '''List the functional classifications of gene models that Ensembl
    associates with a particular species. Useful for restricting the type of
    genes/transcripts retrieved by other endpoints.

    Args:
        species: Species name

    Retturns:
        List of dictionaries satisying the ``species`` provided

    '''
    url = f'{SERVER}/info/biotypes/{species}?'
    return getURL(url)  # type: ignore


def lookUpID(
    eid: str,
    is_protein_coding: bool = True,
) -> Dict:
    '''Find the species and database for a single identifier e.g. gene,
    transcript, protein. DOES not work for an exon ID.
    Filter out non-protein coding exons optionally

    Args:
        eid: Exon ID
        is_protein_coding: If True, look for Transcripts

    Returns:
        Dictionary given the Exon ID

    '''
    global ensembl_cache
    url = SERVER + '/lookup/id/%s?expand=1;utr=1;phenotypes=1' % (eid)
    res = getCache(url, eid, ensembl_cache)

    # is it a gene?
    if is_protein_coding and 'G' in eid:
        out = []
        for item in res['Transcript']:
            if item['biotype'] == 'protein_coding':
                out.append(item)
        res['Transcript'] = out
    return res


def lookUpIDList(
    id_list: List[str],
    species: str,
    is_protein_coding: bool = True,
) -> Dict[str, Dict]:
    '''Get the list of dictionaries for the Exon ``id_list``

    Args:
        id_list: List of Exon IDs
        species: Species name
        is_protein_coding: If True, look for Transcripts

    Returns:
        Dictionary of dictionaries given the Exon IDs

    '''
    global ensembl_cache

    url = f'{SERVER}/lookup/id/{species}'

    data: Dict = {
        'ids': id_list,
        'expand': 1,
        'utr': 1,
        'phenotypes': 1,
    }
    res = getCacheList(url, data, id_list, ensembl_cache)

    if is_protein_coding:
        for eid, lookup in res.items():
            if 'G' in eid:  # not a transcript
                out = []
                for item in lookup['Transcript']:
                    if item['biotype'] == 'protein_coding':
                        out.append(item)
                lookup['Transcript'] = out
    return res


def lookUpSymbol(
    species: str,
    symbol: str,
    is_protein_coding: bool = True,
) -> Dict:
    '''Find the species and database for a single identifier e.g. gene,
    transcript, protein. DOES not work for an exon ID.
    Filter out non-protein coding exons optionally

    Args:
        species: Species name
        symbol: Symbol to reference
        is_protein_coding: If True, look for Transcripts

    Returns:
        Dictionaries given the arguments

    '''
    global ensembl_cache
    url = f'{SERVER}/lookup/symbol/{species}/{symbol}?expand=1'
    res = getCache(url, symbol, ensembl_cache[species])

    if is_protein_coding:
        out = []
        for item in res['Transcript']:
            if item['biotype'] == 'protein_coding':
                out.append(item)
        res['Transcript'] = out
    return res


def lookUpSymbolList(
    species: str,
    symbol_list: List[str],
    is_protein_coding: bool = True,
) -> Dict[str, Any]:
    '''Find the species and database for a single identifier e.g. gene,
    transcript, protein. DOES not work for an exon ID.
    Filter out non-protein coding exons optionally

    Args:
        species: Species name
        symbol: Symbol to reference
        is_protein_coding: If True, look for Transcripts

    Returns:
        Dictionary of dictionaries given the arguments

    '''
    global ensembl_cache
    url = f'{SERVER}/lookup/symbol/{species}'
    data = {
        'symbols': symbol_list,
        'expand': 1,
        'utr': 1,
        'phenotypes': 1,
    }
    res = getCacheList(
        url,
        data,
        symbol_list,
        ensembl_cache[species],
    )

    if is_protein_coding:
        for _, lookup in res.items():
            out = []
            for item in lookup['Transcript']:
                if item['biotype'] == 'protein_coding':
                    out.append(item)
            lookup['Transcript'] = out
    return res


class TranscriptAndUTR(NamedTuple):
    transcript_id: str
    utr_id: str


def getThreePrimeUTRs(
    species: str,
    symbols: List[str],
) -> List[TranscriptAndUTR]:
    '''
    Args:
        species: Species name
        symbol: Symbol to reference

    Returns:
        List of Tuples of the form::

            <canonical transcript ID>, <three prime UTR ID>

    Raises:
        ValueError: Could not find canonical transcript for symbol

    '''

    res: Dict = lookUpSymbolList(species, symbols)
    out: List[TranscriptAndUTR] = []
    for symbol in symbols:
        try:
            item_dict = res[symbol]
        except KeyError:
            print(res)
            raise
        found_utr: bool = False
        for transcript in item_dict['Transcript']:
            if transcript['is_canonical'] == 1:
                three_prime_utr_id: str = transcript['Exon'][-1]['id']
                out.append(
                    TranscriptAndUTR(
                        transcript['id'],
                        three_prime_utr_id,
                    ),
                )
                found_utr = True
        if not found_utr:
            raise ValueError(
                f'Could not find canonical transcript for {symbol}',
            )
    return out


def convertCDNA2Genome(
        transcript_id: str,
        idxs: Tuple[int, int],
) -> Dict:
    '''Indices are relevant to the start of the transript
    i.e. idx 1 is the first base of the transcipt and would be 15881264
    for CALB1

    Args:
        transcript_id: Transcript ID
        idxs:  Indexes to search between

    Returns:
        Dictionary given the arguments

    '''
    url = f'{SERVER}/map/cdna/{transcript_id}/{idxs[0]}..{idxs[1]}?'
    res = getURL(url)
    return res


def convertCDS2Genome(
        transcript_id: str,
        idxs: Tuple[int, int],
) -> Dict:
    '''Indices are relevant to the start of the transript
    i.e. idx 1 is the first base of the transcipt and would be 15881264
    for CALB1

    Args:
        transcript_id: Transcript ID
        idxs:  Indexes to search between

    Returns:
        Dictionary given the arguments

    '''
    url = f'{SERVER}/map/cds/{transcript_id}/{idxs[0]}..{idxs[1]}?'
    res = getURL(url)
    return res


def getSequence(
        eid: str,
        seq_type: str = 'cdna',
) -> str:
    '''
    Args:
        eid: An Ensembl stable ID
        seq_type: Enum(genomic,cds,cdna,protein), default is cdna

    Returns:
        Sequence string matching Ensembl stable ID and type

    '''
    global ensembl_cache
    query = f'/sequence/id/{eid}?;type={seq_type}'
    url = f'{SERVER}{query}'
    return getCache(
        url,
        query,
        ensembl_cache,
        content_type='text/plain',
    )


def getRegionSequence(
        species: str,
        chromosome: str,
        start_idx: int,
        end_idx: int,
        strand: Optional[int] = None,
        is_rev: Optional[bool] = None,
) -> str:
    '''
    Args:
        species: Species name
        chromosome: Chromosome name
        start_idx: Start index
        end_idx: End index
        strand: -1 or 1 strand
        is_rev: Optional selection of strand using a boolean

    Returns:
        sequence string matching query

    Raises:
        ValueError: strand and is_rev arguments cannot both be None
        ValueError: strand argument needs to be -1 or 1

    '''
    global ensembl_cache
    if strand is None:
        if is_rev is None:
            raise ValueError(
                'strand and is_rev arguments cannot both be None',
            )
        strand = -1 if is_rev else 1
    elif strand not in (-1, 1):
        raise ValueError('strand argument needs to be -1 or 1')

    region = '%s:%d..%d:%d' % (chromosome, start_idx, end_idx, strand)
    query = f'/sequence/region/{species}/{region}'
    url = f'{SERVER}{query}'

    res = getCache(
        url,
        query,
        ensembl_cache,
        content_type='text/plain',
    )
    assert (len(res) == (end_idx - start_idx + 1))
    return res


def filterRegionSequence(
        query_seq: str,
        query_strand: int,
        transcript_id: str,
        transcript: Optional[Dict] = None,
        reference_seq: Optional[str] = None,
) -> Tuple[str, bool]:
    '''Confirm sequence exists in the transcript and return the aligned to the
    strand direction of the transcript sequence.  NOTE: Sometimes there are
    errors in probes so be sure to validate all sequence lookups!!!

    Args:
        query_seq: Query sequence
        query_strand: Query strand -1 or 1
        transcript_id: Transcript ID
        transcript_dict: Default is None.  If provided omit `lookUp` call
        reference_seq: Default is None.  If provided omit `getSequence` call

    Returns:
        Tuple of the form::

            sequence, was_rc

        sequence corresponding to the query.  If transcript_id is provided
        the sequence will exist in the transcript and get reverse complemented
        as necessary and was_rc should be checked

    Raises:
        ValueError on sequence not found in the target reference sequence

    '''
    was_rc: bool = False
    query_seq_out: str = query_seq
    if transcript is None:
        transcript = lookUpID(transcript_id)
    if reference_seq is None:
        reference_seq = getSequence(transcript_id)
    if transcript['strand'] != query_strand:
        query_seq_out = reverseComplement(query_seq)
        was_rc = True
    if query_seq_out not in reference_seq:
        raise ValueError(
            f'Region sequence not in transcript_id: {transcript_id}: '
            f'{query_strand}, rc: {was_rc}',
        )
    return query_seq_out, was_rc


def overlap(
        eid: str,
        feature: str = 'variation',
) -> List[Dict[str, Any]]:
    '''
    Args:
        eid: An Ensembl stable ID
        feature: Enum(band, gene, transcript, cds, exon, repeat, simple, misc,
            variation, somatic_variation, structural_variation,
            somatic_structural_variation, constrained, regulatory, motif,
            chipseq, array_probe)

    Returns:
        Overlap dictionary given the feature

    '''
    url = f'{SERVER}/overlap/id/{eid}?;feature={feature}'
    return getCache(url, eid + feature, ensembl_cache)


def excludeVariantsRegion(region: Base) -> List[Tuple[int, int]]:
    '''
    Args:
        region: the region (e.g. an Exon) Ensembl ID

    Returns:
        list of non-variant slices of the region of the transcript.  the second
        index in the tuple is non-inclusive so (1, 10) means every index from 1
        to 9 is included.  Good for python slicing

    '''
    res = overlap(
        region.id,
        feature='variation',
    )
    next_idx, last = region.idxs
    out = []
    for item in res:
        start = item['start']
        end = item['end']
        if start != next_idx:
            out.append((next_idx, start))
        next_idx = end + 1
    if next_idx != last:
        out.append((next_idx, last))
    return out


def excludeVariantsAllRegions(
        transcript: Transcript,
) -> List[List[Tuple[int, int]]]:
    '''Exons are sorted by index in a transcript
    must look only in exons
    could get all regions at once by scanning transcript but then I would need
    to write a parser.  Might as well just deal with the hit.
    Could write an async version to speed it up

    Args:
        transcript: :class:`Transcript` to perfeorm exclusion

    Returns:
        List of lists of region tuple pairs

    '''
    regions = []
    for item in transcript.Exon:
        regions.append(excludeVariantsRegion(Exon(item)))
    return regions


def idxLo2Hi(item: Dict) -> Tuple[int, int]:
    '''Sometimes start and end are reversed see rs214083637 in
    ENSMUSE00000339485 whereby start is 15881352 and end is 15881351
    despite strand being 1

    Args:
        item: Item dictionary

    Returns:
        Tuple for the item of the form:

            <start index>, <end index>

    '''
    start, end = item['start'], item['end']
    if end < start:
        start, end = end, start
    return start, end


def excludeVariantsSet(
        region: Base,
) -> Set[int]:
    '''
    Args:
        region: the region (e.g. an Exon) Ensembl ID

    Returns:
        List of non-variant slices of the region of the transcript.  the second
        index in the tuple is non-inclusive so (1, 10) means every index from 1
        to 9 is included.  Good for python slicing

    '''
    res = overlap(
        region.id,
        feature='variation',
    )
    out = set()
    for item in res:
        start, end = idxLo2Hi(item)
        out |= set(range(start, end + 1))
    return out


def permittedRegions(
        transcript: Transcript,
) -> List[List[Tuple[int, int]]]:
    '''
    Args:
        transcript:  :class:`Transcript` as calculated through a call to
            :meth:`lookUpID` and wrapped in the class

    Returns:
        List of Lists of permitted non-variant regions per exon

    '''
    # capture a set to make sure we exclude all the right indices
    excluded_indices_set = excludeVariantsSet(transcript)
    out = []
    for item in transcript.Exon:
        # assume Exon indices go from low to high and are NOT reversed
        idx_l, idx_h = idxLo2Hi(item)
        exon_set = set(range(idx_l, idx_h + 1))

        delta_set = exon_set.difference(excluded_indices_set)
        indices = sorted(list(delta_set))
        i = indices[0]
        k = i
        slice_list = []
        for j in indices[1:]:
            if j > k + 1:
                slice_list.append((i, k + 1))  # k + 1 because python slicing
                i = j
            k = j
        slice_list.append((i, k + 1))
        out.append(slice_list)
    return out


def slicedSequence(
    exon_id: str,
    offset: int,
    is_rev: bool,
    regions: List[Tuple[int, int]],
) -> List[Tuple[int, str]]:
    '''For reverse strands we must rc the strand.  All indices store are in
    terms of the forward strand and reversing things only matters when dealing
    with sequence.  Exons are listed in reverse order by index for reverse
    strands but it doesn't matter since we deal with exons one at a time.
    Exon indices seem to always be::

         ``start: low index, end: high index``

    but variants don't seem to always enforce that.  see :meth:`idxLo2Hi` for
    an example of this

    Args:
        exon_id:  Ensembl ID of the exon
        offset: the index to subtract to get the slicing correct.  This
            should be the starting index of the exon itself
        is_rev: whether this needs to be dealt with as a reverse strand
        regions:  tuples of slices to partion the sequence string

    Returns:
        the list of genomic start index on the forward strand and
            partioned sequence string of the form::

            idx, sequence

            Keep in mind for the reverse strand that this index is the end
            of the the mRNA sequence, and for the forward strand it is the
            beginning

    '''
    seq_total = getSequence(exon_id)
    if is_rev:
        seq_total = seq_total[::-1]
    out = []
    # base_count = 0
    for region in regions:
        # print(region[0] - offset, region[1] - offset)
        segment = seq_total[
            region[0] - offset:
            region[1] - offset
        ]
        # base_count += len(segment)
        if is_rev:  # reverse the segment back
            segment = segment[::-1]
        out.append((region[0], segment))
    if is_rev:
        return out[::-1]
    else:
        return out


def getArrayProbes(eid: str) -> List[Dict]:
    '''
    Args:
        eid:  ensembl ID

    Returns:
        Uniquified list of probes from a call to the `overlap` REST request

    '''
    res = overlap(
        eid,
        feature='array_probe',
    )
    return _uniqueProbes(res)


def _probe_dict_2_str(probe: Dict) -> str:
    '''
    Args:
        probe: probe dictionary to create string from

    Returns:
        string form of a probe dictionary given ``PROBE_KEYS``

    '''
    global PROBE_KEYS
    return ''.join(str(probe[x]) for x in PROBE_KEYS)


def _uniqueProbes(probe_list: List[Dict]) -> List[Dict]:
    '''
    Args:
        probe_list: List of probe dictionaries to create string for

    Returns:
        List of dictionaries string form of a probe dictionary given
        ``PROBE_KEYS``

    '''
    probe_cache: Set[str] = set()
    out: List[Dict] = []
    for probe in probe_list:
        probe_hash: str = _probe_dict_2_str(probe)
        if probe_hash not in probe_cache:
            out.append(probe)
            probe_cache.add(probe_hash)
    return out


def probeListGroupByProbeName(
        probe_list: List[Dict],
) -> List[Dict]:
    '''Get a dict per probe with a List of microarrays that the probe is in

    Args:
        probe_list: result from a call to :func:`getArrayProbes`

    Returns:
        Lists of probes with only one entry per `probe_name` key

    '''
    probe_name_idx_dict: Dict[str, int] = {}
    idx: int = 0
    out: List[Dict] = []
    for probe in probe_list:
        probe_name: str = probe['probe_name']
        if probe_name not in probe_name_idx_dict:
            probe_copy = copy.copy(probe)
            probe_copy['microarrays'] = [probe_copy['microarray']]
            del probe_copy['microarray']
            out.append(probe_copy)

            # Store the lookup
            probe_name_idx_dict[probe_name] = idx
            idx += 1
        else:
            j: int = probe_name_idx_dict[probe_name]
            ref_probe: Dict = out[j]
            ref_probe['microarrays'].append(probe['microarray'])
    return out


def getProbeFromList(
        probe_name: str,
        probe_list: List[Dict],
) -> Dict:
    '''Look up a probe dictionary in a `probe_list` which was a result from
    a call to  :func:`getArrayProbes`

    Args:
        probe_name: Probe name
        probe_list: Probe list to look up item in

    Returns:
        Dictionary of the probe

    Raises:
        ValueError: probe_name not found in list

    '''
    for x in probe_list:
        if x['probe_name'] == probe_name:
            return x
    raise ValueError(f'probe: {probe_name} not found in list')


def getProbesForID(
        eid: str,
        keep_n: int = 0,
) -> pd.DataFrame:
    '''Fet a dataframe of overlapping probes for a given ensembl ID

    Args:
        eid: ensembl ID of the item (Transcript, Exon, etc)
        keep_n: keep the n most frequent in the arrays founds

    Returns:
        Dataframe of the probes with columns::

            [
                'probe_name',
                'start',
                'end',
                'strand',
                'probe_length',
                'seq_region_name'
            ]

    '''
    out_list: List[Dict] = getArrayProbes(eid)
    if len(out_list) == 0:
        raise ValueError

    COLUMNS: List[str] = [
        'probe_name',
        'start',
        'end',
        'probe_length',
        'feature_type',
        'probe_set',
        'seq_region_name',
        'strand',
        'microarray',
    ]
    df: pd.DataFrame = pd.DataFrame(out_list, columns=COLUMNS)
    probe_name_series: pd.Series = df['probe_name']
    count_of_probe_use: pd.Series = probe_name_series.value_counts()

    if keep_n > 0:
        largest = count_of_probe_use[:keep_n]
        # TODO: Fix the typing of this line
        filtered_probes = df[     # type: ignore
            probe_name_series.isin(
                largest.index.values,
            )
        ]
    else:
        filtered_probes = df    # type: ignore

    columns_to_keep: List[str] = [
        'probe_name',
        'start',
        'end',
        'strand',
        'probe_length',
        'seq_region_name',
    ]
    filtered_probes = filtered_probes.loc[:, columns_to_keep].drop_duplicates()

    # Filter out probes where the length doesn't match the index delta
    filtered_probes = filtered_probes[
        filtered_probes.end -
        filtered_probes.start + 1 == filtered_probes.probe_length
    ]

    # Add an array frequency column
    array_freq: List[int] = [
        count_of_probe_use.loc[filtered_probes['probe_name'].iloc[i]]
        for i in range(len(filtered_probes))
    ]
    filtered_probes = filtered_probes.assign(array_freq=array_freq)

    return filtered_probes.reset_index(drop=True)  # type: ignore


def permittedSequences(
        transcript: Transcript,
        exon_id: Optional[str] = None,
) -> List[List[Tuple[int, str]]]:
    '''Get a list of start indices and their associated Exon dictionary

    Args:
        transcript:  :class:`Transcript` as calculated through a call to
            :meth:`lookUpID` and wrapped in the class
        exon_id: optional to get a specific exon's pemitted sequences.

    Returns:
        List of list of the genomic start indices and the permitted
            sequences of the exons

    '''
    slices: List[List[Tuple[int, int]]] = permittedRegions(transcript)
    is_rev: bool = transcript.strand == -1
    if exon_id is not None:
        regions = []
        for i, item in enumerate(transcript.Exon):
            if item['id'] == exon_id:
                regions = slices[i]
                break
        if len(regions) > 0:
            idx_l, idx_h = idxLo2Hi(item)
            return [slicedSequence(exon_id, idx_l, is_rev, regions)]
        else:
            raise ValueError('exon_id {} not there'.format(exon_id))
    else:
        out = []
        for item, regions in zip(transcript.Exon, slices):
            idx_l, idx_h = idxLo2Hi(item)
            # print(idx_h - idx_l)
            out.append(slicedSequence(item['id'], idx_l, is_rev, regions))
        return out


if __name__ == '__main__':
    CALB1 = 'ENSMUSG00000028222'
    # convert to genome look up:
    # res = lookUpID(CALB1)
    # lookup = LookUp(res)
    # transcript = Transcript(lookup.transcripts[0])
    # pprint(permittedSequences(transcript, 'ENSMUSE00000339485'))
    # pprint(permittedSequences(transcript))

    # GFAP = 'ENSMUSG00000020932'
    # res = lookUpID(GFAP)
    # lookup = LookUp(res)
    # transcript = Transcript(lookup.transcripts[0])
    # print(transcript.id)
    # # pprint(permittedSequences(transcript, 'ENSMUSE00000513876'))
    # pprint(permittedSequences(transcript))
    pprint(lookUpSymbolList('mouse', ['GFAP', 'SST']))

'''
a = {'Transcript': [{'Exon': [{'assembly_name': 'GRCm38',
                           'db_type': 'core',
                           'end': 15881477,
                           'id': 'ENSMUSE00000339485',
                           'object_type': 'Exon',
                           'seq_region_name': '4',
                           'species': 'mus_musculus',
                           'start': 15881264,
                           'strand': 1,
                           'version': 2},
                          {'assembly_name': 'GRCm38',
                           'db_type': 'core',
                           'end': 15882034,
                           'id': 'ENSMUSE00001218329',
                           'object_type': 'Exon',
                           'seq_region_name': '4',
                           'species': 'mus_musculus',
                           'start': 15881958,
                           'strand': 1,
                           'version': 1}],
                 'Parent': 'ENSMUSG00000028222',
                 'Translation': {'Parent': 'ENSMUST00000029876',
                                 'db_type': 'core',
                                 'end': 15904763,
                                 'id': 'ENSMUSP00000029876',
                                 'length': 261,
                                 'object_type': 'Translation',
                                 'species': 'mus_musculus',
                                 'start': 15881399},
                 'UTR': [{'Parent': 'ENSMUST00000029876',
                          'assembly_name': 'GRCm38',
                          'db_type': 'core',
                          'end': 15881398,
                          'id': 'ENSMUST00000029876',
                          'object_type': 'five_prime_UTR',
                          'seq_region_name': '4',
                          'source': 'ensembl_havana',
                          'species': 'mus_musculus',
                          'start': 15881264,
                          'strand': 1},
                         {'Parent': 'ENSMUST00000029876',
                          'assembly_name': 'GRCm38',
                          'db_type': 'core',
                          'end': 15908064,
                          'id': 'ENSMUST00000029876',
                          'object_type': 'three_prime_UTR',
                          'seq_region_name': '4',
                          'source': 'ensembl_havana',
                          'species': 'mus_musculus',
                          'start': 15904764,
                          'strand': 1}],
                 'assembly_name': 'GRCm38',
                 'biotype': 'protein_coding',
                 'db_type': 'core',
                 'display_name': 'Calb1-201',
                 'end': 15908064,
                 'id': 'ENSMUST00000029876',
                 'is_canonical': 1,
                 'logic_name': 'ensembl_havana_transcript',
                 'object_type': 'Transcript',
                 'seq_region_name': '4',
                 'source': 'ensembl_havana',
                 'species': 'mus_musculus',
                 'start': 15881264,
                 'strand': 1,
                 'version': 1}],
 'assembly_name': 'GRCm38',
 'biotype': 'protein_coding',
 'db_type': 'core',
 'description': 'Calbindin  [Source:UniProtKB/Swiss-Prot;Acc:P12658]',
 'display_name': 'Calb1',
 'end': 15908064,
 'id': 'ENSMUSG00000028222',
 'logic_name': 'ensembl_havana_gene',
 'object_type': 'Gene',
 'seq_region_name': '4',
 'source': 'ensembl_havana',
 'species': 'mus_musculus',
 'start': 15881264,
 'strand': 1,
 'version': 2}
 '''
