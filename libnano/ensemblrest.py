# -*- coding: utf-8 -*-
import sys
import io
import os.path
import atexit
from pathlib import Path
from typing import (
    List,
    Tuple,
    Set,
    Any,
    Union
)
from pprint import pprint
import json
import pickle

import requests

USE_CACHE: bool = True
TIMEOUT_FAST: float = 2.0
TIMEOUT_SLOW: float = 6.0
TIMEOUT_REQU: float = TIMEOUT_SLOW

SERVER: str = "https://rest.ensembl.org"
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
    'version'
]
TRANSLATION_KEYS: List[str] = [
    'Parent',
    'db_type',
    'end',
    'id',
    'length',
    'object_type',
    'species',
    'start'
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
    'version'
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
    'version'
]

SPECIES_NAMES: dict = {'mouse': ['mus_musculus']}
ASSEMBLY_NAMES: dict = {'mouse': ['GRCm38']}

DB_NAMES: List[str] = ['EntrezGene', 'MGI', 'Uniprot_gn', 'WikiGene']

OBJECT_TYPES: List[str]  = [
    'Translation',
    'Exon',
    'five_prime_UTR',
    'three_prime_UTR',
    'Transcript'
]
# From overlap(x, feature='variation') call
CONSEQUENCE_TYPES: List[str] = [
    '3_prime_UTR_variant',
    'intron_variant',
    'synonymous_variant',
    'non_coding_transcript_exon_variant',
    'splice_region_variant',
    '5_prime_UTR_variant'
]
POST_JSON: dict ={ "Content-Type" : "application/json", "Accept" : "application/json"}

home_path: str = str(Path.home())
CACHE_FILE: str = os.path.join(home_path, '.ENSEMBLCACHE.pickle')

SPECIES_LIST: List[str] = ['mouse', 'human']
_cache_dirty: bool
ensembl_cache: dict
THE_FILE: str = os.path.basename(__file__)

def makeCache(species_list: List[str] = SPECIES_LIST) -> dict:
    global _cache_dirty
    global SPECIES_LIST
    d = { species: {} for species in species_list}
    d['species_list'] = species_list
    SPECIES_LIST = species_list
    _cache_dirty = False
    return d
# end def

if Path(CACHE_FILE).exists():
    try:
        with io.open(CACHE_FILE, 'rb') as fd:
            ensembl_cache = pickle.load(fd)
            SPECIES_LIST = ensembl_cache['species_list']
            _cache_dirty = False
            print("LOADED cache for {}: {}".format(THE_FILE, CACHE_FILE))
    except:
        print("Couldn't load cache for {}: {}".format(THE_FILE, CACHE_FILE))
        ensembl_cache = makeCache()
else:
    ensembl_cache = makeCache()

def clearCache(species_list: List[str] = None):
    global _cache_dirty
    global ensembl_cache
    _cache_dirty = False
    ensembl_cache = makeCache(species_list=species_list)
# end def

def addSpecies(species: str):
    global SPECIES_LIST
    if species not in SPECIES_LIST:
        SPECIES_LIST.append(species)
        ensembl_cache[species] = {}
        _cache_dirty = True
# end def

def _closeCache():
    global _cache_dirty
    global ensembl_cache
    if _cache_dirty:
        with io.open(CACHE_FILE, 'wb') as fd:
            pickle.dump(ensembl_cache, fd)
            print("UPDATED {} cache".format(THE_FILE))
atexit.register(_closeCache)

def getCache(   url: str,
                arg: str,
                cache: dict,
                content_type: str = "application/json") -> Union[dict, str]:
    global _cache_dirty
    global USE_CACHE

    res = cache.get(arg) if USE_CACHE else None

    if res is None:
        res = getURL(url, content_type=content_type)
        cache[arg] = res
        _cache_dirty = True
    return res
# end def

def getCacheList(url: str, data: dict, arg_list: List[str], cache: dict) -> dict:
    global _cache_dirty
    global USE_CACHE
    do_post: bool = False
    res: dict = {}
    if USE_CACHE:
        for x in arg_list:
            try:
                res[x] = cache[x]
            except:
                do_post = True
                break
    else:
        do_post = True
    if do_post:
        res: dict = postURL(url, data=data)
        cache.update(res)
        _cache_dirty = True
    return res
# end def

class Base:
    key_list: List[str] = []

    def __init__(self, d: dict):
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
# end class

class Exon(Base):
    key_list: List[str] = EXON_KEYS
# end class

class Transcript(Base):
    key_list: List[str] = TRANSCRIPT_KEYS
# end class

class LookUp(Base):
    key_list: List[str] = LOOKUP_KEYS

    @property
    def transcripts(self) -> List[Transcript]:
        return self.d['Transcript']
# end class


def getURL(url: str, content_type="application/json") -> dict:
    r = requests.get(url,
                    headers={ "Content-Type" : content_type},
                    timeout=TIMEOUT_REQU)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    if content_type == 'text/plain':
        return r.text
    else:
        return r.json()
# end def

def postURL(url: str, data: dict,
            headers: dict = POST_JSON) -> dict:
    data_send = json.dumps(data)
    r = requests.post(url,
                        headers=headers,
                        data=data_send,
                        timeout=TIMEOUT_REQU)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.json()
# end def

def archiveID(eid: str) -> dict:
    '''Uses the given identifier to return the archived sequence
    '''
    url = SERVER + "/archive/id/%s?" % (eid)
    return getURL(url)
# end def

def xRefName(species: str, name: str) -> List[dict]:
    '''Performs a lookup based upon the primary accession or display label of an
     external reference and returning the information we hold about the entry
    '''
    url = SERVER + "/xrefs/name/%s/%s?" % (species, name)
    return getURL(url)
# end def

def externalXRefs(species: str, name: str) -> dict:
    res: List[dict] = xRefName(species, name)
    out = {}
    for item in res:
        dbname = item['dbname']
        description = item['description']
        info_type = item['info_type']
        primary_id: item['primary_id']
        synonyms: item['synonyms']
        version = item['version']
        if out.get(dbname) is None:
            out[dbname] = []
        out[dbname].append({
            'description' : item['description'],
            'display_id': item['display_id'],
            'info_type' : item['info_type'],
            'info_text' : item['info_text'],
            'primary_id': item['primary_id'],
            'synonyms': item['synonyms'],
            'version' : item['version']
        })
    return out
# end def

def xRefSymbol(species: str, name: str) -> List[dict]:
    '''Looks up an external symbol and returns all Ensembl objects linked to
    it. This can be a display name for a gene/transcript/translation, a synonym
     or an externally linked reference. If a gene's transcript is linked to the
    supplied symbol the service will return both gene and transcript (it
    supports transient links).
    '''
    url = SERVER + "/xrefs/symbol/%s/%s?" % (species, name)
    return getURL(url)
# end def

def getBiotypes(species: str) -> dict:
    '''List the functional classifications of gene models that Ensembl
    associates with a particular species. Useful for restricting the type of
    genes/transcripts retrieved by other endpoints.'''
    url = SERVER + "/info/biotypes/%s?" % (species)
    return getURL(url)
# end def

def lookUpID(   eid: str,
                is_protein_coding: bool = True) -> dict:
    '''Find the species and database for a single identifier e.g. gene,
    transcript, protein. DOES not work for an exon ID.
    Filter out non-protein coding exons optionally
    '''
    global ensembl_cache
    url = SERVER + "/lookup/id/%s?expand=1;utr=1;phenotypes=1" % (eid)
    res = getCache(url, eid, ensembl_cache)

    if is_protein_coding and 'G' in eid:
        out = []
        for item in res['Transcript']:
            if item['biotype'] == 'protein_coding':
                out.append(item)
        res['Transcript'] = out
    return res
# end def

def lookUpIDList(id_list: List[str],
                is_protein_coding: bool = True) -> dict:
    global ensembl_cache
    url: str  = SERVER + '/lookup/id/%s' % (species)
    data: dict ={   'ids': id_list,
                    'expand': 1,
                    'utr': 1,
                    'phenotypes': 1}
    res = getCacheList(url, data, id_list, ensembl_cache)

    if is_protein_coding:
        for eid, lookup in res.items():
            if 'G' in eid: #not a transcript
                out = []
                for item in lookup['Transcript']:
                    if item['biotype'] == 'protein_coding':
                        out.append(item)
                lookup['Transcript'] = out
    return res
# end def


def lookUpSymbol(   species: str,
                    symbol: str,
                    is_protein_coding: bool = True) -> dict:
    '''Find the species and database for a single identifier e.g. gene,
    transcript, protein. DOES not work for an exon ID.
    Filter out non-protein coding exons optionally
    '''
    global ensembl_cache
    url = SERVER + "/lookup/symbol/%s/%s?expand=1" % (species, symbol)
    res = getCache(url, symbol, ensembl_cache[species])

    if is_protein_coding:
        out = []
        for item in res['Transcript']:
            if item['biotype'] == 'protein_coding':
                out.append(item)
        res['Transcript'] = out
    return res
# end def

def lookUpSymbolList(   species: str,
                        symbol_list: List[str],
                        is_protein_coding: bool = True) -> dict:
    global ensembl_cache
    url: str = SERVER + '/lookup/symbol/%s' % (species)
    data: dict = {  'symbols': symbol_list,
                    'expand': 1,
                    'utr': 1,
                    'phenotypes': 1 }
    res = getCacheList(url, data, symbol_list, ensembl_cache[species])

    if is_protein_coding:
        for symbol, lookup in res.items():
            out = []
            for item in lookup['Transcript']:
                if item['biotype'] == 'protein_coding':
                    out.append(item)
            lookup['Transcript'] = out
    return res
# end def

def convertCDNA2Genome(transcript_id:str, idxs: Tuple[int, int]) -> dict:
    '''indices are relevant to the start of the transript
    i.e. idx 1 is the first base of the transcipt and would be 15881264
    for CALB1
    '''
    url = SERVER + "/map/cdna/%s/%d..%d?" % (transcript_id, idxs[0], idxs[1])
    res = getURL(url)
    return res
# end def

def convertCDS2Genome(transcript_id: str, idxs: Tuple[int, int]) -> dict:
    '''indices are relevant to the start of the transript
    i.e. idx 1 is the first base of the transcipt and would be 15881264
    for CALB1
    '''
    url = SERVER + "/map/cdna/%s/%d..%d?" % (transcript_id, idxs[0], idxs[1])
    res = getURL(url)
    return res
# end def

def getSequence(eid: str, seq_type: str='cdna') -> str:
    '''
    Args:
        eid: An Ensembl stable ID
        seq_type: Enum(genomic,cds,cdna,protein), default is cdna
    '''
    global ensembl_cache
    url = SERVER + "/sequence/id/%s?;type=%s" % (eid, seq_type)
    return getCache(url, eid + 'sequence', ensembl_cache, content_type="text/plain")
# end def

def overlap(eid: str, feature: str = 'variation') -> dict:
    '''
    Args:
        eid: An Ensembl stable ID
        feature: Enum(band, gene, transcript, cds, exon, repeat, simple, misc,
            variation, somatic_variation, structural_variation,
            somatic_structural_variation, constrained, regulatory, motif,
            chipseq, array_probe)
    '''
    url = SERVER + "/overlap/id/%s?;feature=%s" % (eid, feature)
    return getCache(url, eid+feature, ensembl_cache)
# end def

def excludeVariantsRegion(region: Base) -> List[Tuple[int,int]]:
    '''
    Args:
        region: the region (e.g. an Exon) Ensembl ID

    Returns:
        list of non-variant slices of the region of the transcript.  the second
        index in the tuple is non-inclusive so (1, 10) means every index from 1
        to 9 is included.  Good for python slicing
    '''
    res = overlap(region.id, feature='variation')
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
# end def

def excludeVariantsAllRegions(transcript: Transcript) -> List[List[Tuple[int, int]]]:
    '''exons are sorted by index in a transcript
    must look only in exons
    could get all regions at once by scanning transcript but then I would need
    to write a parser.  Might as well just deal with the hit.
    Could write an async version to speed it up
    '''
    regions = []
    for item in transcript.Exon:
        regions.append(excludeVariantsRegion(Exon(item)))
    return regions
# end def

def idxLo2Hi(item):
    '''sometimes start and end are reversed see rs214083637 in
    ENSMUSE00000339485 whereby start is 15881352 and end is 15881351
    despite strand being 1
    '''
    start, end = item['start'], item['end']
    if end < start:
        start, end = end, start
    return start, end
# end def

def excludeVariantsSet(region: Base) -> Set[int]:
    '''
    Args:
        region: the region (e.g. an Exon) Ensembl ID

    Returns:
        list of non-variant slices of the region of the transcript.  the second
        index in the tuple is non-inclusive so (1, 10) means every index from 1
        to 9 is included.  Good for python slicing
    '''
    res = overlap(region.id, feature='variation')
    out = set()
    for item in res:
        start, end = idxLo2Hi(item)
        out |= set(range(start, end + 1))
    return out
# end def

def permittedRegions(transcript: Transcript) -> List[List[Tuple[int, int]]]:
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
                slice_list.append((i, k + 1)) # k + 1 because python slicing
                i = j
            k = j
        slice_list.append((i, k + 1))
        out.append(slice_list)
    return out
# end def

def slicedSequence( exon_id: str,
                    offset: int,
                    is_rev: bool,
                    regions: List[Tuple[int, int]]) -> List[Tuple[int, str]]:
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
        segment = seq_total[region[0] - offset:region[1] - offset]
        # base_count += len(segment)
        if is_rev: # reverse the segment back
            segment = segment[::-1]
        out.append( (region[0], segment) )
    if is_rev:
        return out[::-1]
    else:
        return out
# end def

def permittedSequences( transcript: Transcript,
                        exon_id: str = None) -> List[List[Tuple[int, str]]]:
    '''
    Args:
        transcript:  :class:`Transcript` as calculated through a call to
            :meth:`lookUpID` and wrapped in the class
        exon_id: optional to get a specific exon's pemitted sequences.

    Returns:
        List of list of the genomic start indices and the permitted
            sequences of the exons
    '''
    slices: List[List[Tuple[int, int]]] = permittedRegions(transcript)
    is_rev = transcript.strand == -1
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
            raise ValueError("exon_id {} not there".format(exon_id))
    else:
        out = []
        for item, regions in zip(transcript.Exon, slices):
            idx_l, idx_h = idxLo2Hi(item)
            out.append(slicedSequence(item['id'], idx_l, is_rev, regions))
        return out
# end def

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