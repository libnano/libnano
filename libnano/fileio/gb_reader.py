# -*- coding: utf-8 -*-
"""
http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
http://www.insdc.org/documents/feature_table.html

All keys are native strings, as are values, except the origin which is always
a python 3 byte string (not unicode)
"""

import re
import io
from collections import OrderedDict
from typing import List

def parse(filepath: str, is_ordered: bool = False) -> dict:
    """
    Args:
        filepath:
        is_ordered: default False. if True will retain the order of the qualifiers
    """
    d = {'info': {}}
    d_info = d['info']
    with io.open(filepath, 'r', encoding='utf-8') as fd:
        raw = fd.read()
    start, _, origin = raw.partition("ORIGIN")
    start, _, features = start.partition("FEATURES             Location/Qualifiers\n")

    d_info.update(parseLocus(start))
    d_info.update(parseDefinition(start))
    d_info.update(parseAccession(start))
    d_info.update(parseVersion(start))
    d_info.update(parseDBLink(start))
    d_info.update(parseKeywords(start))
    d_info.update(parseSource(start))
    d_info.update(parseOrganism(start))

    d_info['references'] = parseReference(start)
    _, _, comment = start.partition("COMMENT     ")
    parseComment(d_info, comment)
    d['features'] = parseFeatures(features, is_ordered)
    d['seq'] = parseOrigin(origin)
    return d
# end def

def parseComment(d: dict, comment: str):
    if comment != '':
        # get rid of ApE empty comment
        if comment.startswith("\nCOMMENT     "):
            comment = comment[13:]
        idx_genome_asm_data = -1
        genome_asm_data_newline_count = 0
        lines = comment.split('\n')
        lines = [line.strip() for line in lines]
        # print(lines)
        # handle ##Genome-Assembly-Data-START## edge case
        for i, line in enumerate(lines):
            if line == '':
                genome_asm_data_newline_count += 1
            elif genome_asm_data_newline_count == 2:
                idx_genome_asm_data = i
                genome_asm_data_newline_count = 0
            else:
                genome_asm_data_newline_count = 0
        # end for
        if idx_genome_asm_data < 0:
            d['comment'] = " ".join(lines)
        else:
            d['comment'] = [" ".join(lines[:idx_genome_asm_data-2]),
                    lines[idx_genome_asm_data:-1]]
# end def

re_locus: List[str] = [
    "^LOCUS",                                   # field
    " +(?P<name>[\w|.]+)",                      # name
    " +(?P<length>[0-9]+) bp",                  # sequence length
    "(?: +(?P<stranded>[a-z]{2})-)?",           # opt: ss, ds, ms
    " *(?P<molecule_type>[a-z|A-Z|-|]{2,6})",   # molecule type
    " +(?P<form>[\w]{6,8})?",                   # linear or circular
    " +(?P<gb_division>[a-z|A-Z]{3})?",         # Genbank division
    " +(?P<mod_date>[0-9]+-[A-Z]+-[0-9]+)",     # modification date
    ".*\n"                                      # match line end
]
RE_LOCUS: str = "".join(re_locus)

re_definition: List[str] = [
    "^DEFINITION",                          # field
    " +(?P<definition>(?:.*\n)(?: .*\n)*)"  # look ahead assertion for multiline
]
RE_DEFINITION: str = "".join(re_definition)

re_accession: List[str] = [
                    "^ACCESSION",                           # field
                    " +(?P<accession>[\w|.]*)"   # look ahead assertion for multiline
                    ".*\n"                                  # match line end
]
RE_ACCESSION: str = "".join(re_accession)

re_version: List[str] = [
    "^VERSION",                     # field
    " +(?P<version>[\w|.]+)",    # version
    " +GI:(?P<GI>[\w|.]+)"       # gi field
    ".*\n"                       # match line end
]

RE_DBLINK: str = "^DBLINK +(?P<dblink>[\w|:| |.]+)\n"

RE_VERSION: str = "".join(re_version)

re_keywords: List[str] = [
    "^KEYWORDS",
    " +(?P<keywords>[\w|.]*)"
    ".*\n"
]
RE_KEYWORDS= "".join(re_keywords)

re_source: List[str] = [
    "^SOURCE",
    " +(?P<source>.*)",
    "\n"
]
RE_SOURCE: str = "".join(re_source)

re_organism: List[str] =  [
    "^  ORGANISM",                          # field
    "(?: +(?P<organism0>(?:.*\n))?",
    "(?: +(?P<organism1>(?:.*\n)(?: .*\n)*))?)"  # multiline
]
RE_ORGANISM: str = "".join(re_organism)

RE_PREAMBLE: str = (RE_LOCUS +
                    RE_DEFINITION +
                    RE_ACCESSION +
                    RE_VERSION +
                    RE_DBLINK +
                    RE_KEYWORDS +
                    RE_SOURCE +
                    RE_ORGANISM
)
RE_COMP_PREABLE: '_sre.SRE_Pattern' = re.compile(RE_PREAMBLE, flags=re.M)

RE_COMP_LOCUS: '_sre.SRE_Pattern' = re.compile(RE_LOCUS, flags=re.M)

def parseLocus(raw: str, isbytes: str = False) -> dict:
    m = re.match(RE_COMP_LOCUS, raw)
    d = m.groupdict()
    d['length'] = int(d['length'])
    return d
#end def

RE_COMP_DEFINITION = re.compile(RE_DEFINITION, flags=re.M)
def parseDefinition(raw: str, isbytes: bool = False) -> dict:
    m = re.search(RE_COMP_DEFINITION, raw)
    if m is None:
        return {'definition': None}
    d = m.groupdict()
    if d['definition'] is not None:
        temp_l = d['definition'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        d['definition'] = " ".join(temp_l)[:-1]
    return d
#end def

RE_COMP_ACCESSION = re.compile(RE_ACCESSION, flags=re.M)
def parseAccession(raw: str) -> dict:
    m = re.search(RE_COMP_ACCESSION, raw)
    if m is None:
        return {'accession': None}
    d = m.groupdict()
    return d
# end def

RE_COMP_VERSION = re.compile(RE_VERSION, flags=re.M)
def parseVersion(raw: str) -> dict:
    m = re.search(RE_COMP_VERSION, raw)
    if m is None:
        # print("Version none")
        return {'version': None}
    d = m.groupdict()
    # print("version", d)
    return d
# end def

RE_COMP_DBLINK = re.compile(RE_DBLINK, flags=re.M)
def parseDBLink(raw: str) -> dict:
    m = re.search(RE_COMP_DBLINK, raw)
    if m is None:
        return {'dblink': None}
    d = m.groupdict()
    return d
# end def

RE_COMP_KEYWORDS = re.compile(RE_KEYWORDS, flags=re.M)
def parseKeywords(raw: str) -> dict:
    m = re.search(RE_COMP_KEYWORDS, raw)
    if m is None:
        return {'keywords': None}
    d = m.groupdict()
    return d
# end def

RE_COMP_SOURCE = re.compile(RE_SOURCE, flags=re.M)
def parseSource(raw: str) -> dict:
    m = re.search(RE_COMP_SOURCE, raw)
    if m is None:
        return {'source': None}
    d = m.groupdict()
    return d
# end def

RE_COMP_ORGANISM: '_sre.SRE_Pattern' = re.compile(RE_ORGANISM, flags=re.M)
def parseOrganism(raw: str) -> dict:
    m = re.search(RE_COMP_ORGANISM, raw)
    if m is None:
        return {'organism': [None, None]}
    d = m.groupdict()

    temp_l = d['organism0'].split('\n')
    temp_l = [x.strip() for x in temp_l]
    org0 = " ".join(temp_l)[:-1]

    org1 = None
    if d['organism1'] is not None:
        temp_l = d['organism1'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        org1 = " ".join(temp_l)[:-1]

    del d['organism0']
    del d['organism1']

    d['organism'] = [org0, org1]
    return d
# end def

"""
REFERENCE   1  (bases 1 to 5028)
  AUTHORS   Torpey,L.E., Gibbs,P.E., Nelson,J. and Lawrence,C.W.
  TITLE     Cloning and sequence of REV7, a gene whose function is required for
            DNA damage-induced mutagenesis in Saccharomyces cerevisiae
  JOURNAL   Yeast 10 (11), 1503-1509 (1994)
  PUBMED    7871890
"""
re_reference: List[str] = [
    "^REFERENCE",
    " +(?P<r_index>[0-9]+)(?: +\(bases (?P<start_idx>[0-9]+) to (?P<end_idx>[0-9]+)\)){0,1}",
    ".*\n",
    "^  AUTHORS",
    " +(?P<authors>.+)\n",
    "^  TITLE",                             # field
    " +(?P<title>(?:.*\n)(?: .*\n)*)",      # multiline
    "^  JOURNAL",
    " +(?P<journal_info>.+\n(?: {12}.+\n)*)",
    "(?:^  PUBMED +(?P<pubmed>[0-9]+)\n){0,1}"
]
RE_REFERENCE: str = "".join(re_reference)
RE_COMP_REF: '_sre.SRE_Pattern' = re.compile(RE_REFERENCE, flags=re.M)

def parseReference(raw: str) -> List[dict]:
    ref_list = []

    for m in re.finditer(RE_COMP_REF, raw):
        d = m.groupdict()
        temp_l = d['title'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        d['title'] = " ".join(temp_l)[:-1]

        temp_l = d['journal_info'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        d['journal_info'] = " ".join(temp_l)[:-1]

        d['r_index'] = int(d['r_index'])
        if d['start_idx'] is not None:
            d['start_idx'] = int(d['start_idx'])
        if d['end_idx'] is not None:
            d['end_idx'] = int(d['end_idx'])
        ref_list.append(d)
    # end for
    return ref_list
#end def

def addMultivalue(d, key, val):
    if key in d:
        old_val = d[key]
        if isinstance(old_val, list):
            old_val.append(val)
        else:
            d[key] = [old_val, val]
    else:
        d[key] = val
# end def

"""
see section 3.4 Location
"""
re_feature: List[str] = [
    "^ {5}(?P<feature_key>[\w]+)",
    " +(?P<location>.+)\n",
    "(?P<qualifiers>(?:^ {21}.*\n)*)"
]
RE_FEATURE: str = "".join(re_feature)
RE_COMP_FEATURE = re.compile(RE_FEATURE, flags=re.M)


# Qualifers can have tags with /'s in the value so it's tough to escape them
# for now we need to split on "                     /"

def parseFeatures(raw: str, is_ordered: bool = False) -> List[dict]:
    features_list = []

    for feature_match in re.finditer(RE_COMP_FEATURE, raw):
        feature = feature_match.groupdict()
        if 'qualifiers' not in feature:
            print(feature)
            raise IOError("bad feature")

        d = {'type': feature['feature_key'],
            'location': feature['location'],
            # 'partials': (feature['partial5'], feature['partial3']),
            'qualifiers': OrderedDict() if is_ordered else {}
        }
        qs = d['qualifiers']
        # prevent splitting on </tags>
        qs_list = feature['qualifiers'].split('                     /')
        for qualifier in qs_list[1:]:
            # heal line breaks
            """
            Need to address the multi value key problem
            i.e. more than one value in the list of qualifiers
            """
            q_list = qualifier.split('=')
            key = q_list[0]
            yes_val = True
            try:
                q_list = q_list[1].split('\n')
                if q_list[-1] == '':
                    q_list.pop() # remove ending '' item
            except:
                q_list = ['']
                yes_val = False
            q_list = [x.lstrip() for x in q_list]
            is_str = True
            if key == 'translation':
                temp = "".join(q_list)
            elif key in ('codon_start', 'transl_table'):
                is_str = False
                temp = " ".join(q_list)
            else:
                temp = " ".join(q_list)
            if yes_val and temp[0] == '\"' and temp[-1] == '\"':
                value_to_add = temp[1:-1] if is_str else int(temp[1:-1])
            elif not yes_val:
                value_to_add = None
            else:
                value_to_add = temp if is_str else int(temp)
            addMultivalue(qs, key, value_to_add)
        features_list.append(d)
    # end for
    return features_list
#end def

def parseOrigin(raw: str) -> str:
    out_list = []

    all_lines = raw.split('\n')
    start = 1 if all_lines[0].strip() == '' else 0
    for line in all_lines[start:-1]:
        temp = line.split()
        out_list += temp[1:]
    seq = "".join(out_list)

    assert(seq.isalpha())
    return seq
#end def

if __name__ == "__main__":
    import os.path as opath
    path = opath.dirname(opath.dirname(opath.dirname(opath.abspath(__file__))))
    fn = opath.join(path, "tests", "test_data", "failed.gb")
    def main():
        d = parse(fn)
        return d
    # end def
    print(main())
