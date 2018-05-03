# -*- coding: utf-8 -*-
"""
http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
http://www.insdc.org/documents/feature_table.html

All keys are native strings, as are values, except the origin which is always
a python 3 byte string (not unicode)
"""

import re
import io
import sys
from collections import OrderedDict

try:
    import six
except ImportError:
    from libnano.helpers import six


NEWLINE_STR = '\r\n' if sys.platform == 'win32' else '\n'
NEWLINE_BYT = b'\r\n' if sys.platform == 'win32' else b'\n'

def parse(filepath, is_ordered=False):
    """
    is_ordered == True will retain the order of the qualifiers
    """
    _b = six.b
    d = {b'info': {}}
    d_info = d[b'info']
    # with io.open(filepath, 'r', encoding='utf-8') as fd:
    with io.open(filepath, 'rb') as fd:
        raw = fd.read()
    start, _, origin = raw.partition(b"ORIGIN")
    start, _, features = start.partition(b"FEATURES             Location/Qualifiers%s" % (NEWLINE_BYT))

    parseLocus(start, d_info)
    parseDefinition(start, d_info)
    parseAccession(start, d_info)
    parseVersion(start, d_info)
    parseDBLink(start, d_info)
    parseKeywords(start, d_info)
    parseSource(start, d_info)
    parseOrganism(start, d_info)

    d_info[b'references'] = parseReference(start)
    _, _, comment = start.partition(b"COMMENT     ")
    parseComment(d_info, comment)
    d[b'features'] = parseFeatures(features, is_ordered)
    d[b'seq'] = parseOrigin(origin)
    return d
# end def

def parseComment(d, comment):
    if comment != b'':
        # get rid of ApE empty comment
        if comment.startswith(b"%sCOMMENT     " % (NEWLINE_BYT)):
            comment = comment[13:]
        idx_genome_asm_data = -1
        genome_asm_data_newline_count = 0
        lines = comment.split(NEWLINE_BYT)
        lines = [line.strip() for line in lines]
        # print(lines)
        # handle ##Genome-Assembly-Data-START## edge case
        for i, line in enumerate(lines):
            if line == b'':
                genome_asm_data_newline_count += 1
            elif genome_asm_data_newline_count == 2:
                idx_genome_asm_data = i
                genome_asm_data_newline_count = 0
            else:
                genome_asm_data_newline_count = 0
        # end for
        if idx_genome_asm_data < 0:
            d[b'comment'] = b" ".join(lines)
        else:
            d[b'comment'] = [b" ".join(lines[:idx_genome_asm_data-2]),
                    lines[idx_genome_asm_data:-1]]
# end def

re_locus = [
                "^LOCUS",                                   # field
                " +(?P<name>[\w|.]+)",                      # name
                " +(?P<length>[0-9]+) bp",                  # sequence length
                "(?: +(?P<stranded>[a-z]{2})-)?",           # opt: ss, ds, ms
                " *(?P<molecule_type>[a-z|A-Z|-|]{2,6})",   # molecule type
                " +(?P<form>[\w]{6,8})?",                   # linear or circular
                " +(?P<gb_division>[a-z|A-Z]{3})?",         # Genbank division
                " +(?P<mod_date>[0-9]+-[A-Z]+-[0-9]+)",     # modification date
                ".*%s" % (NEWLINE_STR)                      # match line end
            ]
re_locus = six.b("".join(re_locus))

re_definition = [   "^DEFINITION",                         # field
                    " +(?P<definition>(?:.*%s)(?: .*%s)*)" % (NEWLINE_STR, NEWLINE_STR) # look ahead assertion for multiline
                ]
re_definition = six.b("".join(re_definition))

re_accession = [   "^ACCESSION",                 # field
                    " +(?P<accession>[\w|.]*)"   # look ahead assertion for multiline
                    ".*",                        # match line end
                    NEWLINE_STR
                ]
re_accession = six.b("".join(re_accession))

re_version = [   "^VERSION",                     # field
                    " +(?P<version>[\w|.]+)",    # version
                    " +GI:(?P<GI>[\w|.]+)"       # gi field
                    ".*",                        # match line end
                    NEWLINE_STR
                ]
re_version= six.b("".join(re_version))

re_dblink = b"^DBLINK +(?P<dblink>[\w|:| |.]+)" + NEWLINE_BYT

re_keywords = ["^KEYWORDS",
                " +(?P<keywords>[\w|.]*)"
                ".*",
                NEWLINE_STR
            ]
re_keywords= six.b("".join(re_keywords))

re_source = ["^SOURCE",
            " +(?P<source>.*)",
            NEWLINE_STR
]
re_source = six.b("".join(re_source))

re_organism =  [   "^  ORGANISM",                          # field
                    "(?: +(?P<organism0>(?:.*%s))?" % NEWLINE_STR,
                    "(?: +(?P<organism1>(?:.*%s)(?: .*%s)*))?)" % (NEWLINE_STR, NEWLINE_STR) # multiline
                ]
re_organism = six.b("".join(re_organism))


re_comp_locus = re.compile(re_locus, flags=re.M)
def parseLocus(raw, d_out):
    m = re.match(re_comp_locus, raw)
    d = m.groupdict()
    d['length'] = int(d['length'])
    _b = six.b
    for k, v in d.items():
        d_out[_b(k)] = v
#end def

re_comp_definition = re.compile(re_definition, flags=re.M)
def parseDefinition(raw, d_out):
    m = re.search(re_comp_definition, raw)
    if m is None:
        d_out[b'definition'] = None
    else:
        d = m.groupdict()
        if d['definition'] is not None:
            temp_l = d['definition'].split(NEWLINE_BYT)
            temp_l = [x.strip() for x in temp_l]
            d_out[b'definition'] = b" ".join(temp_l)[:-1]
        else:
            d_out[b'definition'] = None
#end def

re_comp_accession = re.compile(re_accession, flags=re.M)
def parseAccession(raw, d_out):
    m = re.search(re_comp_accession, raw)
    if m is None:
        d_out[b'accession'] = None
    else:
        d = m.groupdict()
        d_out[b'accession'] = d['accession']
# end def

re_comp_version = re.compile(re_version, flags=re.M)
def parseVersion(raw, d_out):
    m = re.search(re_comp_version, raw)
    if m is None:
        d_out[b'version'] = None
    else:
        d = m.groupdict()
        d_out[b'version'] = d['version']
        d_out[b'GI'] = d['GI']
# end def

re_comp_dblink = re.compile(re_dblink, flags=re.M)
def parseDBLink(raw, d_out):
    m = re.search(re_comp_dblink, raw)
    if m is None:
        d_out[b'dblink'] = None
    else:
        d = m.groupdict()
        d_out[b'dblink'] = d['dblink']
# end def

re_comp_keywords = re.compile(re_keywords, flags=re.M)
def parseKeywords(raw, d_out):
    m = re.search(re_comp_keywords, raw)
    if m is None:
        d_out[b'keywords'] = None
    else:
        d = m.groupdict()
        d_out[b'keywords'] = d['keywords']
# end def

re_comp_source = re.compile(re_source, flags=re.M)
def parseSource(raw, d_out):
    m = re.search(re_comp_source, raw)
    if m is None:
        d_out[b'source'] = None
    else:
        d = m.groupdict()
        d_out[b'source'] = d['source']
# end def

re_comp_organism = re.compile(re_organism, flags=re.M)
def parseOrganism(raw, d_out):
    m = re.search(re_comp_organism, raw)
    if m is None:
        d_out[b'organism'] = [None, None]
    else:
        d = m.groupdict()

        temp_l = d['organism0'].split(NEWLINE_BYT)
        temp_l = [x.strip() for x in temp_l]
        org0 = b" ".join(temp_l)[:-1]

        org1 = None
        if d['organism1'] is not None:
            temp_l = d['organism1'].split(NEWLINE_BYT)
            temp_l = [x.strip() for x in temp_l]
            org1 = b" ".join(temp_l)[:-1]

        d_out[b'organism'] = [org0, org1]
# end def

"""
REFERENCE   1  (bases 1 to 5028)
  AUTHORS   Torpey,L.E., Gibbs,P.E., Nelson,J. and Lawrence,C.W.
  TITLE     Cloning and sequence of REV7, a gene whose function is required for
            DNA damage-induced mutagenesis in Saccharomyces cerevisiae
  JOURNAL   Yeast 10 (11), 1503-1509 (1994)
  PUBMED    7871890
"""
re_reference = [    "^REFERENCE",
                    " +(?P<r_index>[0-9]+)(?: +\(bases (?P<start_idx>[0-9]+) to (?P<end_idx>[0-9]+)\)){0,1}",
                    ".*",
                    NEWLINE_STR,
                    "^  AUTHORS",
                    " +(?P<authors>.+)",
                    NEWLINE_STR,
                    "^  TITLE",                             # field
                    " +(?P<title>(?:.*%s)(?: .*%s)*)" % (NEWLINE_STR, NEWLINE_STR),      # multiline
                    "^  JOURNAL",
                    " +(?P<journal_info>.+%s(?: {12}.+%s)*)" % (NEWLINE_STR, NEWLINE_STR),
                    "(?:^  PUBMED +(?P<pubmed>[0-9]+)%s){0,1}" % (NEWLINE_STR)
]
re_reference = six.b("".join(re_reference))
re_comp_ref = re.compile(re_reference, flags=re.M)


def parseReference(raw):
    ref_list = []

    for m in re.finditer(re_comp_ref, raw):
        d_temp = {}
        d = m.groupdict()
        temp_l = d['title'].split(NEWLINE_BYT)
        temp_l = [x.strip() for x in temp_l]
        d_temp[b'title'] = b" ".join(temp_l)[:-1]

        temp_l = d['journal_info'].split(NEWLINE_BYT)
        temp_l = [x.strip() for x in temp_l]
        d_temp[b'journal_info'] = b" ".join(temp_l)[:-1]

        d_temp[b'r_index'] = int(d['r_index'])
        if d['start_idx'] is not None:
            d_temp[b'start_idx'] = int(d['start_idx'])
        else:
            d_temp[b'start_idx'] = None
        if d['end_idx'] is not None:
            d_temp[b'end_idx'] = int(d['end_idx'])
        else:
            d_temp[b'end_idx'] = None
        d_temp[b'authors'] = d['authors']
        d_temp[b'pubmed'] = d['pubmed']
        ref_list.append(d_temp)
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
re_feature = [  "^ {5}(?P<feature_key>[\w]+)",
                " +(?P<location>.+)",
                NEWLINE_STR,
                "(?P<qualifiers>(?:^ {21}.*%s)*)" % (NEWLINE_STR)
]
re_feature = six.b("".join(re_feature))
re_comp_feature = re.compile(re_feature, flags=re.M)


# Qualifers can have tags with /'s in the value so it's tough to escape them
# for now we need to split on "                     /"

QUOTE_BYTE = b'\"'[0]
def parseFeatures(raw, is_ordered=False):
    features_list = []

    for feature_match in re.finditer(re_comp_feature, raw):
        feature = feature_match.groupdict()
        if 'qualifiers' not in feature:
            print(feature)
            raise IOError("bad feature")

        d = {b'type': feature['feature_key'],
            b'location': feature['location'],
            # 'partials': (feature['partial5'], feature['partial3']),
            b'qualifiers': OrderedDict() if is_ordered else {}
        }
        qs = d[b'qualifiers']
        # prevent splitting on </tags>
        qs_list = feature['qualifiers'].split(b'                     /')
        for qualifier in qs_list[1:]:
            # heal line breaks
            """
            Need to address the multi value key problem
            i.e. more than one value in the list of qualifiers
            """
            q_list = qualifier.split(b'=')
            key = q_list[0]
            yes_val = True
            try:
                q_list = q_list[1].split(NEWLINE_BYT)
                if q_list[-1] == b'':
                    q_list.pop() # remove ending '' item
            except:
                q_list = [b'']
                yes_val = False
            q_list = [x.lstrip() for x in q_list]
            is_str = True
            if key == b'translation':
                temp = b"".join(q_list)
            elif key in (b'codon_start', b'transl_table'):
                is_str = False
                temp = b" ".join(q_list)
            else:
                temp = b" ".join(q_list)
            if yes_val and temp[0] == QUOTE_BYTE and temp[-1] == QUOTE_BYTE:
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

def parseOrigin(raw):
    out_list = []

    all_lines = raw.split(NEWLINE_BYT)
    start = 1 if all_lines[0].strip() == b'' else 0
    for line in all_lines[start:-1]:
        temp = line.split()
        out_list += temp[1:]
    seq = b"".join(out_list)

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
