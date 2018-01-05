import re
import io
from collections import OrderedDict

"""
http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
http://www.insdc.org/documents/feature_table.html

All keys are native strings, as are values, except the origin which is always
a python2/3 byte string (not unicode)
"""

def parse(filepath, is_ordered=False):
    """
    is_ordered == True will retain the order of the qualifiers
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

def parseComment(d, comment):
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

re_locus = [
                "^LOCUS",                                   # field
                " +(?P<name>[\w|.]+)",                      # name
                " +(?P<length>[0-9]+) bp",                  # sequence length
                "(?: +(?P<stranded>[a-z]{2})-)?",           # opt: ss, ds, ms
                " *(?P<molecule_type>[a-z|A-Z|-|]{2,6})",     # molecule type
                " +(?P<form>[\w]{6,8})?",               # linear or circular
                " +(?P<gb_division>[a-z|A-Z]{3})?",          # Genbank division
                " +(?P<mod_date>[0-9]+-[A-Z]+-[0-9]+)",     # modification date
                ".*\n"                                      # match line end
            ]
re_locus = "".join(re_locus)

re_definition = [   "^DEFINITION",                         # field
                    " +(?P<definition>(?:.*\n)(?: .*\n)*)"  # look ahead assertion for multiline
                ]
re_definition = "".join(re_definition)

re_accession = [   "^ACCESSION",                           # field
                    " +(?P<accession>[\w|.]*)"   # look ahead assertion for multiline
                    ".*\n"                                  # match line end
                ]
re_accession = "".join(re_accession)

re_version = [   "^VERSION",                     # field
                    " +(?P<version>[\w|.]+)",    # version
                    " +GI:(?P<GI>[\w|.]+)"       # gi field
                    ".*\n"                       # match line end
                ]

re_dblink = "^DBLINK +(?P<dblink>[\w|:| |.]+)\n"
# re_dblink = "^DBLINK +(?P<dblink>[\w|\:| |.]+).*\n"

re_version= "".join(re_version)

re_keywords = ["^KEYWORDS",
                " +(?P<keywords>[\w|.]*)"
                ".*\n"
            ]
re_keywords= "".join(re_keywords)

re_source = ["^SOURCE",
            " +(?P<source>.*)",
            "\n"
]
re_source = "".join(re_source)

re_organism =  [   "^  ORGANISM",                          # field
                    "(?: +(?P<organism0>(?:.*\n))?",
                    "(?: +(?P<organism1>(?:.*\n)(?: .*\n)*))?)"  # multiline
                ]
re_organism = "".join(re_organism)

re_preamble = re_locus + re_definition + \
            re_accession + re_version + re_dblink +\
             re_keywords + re_source + re_organism

re_comp_preable = re.compile(re_preamble, flags=re.M)

re_comp_locus = re.compile(re_locus, flags=re.M)
def parseLocus(raw, isbytes=False):
    m = re.match(re_comp_locus, raw)
    d = m.groupdict()
    d['length'] = int(d['length'])
    return d
#end def

re_comp_definition = re.compile(re_definition, flags=re.M)
def parseDefinition(raw, isbytes=False):
    m = re.search(re_comp_definition, raw)
    if m is None:
        return {'definition': None}
    d = m.groupdict()
    if d['definition'] is not None:
        temp_l = d['definition'].split('\n')
        temp_l = [x.strip() for x in temp_l]
        d['definition'] = " ".join(temp_l)[:-1]
    return d
#end def

re_comp_accession = re.compile(re_accession, flags=re.M)
def parseAccession(raw):
    m = re.search(re_comp_accession, raw)
    if m is None:
        return {'accession': None}
    d = m.groupdict()
    return d
# end def

re_comp_version = re.compile(re_version, flags=re.M)
def parseVersion(raw):
    m = re.search(re_comp_version, raw)
    if m is None:
        # print("Version none")
        return {'version': None}
    d = m.groupdict()
    # print("version", d)
    return d
# end def

re_comp_dblink = re.compile(re_dblink, flags=re.M)
def parseDBLink(raw):
    m = re.search(re_comp_dblink, raw)
    if m is None:
        return {'dblink': None}
    d = m.groupdict()
    return d
# end def

re_comp_keywords = re.compile(re_keywords, flags=re.M)
def parseKeywords(raw):
    m = re.search(re_comp_keywords, raw)
    if m is None:
        return {'keywords': None}
    d = m.groupdict()
    return d
# end def

re_comp_source = re.compile(re_source, flags=re.M)
def parseSource(raw):
    m = re.search(re_comp_source, raw)
    if m is None:
        return {'source': None}
    d = m.groupdict()
    return d
# end def

re_comp_organism = re.compile(re_organism, flags=re.M)
def parseOrganism(raw):
    m = re.search(re_comp_organism, raw)
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
re_reference = [    "^REFERENCE",
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
re_reference = "".join(re_reference)
re_comp_ref = re.compile(re_reference, flags=re.M)


def parseReference(raw):
    ref_list = []

    for m in re.finditer(re_comp_ref, raw):
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
re_feature = [  "^ {5}(?P<feature_key>[\w]+)",
                " +(?P<location>.+)\n",
                "(?P<qualifiers>(?:^ {21}.*\n)*)"
]
re_feature = "".join(re_feature)
re_comp_feature = re.compile(re_feature, flags=re.M)


# Qualifers can have tags with /'s in the value so it's tough to escape them
# for now we need to split on "                     /"

def parseFeatures(raw, is_ordered=False):
    features_list = []

    for feature_match in re.finditer(re_comp_feature, raw):
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

def parseOrigin(raw):
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
