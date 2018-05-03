# -*- coding: utf-8 -*-
import io
import textwrap

def write(fd: '_io.TextIOWrapper', d: dict, order_qualifiers: bool = False):
    d_info = d['info']
    writeLocus(fd, d_info)
    writeDefinition(fd, d_info)
    writeAccession(fd, d_info)
    writeVersion(fd, d_info)
    writeDBLINK(fd, d_info)
    writeKeywords(fd, d_info)
    writeSource(fd, d_info)
    i = 0
    for reference in d_info['references']:
        i += 1
        writeReference(fd, reference, i)
    writeComment(fd, d_info)
    writeFeatures(fd, d, order_qualifiers)
    writeOrigin(fd, d)
# end def

def write_file(filepath: str, d: dict, order_qualifiers: bool = False):
    with io.open(filepath, 'w', encoding='utf-8') as fd:
        write(fd, d, order_qualifiers=order_qualifiers)
# end def

"""
http://www.insdc.org/files/feature_table.html

section 4.1 shows a 79 character feature table format, as well as the LOCUS
being 80 characters but BioPython files output an 80 character version

However section 4.3 Data item positions shows

    22-80              location

as the 80th character is in play.

This code sticks with 79 character limit
"""

def multiline(main_str: str, indent_str: str, lim: int = 67) -> str:
    o_list = [ main_str[i:i+lim] for i in range(0, len(main_str), lim) ]
    return indent_str.join(o_list)

def multiline_spaces(main_str: str, indent_str: str, lim: int = 67) -> str:
    main_list = textwrap.wrap(main_str, lim, drop_whitespace=True,
        break_on_hyphens=False)
    return indent_str.join(main_list)
#end def

def writeLocus(fd: '_io.TextIOWrapper', d: dict):
    """
    Locus in sample file is not valid
    """
    locus_str = "LOCUS       "
    name = d['name']
    molecule_type = d['molecule_type']
    form = d['form']
    if form is None:
        form = ''
    elif form == 'linear':
        form = 'linear  '

    gb_division = d['gb_division']
    if gb_division is None:
        gb_division = ''

    mod_date = d['mod_date']
    if mod_date is None:
        mod_date = ''

    out_str = locus_str + "%-16s %11d bp    %-6s  %8s %s %s\n" %\
        (name, d['length'], molecule_type, form, gb_division, mod_date)
    fd.write(out_str)
# end def

def writeDefinition(fd: '_io.TextIOWrapper', d: dict):
    definition_str = "DEFINITION  "
    indent_str = "\n            "
    out_list = [definition_str,
                multiline_spaces(d['definition'], indent_str),
                '\n']
    fd.write(''.join(out_list))
# end def

def writeAccession(fd: '_io.TextIOWrapper', d: dict):
    accession_str = "ACCESSION   "
    fd.write(''.join([accession_str, d['accession'], '\n']))
# end def

def writeVersion(fd: '_io.TextIOWrapper', d: dict):
    version_str = "VERSION     "
    version = d['version']
    if version is not None:
        gi = d['GI']
        fd.write(''.join([version_str, "%s  GI:%s" % (version, gi), '\n']))
# end def

def writeDBLINK(fd: '_io.TextIOWrapper', d: dict):
    if 'dblink' in d:
        if d['dblink'] is not None:
            fd.write("DBLINK      %s\n" % (d['dblink']))
# end def

def writeKeywords(fd: '_io.TextIOWrapper', d: dict):
    keywords_str = "KEYWORDS    "
    fd.write(''.join([keywords_str, d['keywords'], '\n']))
# end def

def writeSource(fd: '_io.TextIOWrapper', d: dict):
    source_str = "SOURCE      "
    organism_str = "\n  ORGANISM  "
    indent_str = "\n            "
    org = d['organism']
    out_list = [source_str, d['source'],
                organism_str,
                org[0]]
    if org[1] is not None:
        out_list += [   indent_str,
                        multiline_spaces(org[1], indent_str)]
    out_list += ['\n']
    fd.write(''.join(out_list))
# end def


def writeReference(fd: '_io.TextIOWrapper', ref: dict, i: int):
    assert(ref['r_index'] == i)
    reference_str = "REFERENCE   "
    authors_str = "  AUTHORS   "
    title_str = "  TITLE     "
    indent_str = "\n            "
    journal_str = "  JOURNAL   "
    pubmed_str = "  PUBMED    "
    if ref['start_idx'] is not None:
        idx_str ="  (bases %d to %d)" % (ref['start_idx'], ref['end_idx'])
    else:
        idx_str = ''
    out_list = [reference_str,
            "%d%s\n" % (i, idx_str),
            authors_str,
            multiline_spaces(ref['authors'], indent_str), '\n',
            title_str,
            multiline_spaces(ref['title'], indent_str), '\n',
            journal_str,
            multiline_spaces(ref['journal_info'], indent_str), '\n',
            ]
    if ref['pubmed'] is not None:
        out_list += [pubmed_str, ref['pubmed'], '\n']

    fd.write(''.join(out_list))
# end def

def writeComment(fd: '_io.TextIOWrapper', d: dict):
    comment_str = "COMMENT     "
    indent_str = "\n            "
    indent_str_gb_asm = "            "
    if 'comment' in d:
        comment = d['comment']
        if isinstance(comment, list):
            out_list = [comment_str,
                multiline_spaces(comment[0], indent_str),
                indent_str, indent_str, indent_str,
                indent_str.join(comment[1]), '\n'
                ]
        else:
            out_list = [comment_str, multiline_spaces(comment, indent_str), '\n']
        fd.write(''.join(out_list))
# end def


def writeFeatures(fd: '_io.TextIOWrapper', d: dict, order_qualifiers: bool):
    feature_header = "FEATURES             Location/Qualifiers\n"
    feature_type_prefix_str = "     "
    indent_str = "\n                     "
    feature_type_field_size = 16
    fd.write(feature_header)
    for feature in d['features']:
        ftype = feature['type']
        location_str = feature['location']
        out_list = [feature_type_prefix_str,
                    ftype, (feature_type_field_size - len(ftype))*" ",
                    location_str, '\n'
                    ]
        if order_qualifiers:
            quals = feature['qualifiers']
            for key in sorted(quals.keys()):
                value_list = quals[key]
                if not isinstance(value_list, list):
                    value_list = [value_list]
                for value in sorted(value_list): # assumes value_list is all the same type
                    if key in ('codon_start', 'transl_table'):    # codon start is a weird edge case
                        out_list += ["                     /%s=%d\n" % (key, value)]
                    elif key == 'transl_except':
                        out_list += ["                     /%s=%s\n" % (key, value)]
                    elif key == 'translation':
                        qualifier_str = "/%s=\"%s\"" % (key, value)
                        qual_list = multiline(qualifier_str, indent_str, lim=58) # or 58
                        out_list += ["                     " + qual_list, '\n']
                    else:
                        if value is None:
                            qualifier_str = "/%s\n" % (key)
                        else:
                            qualifier_str = "/%s=\"%s\"\n" % (key, value)
                        qual_list = multiline_spaces(qualifier_str, indent_str, lim=58)  #or 58
                        out_list += ["                     " + qual_list, '\n']
        else:
            for key, value_list in feature['qualifiers'].items():
                if not isinstance(value_list, list):
                    value_list = [value_list]
                for value in value_list:
                    if key in ('codon_start', 'transl_table'):    # codon start is a weird edge case
                        out_list += ["                     /%s=%d\n" % (key, value)]
                    elif key == 'transl_except':
                        out_list += ["                     /%s=%s\n" % (key, value)]
                    elif key == 'translation':
                        qualifier_str = "/%s=\"%s\"" % (key, value)
                        qual_list = multiline(qualifier_str, indent_str, lim=58) # or 58
                        out_list += ["                     " + qual_list, '\n']
                    else:
                        if value is None:
                            qualifier_str = "/%s\n" % (key)
                        else:
                            qualifier_str = "/%s=\"%s\"\n" % (key, value)
                        qual_list = multiline_spaces(qualifier_str, indent_str, lim=58)  #or 58
                        out_list += ["                     " + qual_list, '\n']
        fd.write(''.join(out_list))
# end def

def writeOrigin(fd: '_io.TextIOWrapper', d: dict):
    i = 0
    format_string = '%9d %s %s %s %s %s %s\n'
    origin_header = "ORIGIN      \n"
    fd.write(origin_header)
    origin = d['seq']
    count = len(origin) - 60
    while i < count:
        str_tup = (i+1, ) + tuple([origin[j:j+10] for j in range(i, i+60, 10) ])
        fd.write(format_string % str_tup)
        i += 60
    last_str = origin[i:]
    last_str = last_str + (60-(len(last_str)))*" "
    last_list = [last_str[j:j+10] for j in range(0, 60, 10) ]
    last_list = [x.rstrip() for x in last_list if x[0] != " "]
    str_tup = (i+1, ) + tuple(last_list)
    format_string = '%9d' + " %s"*len(last_list) + '\n'
    fd.write(format_string % str_tup)
    fd.write("//\n")
# end def

if __name__ == '__main__':
    import filecmp
    import gb_reader as gbr
    import os.path as opath
    path = opath.dirname(opath.dirname(opath.dirname(opath.abspath(__file__))))
    fn = opath.join(path, "tests", "test_data", "mds42_full.gb")
    # fn = opath.join(path, "tests", "test_data", "failed.gb")
    # fn = opath.join(path, "tests", "test_data", "sample_complex.gb")
    fn_out = opath.join(path, "tests", "test_data", "sample_out.gb")
    gbd = gbr.parse(fn, is_ordered=True)
    # print(gbd['references'])
    write_file(fn_out, gbd)
    print(filecmp.cmp(fn, fn_out))
