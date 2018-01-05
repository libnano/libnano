
try:
    import libnano.helpers.textwrap as textwrap
except:
    import textwrap

from collections import OrderedDict
from libc.stdio cimport FILE, fopen, fclose, fwrite, fflush, fprintf, sprintf
from libnano from helpers cimport c_util
from libc.stdlib cimport malloc, free

cdef inline fwrite_b(FILE* fd, bytes out_str_b):
    """ helper function for fwrite
    """
    cdef char* out_str_c = out_str_b
    cdef int len_out = len(out_str_b)
    fwrite(out_str_c, sizeof(char), len_out, fd)

cdef write(FILE* fd, d):
    d_info = d[b'info']
    writeLocus(fd, d_info)
    writeDefinition(fd, d_info)
    writeAccession(fd, d_info)
    writeVersion(fd, d_info)
    writeDBLINK(fd, d_info)
    writeKeywords(fd, d_info)
    writeSource(fd, d_info)
    i = 0
    for reference in d_info[b'references']:
        i += 1
        writeReference(fd, reference, i)
    writeComment(fd, d_info)
    writeFeatures(fd, d)
    writeOrigin(fd, d)
# end def

def write_file(filepath, d):
    cdef char* filepath_c
    cdef FILE* fd = NULL
    cdef bytes filepath_b = c_util._bytes(filepath)

    try:
        filepath_c = filepath_b
        fd = fopen(filepath_c, 'wb')
        if fd == NULL:
            raise OSError("file %s could not be open", filepath.encode('utf-8'))
        else:
            write(fd, d)
    finally:
        if fd != NULL:
            fclose(fd)
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

def multiline(main_str, indent_str, lim=67):
    o_list = [ main_str[i:i+lim] for i in range(0, len(main_str), lim) ]
    return indent_str.join(o_list)

def multiline_spaces(main_str, indent_str, lim=67):
    """ python textwrap is incompatible with bytestring
    maybe port python2s version?
    """
    #indent_str = indent_str.decode('utf-8')
    #main_list = textwrap.wrap(main_str.decode('utf-8'), lim, drop_whitespace=True,
    #    break_on_hyphens=False)
    #return indent_str.join(main_list).encode('utf-8')

    main_list = textwrap.wrap(main_str, lim, drop_whitespace=True,
        break_on_hyphens=False)
    return indent_str.join(main_list)
#end def

cdef writeLocus(FILE* fd, d):
    """
    Locus in sample file is not valid
    """
    cdef char* name = c_util.obj_to_cstr(d[b'name'])
    cdef char* molecule_type = c_util.obj_to_cstr(d[b'molecule_type'])
    cdef char* gb_division = c_util.obj_to_cstr(d[b'gb_division'])
    cdef char* mod_date = c_util.obj_to_cstr(d[b'mod_date'])
    cdef int length = d[b'length']
    cdef char* form = NULL
    cdef object form_obj = d[b'form']
    cdef const char* out_str = "LOCUS       %-16s %11d bp    %-6s  %8s %s %s\n"

    if form_obj is None:
        form = ''
    elif form_obj == b'linear':
        form = 'linear  '

    fprintf(fd, out_str,
            name, length,
            molecule_type, form,
            gb_division, mod_date)
# end def

cdef writeDefinition(FILE* fd, d):
    cdef bytes out_str
    cdef bytes definition_str = b"DEFINITION  "
    cdef bytes indent_str = b"\n            "

    multiline_spaces(d[b'definition'], indent_str)
    cdef list out_list = [definition_str,
                multiline_spaces(d[b'definition'], indent_str),
                b'\n']
    cdef bytes out_str_b = b''.join(out_list)
    fwrite_b(fd, out_str_b)
# end def

cdef writeAccession(FILE* fd, d):
    cdef char* accession = c_util.obj_to_cstr(d[b'accession'])
    fprintf(fd, "ACCESSION   %s\n", accession)
# end def

cdef writeVersion(FILE* fd, d):
    cdef bytes version = d[b'version']
    cdef char* version_c
    cdef char* gi
    if version is not None:
        version_c = version
        gi = c_util.obj_to_cstr(d[b'GI'])
        fprintf(fd, "VERSION     %s  GI:%s\n", version, gi)
# end def

cdef writeDBLINK(FILE* fd, d):
    cdef char* dblink
    if 'dblink' in d:
        if d['dblink'] is not None:
            dblink = c_util.obj_to_cstr(d[b'dblink'])
            fprintf(fd, "DBLINK      %s\n", dblink)
# end def

cdef writeKeywords(FILE* fd, d):
    cdef char* keywords = c_util.obj_to_cstr(d[b'keywords'])
    fprintf(fd, "KEYWORDS    %s\n", keywords)
# end def

cdef writeSource(FILE* fd, d):
    cdef bytes out_str_b
    source_str = b"SOURCE      "
    organism_str = b"\n  ORGANISM  "
    indent_str = b"\n            "
    org = d[b'organism']

    out_list = [source_str, d[b'source'],
                organism_str,
                org[0]]
    if org[1] is not None:
        out_list += [   indent_str,
                        multiline_spaces(org[1], indent_str)]
    out_list += [b'\n']
    out_str_b = b''.join(out_list)
    fwrite_b(fd, out_str_b)
# end def


cdef writeReference(FILE* fd, ref, i):
    cdef bytes out_str_b
    assert(ref[b'r_index'] == i)
    reference_str = b"REFERENCE   "
    authors_str = b"  AUTHORS   "
    title_str = b"  TITLE     "
    indent_str = b"\n            "
    journal_str = b"  JOURNAL   "
    pubmed_str = b"  PUBMED    "
    if ref[b'start_idx'] is not None:
        idx_str = "%d  (bases %d to %d)\n" % (i, ref[b'start_idx'], ref[b'end_idx'])
    else:
        idx_str = '%d\n'
    out_list = [reference_str,
            idx_str.encode('utf-8'),
            authors_str,
            multiline_spaces(ref[b'authors'], indent_str), b'\n',
            title_str,
            multiline_spaces(ref[b'title'], indent_str), b'\n',
            journal_str,
            multiline_spaces(ref[b'journal_info'], indent_str), b'\n',
            ]
    if ref[b'pubmed'] is not None:
        out_list += [pubmed_str, ref[b'pubmed'], b'\n']

    out_str_b = b''.join(out_list)
    fwrite_b(fd, out_str_b)
# end def

cdef writeComment(FILE* fd, d):
    cdef bytes out_str_b
    comment_str = b"COMMENT     "
    indent_str = b"\n            "
    indent_str_gb_asm = b"            "
    if b'comment' in d:
        comment = d[b'comment']
        if isinstance(comment, list):
            out_list = [comment_str,
                multiline_spaces(comment[0], indent_str),
                indent_str, indent_str, indent_str,
                indent_str.join(comment[1]), b'\n'
                ]
        else:
            out_list = [comment_str, multiline_spaces(comment, indent_str), '\n']

        out_str_b = b''.join(out_list)
        fwrite_b(fd, out_str_b)
# end def


cdef bytes pysprint_key_val(const char* fmt, bytes key, bytes val, int size):
    cdef char fixed_buf[1024]
    cdef char* buf = fixed_buf
    cdef bytes out
    cdef char* key_c = key
    cdef char* val_c = val
    if size > 1024:
        buf = <char*> malloc((size+80)*sizeof(char))
        sprintf(buf, fmt, key_c, val_c)
        out = c_util.cstr_to_obj_nolength(buf, 1)
        free(buf)
    else:
        sprintf(buf, fmt, key_c, val_c)
        out = c_util.cstr_to_obj_nolength(buf, 1)
    return out
# end def

cdef bytes pysprint_key_val_d(const char* fmt, bytes key, int val, int size):
    """ key unused for now
    """
    cdef char buf[256]
    cdef bytes out
    cdef char* key_c = key
    sprintf(buf, fmt, key_c, val)
    out = c_util.cstr_to_obj_nolength(buf, 1)
    return out
# end def

cdef writeFeatures(FILE* fd, d):
    cdef bytes out_str_b
    cdef bytes temp

    cdef bytes feature_header = b"FEATURES             Location/Qualifiers\n"
    fwrite_b(fd, feature_header)
    feature_type_prefix_str = b"     "
    indent_str = b"\n                     "
    feature_type_field_size = 16
    #d_features = d[b'features']
    for feature in d[b'features']:
        ftype = feature[b'type']
        location_str = feature[b'location']
        out_list = [feature_type_prefix_str,
                    ftype, (feature_type_field_size - len(ftype))*b" ",
                    location_str, b'\n'
                    ]
        for key, value_list in feature[b'qualifiers'].items():
            if not isinstance(value_list, list):
                value_list = [value_list]
            for value in value_list:
                if key in (b'codon_start', b'transl_table'):    # codon start is a weird edge case
                    temp = pysprint_key_val_d("                     /%s=%d\n",
                                                key, value, len(key))
                    out_list.append(temp)
                elif key == b'transl_except':
                    temp = pysprint_key_val("                     /%s=%s\n",
                                            key, value, len(key) + len(value))
                    out_list.append(temp)
                elif key == b'translation':
                    qualifier_str = pysprint_key_val("/%s=\"%s\"",
                                                key, value, len(key) + len(value))
                    qual_list = multiline(qualifier_str, indent_str, lim=58) # or 58
                    out_list += [b"                     " + qual_list, b'\n']
                else:
                    if value is None:
                        qualifier_str = pysprint_key_val("/%s%s\n",
                                                    key, b"", len(key))
                    else:
                        qualifier_str = pysprint_key_val("/%s=\"%s\"\n",
                                                key, value, len(key) + len(value))
                    qual_list_str = multiline_spaces(qualifier_str, indent_str, lim=58)  #or 58
                    out_list += [b"                     " + qual_list_str, b'\n']

        out_str_b = b''.join(out_list)
        fwrite_b(fd, out_str_b)
# end def

cdef writeOrigin(FILE* fd, d):
    cdef bytes temp_b
    cdef int k
    cdef int i = 0
    cdef int len_temp
    cdef char* temp
    cdef char* tup_c[6]
    cdef list list_strb = [None]*6
    cdef bytes origin = d[b'seq']
    cdef bytes origin_header = b"ORIGIN      \n"

    fwrite_b(fd, origin_header)

    count = len(origin) - 60

    while i < count:
        #str_tup = (i+1, ) + tuple([origin[j:j+10] for j in range(i, i+60, 10) ])
        k = 0
        for j in range(i, i + 60, 10):
            list_strb[k] = origin[j:j + 10]
            tup_c[k] = c_util.obj_to_cstr(list_strb[k])
            k += 1
        fprintf(fd, "%9d %s %s %s %s %s %s\n",
                i + 1,
                tup_c[0],
                tup_c[1],
                tup_c[2],
                tup_c[3],
                tup_c[4],
                tup_c[5]
                )
        i += 60
    # end while

    # do last line in unicode for ease of using % formatting due to variable
    # line length
    last_str = origin[i:].decode('utf-8')
    last_str = last_str + (60-(len(last_str)))*" "
    last_list = [last_str[j:j+10] for j in range(0, 60, 10) ]
    last_list = [x.rstrip() for x in last_list if x[0] != " "]
    str_tup = (i+1, ) + tuple(last_list)

    format_string = '%9d' + " %s"*len(last_list) + "\n//\n"

    temp_b = (format_string % str_tup).encode('utf-8')
    fwrite_b(fd, temp_b)
# end def

if __name__ == '__main__':
    import filecmp
    import gb_reader_b as gbr
    import os.path as opath
    path = opath.dirname(opath.dirname(opath.dirname(opath.abspath(__file__))))
    # fn = opath.join(path, "tests", "test_data", "sample.gb")
    fn = opath.join(path, "tests", "test_data", "mds42_full.gb")
    # fn = opath.join(path, "tests", "test_data", "failed.gb")
    # fn = opath.join(path, "tests", "test_data", "sample_complex.gb")
    fn_out = opath.join(path, "tests", "test_data", "sample_out.gb")
    gbd = gbr.parse(fn, is_ordered=True)
    # print(gbd['references'])
    write_file(fn_out, gbd)
    print(filecmp.cmp(fn, fn_out))
