# -*- coding: utf-8 -*-
import sys
from typing import (
    List,
    Dict,
    Tuple
)

if sys.platform == 'win32':
    '''On Windows, calling init() will filter ANSI escape sequences out of any
    text sent to stdout or stderr, and replace them with equivalent Win32
    calls. Otherwise no color!!!!'''
    from colorama import init as colorama_init
    colorama_init()
    from pygments.formatters import TerminalFormatter as TermFormatter
else:
    from pygments.formatters import TerminalTrueColorFormatter as TermFormatter
from pygments.style import Style
from pygments.token import Text as pygText
from pygments import highlight
from pygments.lexer import RegexLexer
from ssw import (
    Alignment,
    force_align
)

import __init__
from libnano.seqstr import (
    reverseComplement,
    complement,
    reverse
)

TRANTAB: Dict[int, int] = str.maketrans('acgt', 'ACGT')

class DNALex(RegexLexer):
    name: str = 'DNALex'
    aliases: List[str] = ['dna']
    filenames: List[str] = ['*.dna']

    tokens: dict = {
        'root': [
            (r'[ACGT]', pygText)
        ]
    }
# end class

class MisMatchStyle(Style):
    default_style: str = ""
    styles: dict = {
        pygText: 'bold nounderline #ff0000',
    }
# end class

def align_complement(fwd: str, rev: str) -> Tuple[Alignment, str]:
    rc_rev: str = reverseComplement(rev)
    alignment: Alignment = force_align(rc_rev, fwd)
    return alignment, rc_rev

def string_align_complement(
        fwd: str,
        rev: str,
        do_print: bool = False,
        do_highlight: bool = False) -> Tuple[str, str]:
    rc_rev: str
    alignment: Alignment
    alignment, rc_rev = align_complement(fwd, rev)
    reverse_rev: str = reverse(rev)

    fwd_idx0: int = alignment.reference_start
    fwd_idx1: int = alignment.reference_end
    rev_idx0: int = alignment.read_start
    rev_idx1: int = alignment.read_end

    max_delta_fwd: int = len(fwd) - fwd_idx0
    max_delta_rev: int = len(rev) - rev_idx0
    if max_delta_fwd < max_delta_rev:
        lim_hi: int = max_delta_fwd
    else:
        lim_hi: int = max_delta_rev

    lo_delta: int
    buffer_fwd: str = ''
    buffer_rev: str = ''
    if fwd_idx0 < rev_idx0:
        lo_delta = fwd_idx0
        buffer_fwd = ' '*(rev_idx0 - fwd_idx0)
    else:
        lo_delta = rev_idx0
        buffer_rev = ' '*(fwd_idx0 - rev_idx0)

    if do_highlight:
        highlight_rev_list: List[str] = []
        fwd_lo_idx: int = fwd_idx0 - lo_delta
        rev_lo_idx: int = rev_idx0 - lo_delta
        total_delta: int = lo_delta + lim_hi

        for i in range(total_delta):
            reverse_rev_base = reverse_rev[rev_lo_idx + i]
            if rc_rev[rev_lo_idx + i] != fwd[fwd_lo_idx + i]:
                highlight_rev_list.append(reverse_rev_base.lower())
            else:
                highlight_rev_list.append(reverse_rev_base)


            highlight_unformat_rev: str = (
                buffer_rev +
                reverse_rev[:rev_lo_idx] +
                ''.join(highlight_rev_list) +
                reverse_rev[rev_lo_idx+total_delta:]
            )
            out_rev: str = highlight(   highlight_unformat_rev,
                                    DNALex(),
                                    TermFormatter(style=MisMatchStyle))
            out_rev = out_rev.strip().translate(TRANTAB)
    else:
        out_rev: str = buffer_rev + reverse_rev
    out_fwd: str = buffer_fwd + fwd
    if do_print:
        print(out_fwd)
        print(out_rev)
    return out_fwd, out_rev
# end def

if __name__ == '__main__':
    def printer(x: int, fwd: str, rev:str) -> int:
        print("%d. " % x)
        string_align_complement(fwd, rev, do_print=True, do_highlight=True)
        return x + 1
    i = 0
    i = printer(i, "GGATCCAAA", "TTTGGATC")
    i = printer(i, "GGATCCAAA", "TTGGATC")
    i = printer(i, "GGATCCAAA", "CTTGGATC")
    i = printer(i, "GGATCCAAA", "TTTGGATCAAAAAA")
    i = printer(i, "GGATCCAAA", "TTTGGATCAAAAAA")
    i = printer(i, "GGATCCAACCCCCC", "TTTGGATCAAAAAA")
    i = printer(i, "GGATCCAAA", "GGGGGATTGGATC")
    i = printer(i, "G"*10 +"GGATCCAAA", "TTTGGATC")
    i = printer(i, "G"*10 +"GGATCCAAA", "ATTTGCATC")
    string_align_complement("G"*10 +"GGATCCAAA", "ATTTGCATC", do_print=True)
