# -*- coding: utf-8 -*-
import sys
from typing import (
    List,
    Dict,
    Tuple
)
from enum import IntEnum

IS_WINDOWS: bool = sys.platform == 'win32'

if IS_WINDOWS:
    '''On Windows, calling init() will filter ANSI escape sequences out of any
    text sent to stdout or stderr, and replace them with equivalent Win32
    calls. Otherwise no color!!!!'''
    from colorama import init as colorama_init
    colorama_init()
    from pygments.formatters import TerminalFormatter as TermFormatter
    NEWLINE_STR: str = '\r\n'
else:
    from pygments.formatters import TerminalTrueColorFormatter as TermFormatter
    NEWLINE_STR: str = '\n'
from pygments.style import Style
from pygments.token import Text as PygText
from pygments import highlight
from pygments.lexer import RegexLexer
from ssw import (
    Alignment,
    force_align
)

# import __init__
from libnano.seqstr import (
    reverseComplement,
    complement,
    reverse
)

class PrimeEnum(IntEnum):
    FIVE = 0
    THREE = 2
    BLUNT = 4

PRIME_ENUM_MAP = {
    PrimeEnum.FIVE : 'five',
    PrimeEnum.THREE : 'three',
    PrimeEnum.BLUNT : 'blunt'
}

TRANTAB: Dict[int, int] = str.maketrans('acgt', 'ACGT')

# NOTE: For some reason windows higlighting needs to invert the
# case of the REGEX
WIN_TOKENS: dict = {
    'root': [
        (r'[ACGT]', PygText)
    ]
}

POSIX_TOKENS: dict = {
    'root': [
        (r'[acgt]', PygText)
    ]
}

class DNALex(RegexLexer):
    name: str = 'DNALex'
    aliases: List[str] = ['dna']
    filenames: List[str] = ['*.dna']

    tokens: dict = WIN_TOKENS if IS_WINDOWS else POSIX_TOKENS
# end class

class MisMatchStyle(Style):
    default_style: str = ""
    styles: dict = {
        PygText: 'bold nounderline #ff0000',
    }
# end class

def align_complement(fwd: str, rev: str) -> Tuple[Alignment, str]:
    rc_rev: str = reverseComplement(rev)
    alignment: Alignment = force_align(rc_rev, fwd)
    return alignment, rc_rev
# end def

def five_prime_type(alignment: Alignment, fwd: str, rev: str):
    fwd_idx0: int = alignment.reference_start
    rev_idx0: int = alignment.read_start
    if fwd_idx0 > rev_idx0:
        return PrimeEnum.FIVE, fwd[:fwd_idx0]
    elif fwd_idx0 < rev_idx0:
        return PrimeEnum.THREE, reverse(reverse(rev)[:rev_idx0])
    else:
        return PrimeEnum.BLUNT, ''

def three_prime_type(alignment: Alignment, fwd: str, rev: str):
    # fwd_idx0: int = alignment.reference_start
    # rev_idx0: int = alignment.read_start
    fwd_idx1: int = alignment.reference_end
    rev_idx1: int = alignment.read_end
    max_idx_fwd: int = len(fwd) - 1
    max_idx_rev: int = len(rev) - 1

    if fwd_idx1 < max_idx_fwd:
        return PrimeEnum.THREE, fwd[fwd_idx1+1:]
    elif rev_idx1 < max_idx_rev:
        # print('$$$', rev_idx1, max_idx_rev)
        return PrimeEnum.FIVE, reverse(reverse(rev)[rev_idx1+1:])
    else:
        return PrimeEnum.BLUNT, ''
# end def

def string_align_complement(
        fwd: str,
        rev: str,
        alignment: Alignment = None,
        rc_rev: str = None,
        do_print: bool = False,
        do_highlight: bool = False) -> Tuple[str, str]:
    if alignment is None:
        alignment, rc_rev = align_complement(fwd, rev)
    if rc_rev is None:
        rc_rev = reverseComplement(rev)
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


            highlight_unformat_rev: str = ''.join(highlight_rev_list)
            highlight_rev: str = highlight(   highlight_unformat_rev,
                                    DNALex(),
                                    TermFormatter(style=MisMatchStyle))
            out_rev: str = (
                buffer_rev +
                reverse_rev[:rev_lo_idx] +
                highlight_rev.strip().translate(TRANTAB) +
                reverse_rev[rev_lo_idx+total_delta:]
            )
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
