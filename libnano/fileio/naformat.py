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
libnano.fileio.naformat
~~~~~~~~~~~~~~~~~~~~~~~

'''
import sys
from enum import IntEnum
from typing import (
    Dict,
    List,
    Tuple,
)

from pygments import highlight
from pygments.lexer import RegexLexer
from pygments.style import Style
from pygments.token import Text as PygText
from ssw import (  # type: ignore
    Alignment,
    force_align,
)

from libnano.seqstr import (  # type: ignore
    reverse_complement,
    reverse_seq,
)

IS_WINDOWS: bool = sys.platform == 'win32'

if IS_WINDOWS:
    '''On Windows, calling init() will filter ANSI escape sequences out of any
    text sent to stdout or stderr, and replace them with equivalent Win32
    calls. Otherwise no color!!!!'''
    from colorama import init as colorama_init  # type: ignore
    colorama_init()
    from pygments.formatters import TerminalFormatter as TermFormatter
    NEWLINE_STR = '\r\n'
else:
    from pygments.formatters import TerminalTrueColorFormatter
    TermFormatter = TerminalTrueColorFormatter  # type: ignore
    NEWLINE_STR = '\n'  # type: ignore


class PrimeEnum(IntEnum):
    FIVE = 0
    THREE = 2
    BLUNT = 4


PRIME_ENUM_MAP = {
    PrimeEnum.FIVE: 'five',
    PrimeEnum.THREE: 'three',
    PrimeEnum.BLUNT: 'blunt',
}

TRANTAB_DICT = str.maketrans('acgt', 'ACGT')

# NOTE: For some reason windows higlighting needs to invert the
# case of the REGEX
WIN_TOKENS_DICT = {
    'root': [
        (r'[ACGT]', PygText),
    ],
}

POSIX_TOKENS: dict = {
    'root': [
        (r'[acgt]', PygText),
    ],
}


class DNALex(RegexLexer):
    '''DNALex class

    '''
    name: str = 'DNALex'
    aliases: List[str] = ['dna']
    filenames: List[str] = ['*.dna']

    tokens: dict = WIN_TOKENS_DICT if IS_WINDOWS else POSIX_TOKENS


class MisMatchStyle(Style):
    '''MisMatchStyle class

    '''
    default_style: str = ''
    styles: dict = {
        PygText: 'bold nounderline #ff0000',
    }


def align_complement(
        fwd: str,
        rev: str,
) -> Tuple[Alignment, str]:
    '''Align complement ``fwd`` and ``rev`` strands

    Args:
        fwd: Forward strand
        rev: Reverse strand

    Returns:
        Tuple of the form::

            <Alignment instance>, <reverse strand reverse complement>

    '''
    rc_rev = reverse_complement(rev)
    alignment: Alignment = force_align(rc_rev, fwd)
    return alignment, rc_rev


def five_prime_type(
        alignment: Alignment,
        fwd: str,
        rev: str,
) -> Tuple[int, str]:
    '''Five prime type tuple

    Args:
        alignment: class:`ssw.alignmenttuple.Alignment` instance
        fwd: Forward strand
        rev: Reverse strand

    Returns:
        Tuple of the form::
            <PrimeEnum>, <strand offset>
    '''
    fwd_idx0: int = alignment.reference_start
    rev_idx0: int = alignment.read_start
    if fwd_idx0 > rev_idx0:
        return PrimeEnum.FIVE, fwd[:fwd_idx0]
    elif fwd_idx0 < rev_idx0:
        return PrimeEnum.THREE, reverse_seq(reverse_seq(rev)[:rev_idx0])
    else:
        return PrimeEnum.BLUNT, ''


def three_prime_type(
        alignment: Alignment,
        fwd: str,
        rev: str,
) -> Tuple[int, str]:
    '''Three prime type tuple

    Args:
        alignment: class:`ssw.alignmenttuple.Alignment` instance
        fwd: Forward strand
        rev: Reverse strand

    Returns:
        Tuple of the form::
            <PrimeEnum>, <strand offset>
    '''
    fwd_idx1: int = alignment.reference_end
    rev_idx1: int = alignment.read_end
    max_idx_fwd: int = len(fwd) - 1
    max_idx_rev: int = len(rev) - 1

    if fwd_idx1 < max_idx_fwd:
        return PrimeEnum.THREE, fwd[fwd_idx1 + 1:]
    elif rev_idx1 < max_idx_rev:
        return PrimeEnum.FIVE, reverse_seq(reverse_seq(rev)[rev_idx1 + 1:])
    else:
        return PrimeEnum.BLUNT, ''


def string_align_complement(
        fwd: str,
        rev: str,
        alignment: Alignment = None,
        rc_rev: str = '',
        do_print: bool = False,
        do_highlight: bool = False,
) -> Tuple[str, str]:
    '''String alignment complement

    Args:
        fwd: Forward strand
        rev: Reverse strand
        alignment: class:`ssw.alignmenttuple.Alignment` instance
        rc_rev: Optional reverse complement of the `rev`
        do_print: If True, print results
        do_highlight: If True, highlight print results

    Returns:
        Tuple of the form::
            <PrimeEnum>, <strand offset>
    '''
    if alignment is None:
        alignment, rc_rev = align_complement(fwd, rev)
    if not rc_rev:
        rc_rev = reverse_complement(rev)
    reverse_rev: str = reverse_seq(rev)

    fwd_idx0: int = alignment.reference_start
    # fwd_idx1: int = alignment.reference_end
    rev_idx0: int = alignment.read_start
    # rev_idx1: int = alignment.read_end

    max_delta_fwd: int = len(fwd) - fwd_idx0
    max_delta_rev: int = len(rev) - rev_idx0

    lim_hi = max_delta_rev
    if max_delta_fwd < max_delta_rev:
        lim_hi = max_delta_fwd

    lo_delta: int
    buffer_fwd = ''
    buffer_rev = ''
    if fwd_idx0 < rev_idx0:
        lo_delta = fwd_idx0
        buffer_fwd = ' ' * (rev_idx0 - fwd_idx0)
    else:
        lo_delta = rev_idx0
        buffer_rev = ' ' * (fwd_idx0 - rev_idx0)

    out_rev: str = buffer_rev + reverse_rev
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
            highlight_rev: str = highlight(
                highlight_unformat_rev,
                DNALex(),
                TermFormatter(style=MisMatchStyle),
            )
            out_rev = (
                buffer_rev +
                reverse_rev[:rev_lo_idx] +
                highlight_rev.strip().translate(TRANTAB_DICT) +
                reverse_rev[rev_lo_idx + total_delta:]
            )

    out_fwd: str = buffer_fwd + fwd
    if do_print:
        print(out_fwd)
        print(out_rev)
    return out_fwd, out_rev


if __name__ == '__main__':
    def printer(x: int, fwd: str, rev: str) -> int:
        print('%d. ' % x)
        string_align_complement(fwd, rev, do_print=True, do_highlight=True)
        return x + 1
    i = 0
    i = printer(i, 'GGATCCAAA', 'TTTGGATC')
    i = printer(i, 'GGATCCAAA', 'TTGGATC')
    i = printer(i, 'GGATCCAAA', 'CTTGGATC')
    i = printer(i, 'GGATCCAAA', 'TTTGGATCAAAAAA')
    i = printer(i, 'GGATCCAAA', 'TTTGGATCAAAAAA')
    i = printer(i, 'GGATCCAACCCCCC', 'TTTGGATCAAAAAA')
    i = printer(i, 'GGATCCAAA', 'GGGGGATTGGATC')
    i = printer(i, 'G' * 10 + 'GGATCCAAA', 'TTTGGATC')
    i = printer(i, 'G' * 10 + 'GGATCCAAA', 'ATTTGCATC')
    string_align_complement('G' * 10 + 'GGATCCAAA', 'ATTTGCATC', do_print=True)
