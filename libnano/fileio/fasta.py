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
libnano.fileio.fasta
~~~~~~~~~~~~~~~~~~~~

Read and write FASTA files

'''
import io
import re
from typing import (
    Dict,
    Iterator,
    List,
    Tuple,
    Union,
)

import _io  # type: ignore

_fasta_re: str = r'>(.+)[\n\r]+((?:[^>\S]|[{}]+[\n\r]+)+)'

_ALPHABETS: Dict[str, str] = {
    'DNA': r'ACGTWSMKRYBDHVN',
    'RNA': r'ACGUWSMKRYBDHVN',
    'AMINO_ACID': r'ABCDEFGHIJKLMNOPQRSTUVWYZX*-',
}

Record_T = Union[Tuple[str, str], Tuple[bytes, bytes]]


def sanitizeRecSeq(
        rec_id: str,
        rec_seq: str,
        alphabet: str,
        not_allowed: str,
):
    '''Sanitize sequence record.  Replaces bad characters

    Args:
        rec_id: Record ID
        rec_seq: Record sequence
        alphabet: {'DNA', 'RNA', 'AMINO_ACID'}
        not_allowed: String of characters not allow in the sequence

    Raises:
        ValueError: Illegal character
    '''
    if alphabet is None and not_allowed is not None:
        for idx, c in enumerate(rec_seq):
            if c in not_allowed:
                raise ValueError(
                    f'Record: {rec_id}\n Contains illegal char {c} at '
                    f'position {idx}',
                )
    if alphabet:
        ab = _ALPHABETS.get(alphabet.upper(), alphabet)
        replace_str = ''

        if not_allowed is not None:
            for c in not_allowed:
                ab = ab.replace(c, replace_str)
        for idx, c in enumerate(rec_seq):
            if c not in ab:
                raise ValueError(
                    f'Record: {rec_id}\n Contains illegal char {c} at '
                    f'position {idx}',
                )


def parseFasta(
        fasta_fn: str,
        alphabet: str = 'DNA',
        not_allowed: str = '',
) -> List[Tuple[str, str]]:
    '''Parse a standard fasta file and return a list of tuples containing
    the record id, sequence. Optionally check each sequence against an alphabet
    (DNA, RNA, AMINO_ACID, or a custom alphabet) and/or against a list of
    characters that are not allowed:

        alphabet = 'DNA' and not_allowed = 'N' - check against the DNA alphabet
                                                  but do not allow degeneracy
        alphabet = 'ATGC' - custom alphabet (insure that every character in
                            the sequence is A, T, G, or C)

    If a sequence fails the alphabet / not_allowed check a ValueError is
    raised.

    Args:
        fasta_fn: FASTA file path
        alphabet: {'DNA', 'RNA', 'AMINO_ACID'}
        not_allowed: String of characters not allow in the sequence

    Returns:
        List of tuples of the form::

            <record id>, <sequence>

    '''
    re_comp = re.compile(
        _fasta_re.format(
            _ALPHABETS.get(alphabet, alphabet),
        ),
    )

    with io.open(fasta_fn, 'r', encoding='utf-8') as fd:
        d = [
            (
                match.group(1).strip(),
                ''.join(match.group(2).strip().split()),
            )
            for match in re.finditer(re_comp, fd.read())
        ]
        if alphabet or not_allowed:
            for rec_id, rec_seq in d:
                sanitizeRecSeq(
                    rec_id,
                    rec_seq,
                    alphabet,
                    not_allowed,
                )
    return d


def parseFastaGen(
    fasta_fn: str,
    alphabet: str = 'DNA',
    not_allowed: str = '',
) -> Iterator[Tuple[str, str]]:
    '''Iterator that yields parsed records (ID, sequence) from a FASTA
    file.

    Args:
        fasta_fn: FASTA file path
        alphabet: {'DNA', 'RNA', 'AMINO_ACID'}
        not_allowed: String of characters not allow in the sequence

    Yields:
        Tuples of the form::

            <record id>, <sequence>

    '''

    rec_id = ''
    rec_seq = ''
    start_record_delim = '>'
    end_record_delim = ''

    def join_base():
        return ''
    split_str = ' '

    with io.open(fasta_fn, 'r', encoding='utf-8') as fd:
        for line in fd:
            if start_record_delim in line:
                if rec_id != end_record_delim:
                    if alphabet or not_allowed:
                        sanitizeRecSeq(
                            rec_id,
                            rec_seq,
                            alphabet,
                            not_allowed,
                        )
                    yield rec_id, rec_seq
                rec_id = line.strip().split(split_str)[0][1:]
                rec_seq = join_base()
            else:
                rec_seq += line.strip()
    if alphabet or not_allowed:
        sanitizeRecSeq(
            rec_id,
            rec_seq,
            alphabet,
            not_allowed,
        )
    yield rec_id, rec_seq


def write(
        fasta_fn: str,
        records: Record_T,
        iotype: str,
        wrap_len: int = 60,
) -> None:
    '''Write records to FASTA file

    Args:
        fasta_fn: the filename to write to
        records: a list of tuples of form::

            (rec_id, rec_seq)
        iotype: string 'unicode' or 'bytes'
        wrap_len: length to wrap the texts

    '''
    if not isinstance(iotype, str):
        raise TypeError('Unsupported type of iotype')
    if iotype not in ['unicode', 'bytes']:
        raise TypeError(f'Unsupported iotype {iotype}')
    if not isinstance(records, list):
        raise TypeError('records needs to be a list')
    if iotype == 'unicode':
        with io.open(fasta_fn, 'w', encoding='utf-8') as fd:
            for record in records:
                writeRecord(
                    fd,
                    records[0],
                    records[1],
                    wrap_len,
                )
    else:
        with io.open(fasta_fn, 'wb') as fd:
            for record in records:
                writeRecordB(
                    fd,
                    records[0],
                    records[1],
                    wrap_len,
                )


def writeRecord(
        fd: _io.TextIOWrapper,
        rec_id: str,
        rec_seq: str,
        wrap_len: int,
) -> None:
    '''Write record to file descriptor

    Args:
        fd: File descriptor
        rec_id: Record ID
        rec_seq: Record sequence
        wrap_len: Wrap length

    '''
    fd.write('>' + rec_id + '\n')
    for i in range(0, len(rec_seq), wrap_len):
        fd.write(rec_seq[i:i + wrap_len])
        fd.write('\n')


def writeRecordB(
        fd: _io.TextIOWrapper,
        rec_id: bytes,
        rec_seq: bytes,
        wrap_len: int,
):
    '''Write record to binary file descriptor

    Args:
        fd: File descriptor
        rec_id: Record ID
        rec_seq: Record sequence
        wrap_len: Wrap length

    '''
    fd.write(b'>' + rec_id + b'\n')
    for i in range(0, len(rec_seq), wrap_len):
        fd.write(rec_seq[i:i + wrap_len])
        fd.write(b'\n')
