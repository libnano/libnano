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
libnano.fileio.__init__
~~~~~~~~~~~~~~~~~~~~~~~~
'''
import os.path as op

try:
    import libnano
except (ImportError, ModuleNotFoundError):
    import sys
    PACKAGE_PATH = op.dirname(op.dirname(op.dirname(op.abspath(__file__))))
    sys.path.append(PACKAGE_PATH)

from libnano.fileio import (
    fasta,
    gb_reader,
    gb_writer,
    gbtofasta,
    xmfa,
)


def getSeqFromFile(seq_fp: str) -> str:
    _, ext = op.splitext(seq_fp)
    if ext in ('.gb', '.gbk', '.genbank'):
        rec = gb_reader.parse(seq_fp)
        return str(rec['seq'])
    elif ext in ('.fa', '.fasta'):
        # Just returns the first sequence in the file for now --
        # this may be undesirable as a general behavior
        return fasta.parseFasta(seq_fp)[0][1]
    return ''


__all__ = [
    'fasta', 'gbtofasta', 'xmfa', 'getSeqFromFile',
    'gb_reader', 'gb_writer',
]
