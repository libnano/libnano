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
libnano.fileio.gbtofasta
~~~~~~~~~~~~~~~~~~~~~~~~
Rips the sequence out of a Genbank file and outputs a fasta file:

<filename.gb> -> <filename.fa>

'''
import os
from typing import (
    Dict,
    Optional,
)

from libnano.fileio import gb_reader


def genbank_to_fasta(
        genbank_fp: str,
        seq_record: Optional[Dict] = None,
) -> None:
    '''
    Args:
        genbank_fp: optional Genbank filepath
        seq_record: Optional seq_record dictionary
    '''
    if seq_record:
        primary_rec = seq_record
    else:
        rec = gb_reader.parse(genbank_fp)
        primary_rec = rec

    genbank_fp_no_ext = os.path.splitext(genbank_fp)[0]
    output_fp = genbank_fp_no_ext + '.fa'
    with open(output_fp, 'w', encoding='utf-8') as output_fd:
        output_fd.write('>{}\n'.format(genbank_fp_no_ext.split('/\\')[-1]))
        output_fd.write(str(primary_rec['seq']))
        output_fd.write('\n')


if __name__ == '__main__':
    import sys
    genbank_fp = sys.argv[1]
    genbank_to_fasta(genbank_fp)
