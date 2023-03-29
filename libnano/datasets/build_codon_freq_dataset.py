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
libnano.datasets.build_condon_freq_dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Includes methods for retrieving and parsing codon frequency data from the
Kazusa DNA Research Institute codon usage database.

By default the following species' codon frequencies are included:

    Organism                Class   Organism ID
    --------------------    -----   ----------------------
    E. coli (K12)           bct     83333
    H. sapiens              pri     93487
    S. cerevisiae           pln     14411

Nucleic acid sequences are in RNA space

Format: [<lookup name>, <organism class>, <organism id>]
Lookup name is for user-facing references, not official database look-ups

To find these numbers, look up your organism of interest in the
Kazusa codon database (generally works to Google "<organism name> codon
usage table")

Once you get to the respective page, the ID will be in the url:
http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4922
Here, "4922" is the ID for P. pastoris.

The first line of the page will include the organism class:
Pichia pastoris [gbpln]: 137 CDS's (81301 codons)
Here, "gbpln" is the class, which is shortened to "pln" for lookups
'''

import gzip
import io
import json
import os.path as op
import re
from ftplib import FTP
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
)

LOCAL_DIR: str = op.dirname(
    op.realpath(__file__),
)

_FTP_URLS: Dict[str, Dict[str, str]] = {
    'EBI': {
        'root': 'ftp.ebi.ac.uk',
        'db_dir': 'pub/databases/cutg/',
    },
    'KAZ': {
        'root': 'tp.kazusa.or.jp',
        'db_dir': 'pub/codon/current/',
    },
}


AA_TO_RNA: Dict[str, List[str]] = {
    'F': ['UUU', 'UUC'],
    'S': ['UCU', 'UCA', 'AGU', 'UCG', 'AGC', 'UCC'],
    'Y': ['UAU', 'UAC'],
    'C': ['UGU', 'UGC'],
    'L': ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'H': ['CAU', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'I': ['AUU', 'AUC', 'AUA'],
    'M': ['AUG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'N': ['AAU', 'AAC'],
    'K': ['AAA', 'AAG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'W': ['UGG'],
    '*': ['UAA', 'UGA', 'UAG'],
}

RNA_TO_AA: Dict[str, str] = {v: k for k, vs in AA_TO_RNA.items() for v in vs}
DNA_TO_AA: Dict[str, str] = {
    k.replace('U', 'T'): v for k, v in RNA_TO_AA.items()
}
AA_TO_DNA: Dict[str, List[str]] = {
    k: [v.replace('U', 'T') for v in vs] for k, vs in AA_TO_RNA.items()
}

_ORGANISM_CLASSES: List[str] = [
    'bct',      # bacteria
    'inv',      # invertebrates
    'mam',      # mammals
    'pln',
    'pri',
    'rod',
    'vrl',
    'vrt',
]


ORG_T = Tuple[str, str, int]
_ORGANISM_DEFAULTS: List[ORG_T] = [
    ('E. coli', 'bct', 83333),
    ('H. sapiens', 'pri', 93487),
    ('S. cerevisiae', 'pln', 14411),
    ('P. pastoris', 'pln', 4922),
]

_FREQ_ORDER = None


def _get_ftp_conn(target: str) -> FTP:
    for institute in ('EBI', 'KAZ'):
        try:
            root_url = _FTP_URLS[institute]['root']
            target_path = _FTP_URLS[institute][target]
            ftp = FTP(root_url)
            ftp.login()
            ftp.cwd(target_path)
            return ftp
        except IOError:
            pass
    raise IOError('Could not connect to the EBI or Kazusa FTP databases')


def _get_gzip_data(
        ftp_conn: FTP,
        gz_fn: str,
) -> str:
    buff = io.BytesIO()
    ftp_conn.retrbinary(
        f'RETR {gz_fn}',
        callback=buff.write,
    )
    buff.flush()
    buff.seek(0)
    gz_fh = gzip.GzipFile(
        fileobj=buff,
    )
    data = gz_fh.read().decode('utf-8')
    ftp_conn.close()
    return data


def _get_organism_record_list(
        raw_data: str,
        organism_id: int,
) -> List[int]:
    org_re = f'{organism_id}[^\n]*\n([^\n]*)\n'
    m = re.search(
        org_re,
        raw_data,
    )
    if m is not None:
        return [
            int(g)
            for g in m.group(1).split(' ') if len(g)
        ]
    return []


def _get_frequency_order() -> List[str]:
    global _FREQ_ORDER
    if _FREQ_ORDER is None:
        ftp_conn = _get_ftp_conn('db_dir')
        raw_data = _get_gzip_data(
            ftp_conn,
            'SPSUM_LABEL.gz',
        )
        m = re.search(
            '\n([^\n]+)\n',
            raw_data,
        )
        if m is not None:
            # Ensure codons are in DNA-space rather than RNA-space
            _FREQ_ORDER = [
                c.replace('U', 'T')
                for c in m.group(1).split(' ')
            ]
        else:
            _FREQ_ORDER = []
        ftp_conn.close()
    return _FREQ_ORDER


def _compute_frequencies_dict(
        org_counts: List[int],
) -> Dict[str, Dict[str, float]]:
    '''
    Args:
        org_counts:

    Returns:
        Dictionary mapping amino acids to codons to frequency count values
    '''
    codon_counts = zip(
        _get_frequency_order(),
        org_counts,
    )
    freq_dict: Dict[str, Dict[str, float]] = {}

    for codon, count in codon_counts:
        freq_dict.setdefault(DNA_TO_AA[codon], {})
        freq_dict[DNA_TO_AA[codon]][codon] = count

    for aa, codon_dict in freq_dict.items():
        total_count = sum(codon_dict.values())
        for codon in freq_dict[aa].keys():
            freq_dict[aa][codon] = float(freq_dict[aa][codon]) / total_count
    return freq_dict


def get_organism_frequencies_dict(
        organism_id: int,
        organism_class: str,
) -> Dict[str, Dict[str, float]]:
    '''
    Args:
        organism_id: Organism ID
        organism_class: Organism class

    Returns:
        Dictionary mapping amino acids to codons to frequency count values
    '''
    ftp_conn = _get_ftp_conn('db_dir')
    class_gz_fn = 'gb%s.spsum.gz' % organism_class

    raw_data = _get_gzip_data(
        ftp_conn,
        class_gz_fn,
    )
    org_rec_list = _get_organism_record_list(
        raw_data,
        organism_id,
    )
    return _compute_frequencies_dict(org_rec_list)


def update_codon_frequency_dataset(
        additional_organisms: Optional[List[ORG_T]] = None,
) -> None:
    '''Update the codon frequency dataset (stored in codon_freq_dataset.json)

    Includes the organisms listed in `_ORGANISM_DEFAULTS` in this module.
    Additional organisms may be added by passing a list of lists (in the same
    format as `_ORGANISM_DEFAULTS`) as the `additional_organisms` arg.

    Args:
        additional_organisms: If not None, additional data to search for
    '''
    organism_data_list: List[ORG_T] = additional_organisms or []
    organism_data_list += _ORGANISM_DEFAULTS
    codon_freq_dict = {}
    for org_name, org_class, org_id in organism_data_list:
        codon_freq_dict[org_name] = get_organism_frequencies_dict(
            org_id,
            org_class,
        )
    filepath = op.join(LOCAL_DIR, 'codon_freq_dataset.json')
    with open(filepath, 'w') as fh:
        json.dump(
            codon_freq_dict,
            fh,
        )
