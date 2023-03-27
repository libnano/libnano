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
tests.test_geneprober
~~~~~~~~~~~~~~~~~~~~~

'''
import pytest

from libnano.scripts import geneprober


def test_cli():
    with pytest.raises(SystemExit):
        # Need this because click calls sys.exit(0) and pytest fails for
        # that
        geneprober.cli()


def test_list_details():
    with geneprober.disable_cache():
        geneprober.listDetails(
            'mouse',
            ['DAXX'],
        )


def test_list_detailsfail():
    a_fake_gene_symbol: str = 'DOOF'
    with geneprober.disable_cache():
        with pytest.raises(KeyError):
            geneprober.listDetails(
                'mouse',
                [a_fake_gene_symbol],
            )
