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
libnano.search.restrictionmatch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''

from typing import (
    NamedTuple,
    Tuple,
)

CUT_IDXS_T = Tuple[Tuple[int, int], ...]


class RestrictionMatch(NamedTuple):
    strand_dir: int
    start_idx: int
    end_idx: int
    enzyme: str
    cutside_number: int  # Index into the list of cutsites for the enzyme
    cut_idxs: CUT_IDXS_T
    pair_regex: str
