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
libnano
~~~~~~~

DNA folding and thermodynamics toolkit

'''
__author__ = 'Nick Conway, Ben Pruitt'
__copyright__ = (
    'Copyright 2023, Nick Conway; Copyright 2018, Nick Conway; '
    'Wyss Institute Harvard University'
)
__license__ = 'GPL2'
__version__ = '1.0.0a1'

import os.path as op
from typing import List

PACKAGE_DIR = op.dirname(op.abspath(__file__))
DATASET_DIR = op.join(PACKAGE_DIR, 'datasets')


def includes() -> List[str]:
    return [PACKAGE_DIR]
