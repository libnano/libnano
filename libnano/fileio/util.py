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
libnano.fileio.util
~~~~~~~~~~~~~~~~~~~

'''
import mmap
import os
from typing import (
    Tuple,
    Union,
)


def get_file_contents(
    fp: str,
    max_in_mem_size: int = 200,
) -> Tuple[Union[bytes, mmap.mmap], bool]:
    '''Return memory map of file

    Args:
        fp: filepath to read
        max_in_mem_size: maximum in memory size before a memory map is returned

    Returns:
        The contents of file at `fp` either as a string (in memory) or
        as a memory-mapped file (if the file size is > `max_in_mem_size` in MB)
    '''
    contents_object: Union[bytes, mmap.mmap] = b''
    is_mmap: bool = False
    with open(fp, 'r+b') as fd:  # r+ read mode necessary for mmap
        if os.path.getsize(fp) / 1024 > max_in_mem_size:
            contents_object = mmap.mmap(fd.fileno(), 0)  # 0->mmap entire file
            is_mmap = True
        else:
            contents_object = fd.read()
    return contents_object, is_mmap
