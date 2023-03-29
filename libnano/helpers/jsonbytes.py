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
libnano.helpers.jsonbytes
~~~~~~~~~~~~~~~~~~~~~~~~~

Patched json loading to return dictionary with all key and value strings
as bytes objects, rather than unicode.
'''
import json
# import sys
from collections import OrderedDict
from typing import (
    Any,
    Dict,
    List,
    Tuple,
    Union,
)

import _io  # type: ignore


def base_decode_list(
        data: Union[List[Any], Tuple[Any, ...]],
) -> List[Any]:
    res = []
    for item in data:
        if isinstance(item, str):
            item = item.encode('utf-8')
        elif isinstance(item, (tuple, list)):
            item = base_decode_list(item)
        elif isinstance(item, dict):
            item = base_decode_dict(item)
        res.append(item)
    return res


def base_decode_dict(data: Dict) -> Dict:
    res = {}
    for key, value in data.items():
        if isinstance(key, str):
            key = key.encode('utf-8')
        if isinstance(value, str):
            value = value.encode('utf-8')
        elif isinstance(value, list):
            value = base_decode_list(value)
        elif isinstance(value, dict):
            value = base_decode_dict(value)
        res[key] = value
    return res
# end def


class BOrderedDict(OrderedDict):
    '''``OrderedDict`` with byte strings
    '''

    def __init__(self, *args, **kwds):
        super(BOrderedDict, self).__init__(*args, **kwds)

    def __setitem__(
        self, key, value,
        dict_setitem=dict.__setitem__,
    ):
        if isinstance(key, str):
            key = key.encode('utf-8')
        if isinstance(value, str):
            value = value.encode('utf-8')
        elif isinstance(value, list):
            for i in range(len(value)):
                item = value[i]
                if isinstance(item, str):
                    value[i] = item.encode('utf-8')
        super(BOrderedDict, self).__setitem__(key, value, dict_setitem)


def load(
    fd: _io.TextIOWrapper,
    cls=None,
    parse_float=None,
    parse_int=None,
    parse_constant=None,
    ordered: bool = False,
    **kw,
) -> Union[Dict, BOrderedDict]:
    return loads(
        fd.read(),
        cls=cls,
        parse_float=parse_float, parse_int=parse_int,
        parse_constant=parse_constant, ordered=ordered, **kw,
    )


def loads(
    s,
    encoding=None,
    cls=None,
    parse_float=None,
    parse_int=None,
    parse_constant=None,
    ordered: bool = False,
    **kw,
) -> Union[Dict, BOrderedDict]:
    if ordered:
        return json.loads(
            s,
            cls=cls,
            parse_float=parse_float,
            parse_int=parse_int,
            parse_constant=parse_constant,
            object_pairs_hook=BOrderedDict, **kw,
        )
    else:
        return json.loads(
            s,
            cls=cls,
            object_hook=base_decode_dict,
            parse_float=parse_float,
            parse_int=parse_int,
            parse_constant=parse_constant, **kw,
        )
