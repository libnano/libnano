# -*- coding: utf-8 -*-
'''libnano.helpers.jsonbytes
~~~~~~~~~~~~~~~~~~~~~~~~~

Patched json loading to return dictionary with all key and value strings
as bytes objects, rather than unicode.
'''

import json
import sys
from collections import OrderedDict
from typing import (
    List,
    Union
)

def base_decode_list(data: List) -> List:
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
# end def

def base_decode_dict(data: dict) -> dict:
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
    """ OrderedDict with byte strings
    """
    def __init__(self, *args, **kwds):
        super(BOrderedDict, self).__init__(*args, **kwds)

    def __setitem__(self, key, value,
                    dict_setitem=dict.__setitem__):
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
# end class

def load(   fd: '_io.TextIOWrapper',
            cls=None,
            parse_float=None,
            parse_int=None,
            parse_constant=None,
            ordered: bool = False,
            **kw) -> Union[dict, BOrderedDict]:
    return loads(fd.read(),
        cls=cls,
        parse_float=parse_float, parse_int=parse_int,
        parse_constant=parse_constant, ordered=ordered, **kw)
# end def


def loads(s, encoding=None,
            cls=None,
            parse_float=None,
            parse_int=None,
            parse_constant=None,
            ordered: bool = False,
            **kw) -> Union[dict, BOrderedDict]:
    if ordered:
        return json.loads(s,
            cls=cls,
            parse_float=parse_float,
            parse_int=parse_int,
            parse_constant=parse_constant,
            object_pairs_hook=BOrderedDict, **kw)
    else:
        return json.loads(s,
            cls=cls,
            object_hook=base_decode_dict,
            parse_float=parse_float,
            parse_int=parse_int,
            parse_constant=parse_constant, **kw)
# end def
