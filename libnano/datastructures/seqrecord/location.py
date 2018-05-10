from collections import namedtuple
from typing import (
    List,
    Tuple
)

"""
compound_location = (operator, (remote, idx5p, idx3p, flags))
flags = <IS_PARTIAL_5 | IS_PARTIAL_3 | IS_INTERNAL]
"""
IS_INTERNAL = 1
IS_PARTIAL_3 = 2
IS_PARTIAL_5 = 4

# base 'simple' location
Location = namedtuple('Location', ['remote', 'idx5p', 'idx3p', 'flags'])
Operator = namedtuple('Operator', ['op', 'arg'])

def createLocationStr(location: Location) -> str:
    if location[0] == 'join':
        return 'join(' + createLocationStr(location[1]) + ')'
    elif location[0] == 'order':
        return 'order(' + createLocationStr(location[1]) + ')'
    elif location[0] == 'complement':
        return 'complement(' + createLocationStr(location[1]) + ')'
    elif isinstance(location, list):
        out_list = []
        for item in location:
            out_list.append(createLocationStr(item))
        return ','.join(out_list)
    else:
        return createLocationStrBases(location)
# end def

def createLocationStrBases(loc: Location) -> str:
    """
    location = (operator, (remote, idx5p, idx3p, flags))
    flags = <IS_PARTIAL_5 | IS_PARTIAL_3 | IS_INTERNAL]
    """
    out_list = []
    if loc.remote is not None:
        out_list += [loc.remote, ':']
    if loc.flags & IS_PARTIAL_5:
        out_list.append('<')
    out_list.append(str(loc.idx5p + 1))
    if loc.idx3p != loc.idx5p: #is not None:
        if loc.flags & IS_INTERNAL:
            out_list.append('^')
        else:
            out_list.append('..')
            if loc.flags & IS_PARTIAL_3:
                out_list.append('>')
        out_list.append(str(loc.idx3p + 1))
    return ''.join(out_list)
# end def

def parseLocation(loc_str: str) -> Tuple:
    """
    assumes last character is a matching parenthesis
    """
    if loc_str[0:5] == 'join(':
        return ('join', parseJoinOrder(loc_str[5:-1]))
    elif loc_str[0:6] == 'order(':
        return ('order', parseJoinOrder(loc_str[6:-1]))
    elif loc_str[0:11] == 'complement(':
        return ('complement', parseLocation(loc_str[11:-1]))
    else:
        return parseLocationBases(loc_str)
#end def

def parseJoinOrder(loc_str: str) -> List[Tuple]:
    # find outer most commas
    comma_idx_list = []
    ignore_flag = 0
    for i, c in enumerate(loc_str):
        if c == ',' and ignore_flag == 0:
            comma_idx_list.append(i)
        elif c == '(':
            ignore_flag += 1
        elif c == ')':
            ignore_flag -= 1
    j = 0
    out_list = []
    for i in comma_idx_list:
        out_list.append(parseLocation(loc_str[j:i]))
        j = i + 1
    out_list.append(parseLocation(loc_str[j:]))
    return out_list
#end def

def parseLocationBases(bases_str: str) -> Location:
    # find remote if it exists
    if ':' in bases_str:
        bases_list = split(':')
        remote = bases_list[0]
        bases_str = bases_list[1]
    else:
        remote = None

    if '^' in bases_str:
        bl = bases_str.split('^')
        return Location(remote, int(bl[0])-1, int(bl[1])-1, IS_INTERNAL)
    elif '..' in bases_str:
        flags = 0
        bl = bases_str.split('..')
        if bl[0][0] == '<':
            bl[0] = bl[0][1:]
            flags += IS_PARTIAL_5
        if bl[1][0] == '<':
            bl[1] = bl[1][1:]
            flags += IS_PARTIAL_3
        return Location(remote, int(bl[0]) - 1, int(bl[1]) - 1, flags)
    else: # single base
        base_idx = int(bases_str) - 1
        return Location(remote, base_idx, base_idx, 0)
# end def