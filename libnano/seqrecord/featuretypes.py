# -*- coding: utf-8 -*-
from bisect import (
    bisect_left,
    insort_left
)
from typing import (
    List,
    Tuple,
    Dict
)

class FeatureTypes(object):
    """Data structure to enumerate and track feature types
    """
    def __init__(self):
        self.count: int = 0  # number of features registered
        self.ftype_dict: Dict[str, int] = {}
        self.ftype_descriptions: Dict[str, str] = {}
        self.ftype_LUT: List[str] = [] # lookup table keyed by assigned enumeration per type
        self.recycle_bin: List[int] = []

    def __contains__(self, feature_name: str) -> bool:
        """allow for ``feature_name in feature_types`` construction
        """
        return feature_name in self.ftype_dict

    def getFTID(self, feature_name: str) -> int:
        return self.ftype_dict[feature_name]

    def addFeatureType(self, feature_name: str,
                            description: str = None,
                            ft_id: int = None) -> int:
        """
        ft_id not None is used in a redo kind of operation
        """
        ftype_LUT = self.ftype_LUT
        ftype_dict = self.ftype_dict
        rb = self.recycle_bin
        if feature_name in ftype_dict:
            return
        else:
            if ft_id is not None:
                try:
                    existing_fname = ftype_LUT[ft_id]
                    if existing_fname == None:
                        ftype_LUT[ft_id] = feature_name
                    elif existing_fname == feature_name:
                        msg = "addFeatureType: feature_name" +\
                                "%s already exists" % (feature_name)
                        raise ValueError(msg)
                    else:
                        msg = "addFeatureType: feature_name" +\
                                "%s with ft_id %d collision" % (feature_name, ft_id)
                        raise ValueError(msg)
                except IndexError:
                    raise
            else:
                if len(rb):
                    ft_id = rb.pop(0)
                    ftype_LUT[ft_id] = feature_name
                else:
                    ft_id = len(ftype_LUT)
                    ftype_LUT.append(feature_name)

                ftype_dict[feature_name] = ft_id
                self.ftype_descriptions[feature_name] = description

            self.count += 1
            return ft_id
    # end def

    def removeFeatureType(self, feature_name: str) -> Tuple[int , str]:
        ftype_dict = self.ftype_dict
        desc = self.ftype_descriptions[feature_name]
        del self.ftype_descriptions[feature_name]
        ft_id = ftype_dict[feature_name]
        self.ftype_LUT[ft_id] = None
        del ftype_dict[feature_name]
        insort_left(self.recycle_bin, ft_id)
        self.count -= 1
        return ft_id, desc
    # end def
# end class