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
libnano.seqrecord.featuretypes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''
from bisect import insort_left  # bisect_left,
from typing import (
    Dict,
    List,
    Optional,
    Tuple,
)


class FeatureTypes(object):
    '''Data structure to enumerate and track feature types

    Attributes:
        count: Number of features registered
        ftype_dict: Map of Feature name to Feature ID
        ftype_descriptions: Map of Feature name to Feature Descriptions
        ftype_LUT: Lookup table keyed by assigned enumeration per type
        recycle_bin: Recycle bin for Feature ID re-use
    '''

    def __init__(self):
        self.count = 0
        self.ftype_dict: Dict[str, int] = {}
        self.ftype_descriptions: Dict[str, str] = {}
        self.ftype_LUT: List[str] = []
        self.recycle_bin: List[int] = []

    def __contains__(self, feature_name: str) -> bool:
        '''Allow for ``feature_name in feature_types`` construction

        Args:
            feature_name: Name of the feature
        '''
        return feature_name in self.ftype_dict

    def getFTID(
            self,
            feature_name: str,
    ) -> int:
        '''Get the feature ID by name

        Args:
            feature_name: Name of the feature

        Returns:
            Feature ID

        '''
        return self.ftype_dict[feature_name]

    def add_feature_type(
        self,
        feature_name: str,
        description: str = '',
        ft_id: Optional[int] = None,
    ) -> int:
        '''Add a feature type

        Args:
            feature_name: Name of the feature
            description: Feature type description
            ft_id: If not ``None`` is used in a redo kind of operation

        Returns:
            Feature ID

        '''
        ftype_LUT = self.ftype_LUT
        ftype_dict = self.ftype_dict
        rb = self.recycle_bin
        if feature_name in ftype_dict:
            return ftype_dict[feature_name]
        else:
            if ft_id is not None:
                try:
                    existing_fname = ftype_LUT[ft_id]
                    if existing_fname is None:
                        ftype_LUT[ft_id] = feature_name
                    elif existing_fname == feature_name:
                        raise ValueError(
                            f'add_feature_type: feature_name {feature_name} '
                            'already exists',
                        )
                    else:
                        raise ValueError(
                            f'add_feature_type: feature_name {feature_name} '
                            f'with ft_id {ft_id} collision',
                        )
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

    def remove_feature_type(
            self,
            feature_name: str,
    ) -> Tuple[int, str]:
        '''
        Args:
            feature_name: Name of the feature

        Returns:
            Tuple of the form::

            <Feature ID>, <feature type description>
        '''
        ftype_dict = self.ftype_dict
        desc = self.ftype_descriptions[feature_name]
        del self.ftype_descriptions[feature_name]
        ft_id = ftype_dict[feature_name]
        self.ftype_LUT[ft_id] = ''
        del ftype_dict[feature_name]
        insort_left(self.recycle_bin, ft_id)
        self.count -= 1
        return ft_id, desc
