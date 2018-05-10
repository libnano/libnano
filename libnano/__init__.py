# -*- coding: utf-8 -*-
"""libnano
    ~~~~~~

    DNA folding and thermodynamics toolkit

"""
__author__: str = "Nick Conway, Ben Pruitt"
__copyright__ = 'Copyright 2018, Nick Conway; Wyss Institute Harvard University'
__license__ = 'GPL2'
__version__: str = '0.2.2.0'

import os
from typing import List

PACKAGE_DIR: str = os.path.dirname(os.path.realpath(__file__))
DATASET_DIR: str = os.path.join(PACKAGE_DIR, 'datasets')


def includes() -> List[str]:
    return [PACKAGE_DIR]


