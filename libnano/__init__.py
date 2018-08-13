# -*- coding: utf-8 -*-
"""libnano
    ~~~~~~

    DNA folding and thermodynamics toolkit

"""
__author__: str = "Nick Conway, Ben Pruitt"
__copyright__ = 'Copyright 2018, Nick Conway; Wyss Institute Harvard University'
__license__ = 'GPL2'
__version__: str = '0.2.3.0'

from typing import List

from os.path import (
    dirname,
    realpath,
    join
)
import sys

PACKAGE_DIR: str = dirname(realpath(__file__))
DATASET_DIR: str = join(PACKAGE_DIR, 'datasets')

sys.path = [dirname(PACKAGE_DIR)] + sys.path

def includes() -> List[str]:
    return [PACKAGE_DIR]




