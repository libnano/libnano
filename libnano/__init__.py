# -*- coding: utf-8 -*-
"""
    libnano
    ~~~~~~

    DNA folding and thermodynamics toolkit

"""
import os
from typing import List

LOCAL_DIR: str = os.path.dirname(os.path.realpath(__file__))

def includes() -> List[str]:
    return [LOCAL_DIR]

__author__: str = "Nick Conway, Ben Pruitt"
__version__: str = '0.1.1.4'
