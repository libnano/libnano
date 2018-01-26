"""
    libnano
    ~~~~~~

    DNA folding and thermodynamics toolkit

"""
import os
LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))

def includes():
    return [LOCAL_DIR]

__author__ = "Nick Conway, Ben Pruitt"
__all__ = ['seqstr', 'seqint', 'seqscreen']
__version__ = '0.1.1.3'
