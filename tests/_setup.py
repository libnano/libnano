import sys
from os.path import join, abspath, dirname

import numpy

print(numpy.__file__, numpy.__version__)

sys.path = [abspath(join(dirname(__file__), '..'))] + sys.path

print(sys.path)
from libnano.core import seqstr
print(dir(seqstr))
