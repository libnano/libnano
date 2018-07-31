import sys
from os.path import join, abspath, dirname

root_dir = dirname(dirname(abspath(__file__)))
sys.path = [root_dir] + sys.path
