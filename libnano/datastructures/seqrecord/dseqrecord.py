import sys
from os.path import abspath, dirname
LIBNANO_PATH = abspath(dirname(dirname(dirname(dirname(__file__)))))
print(LIBNANO_PATH)
sys.path = [LIBNANO_PATH] + sys.path

from libnano.datastructures.seqrecord.seqrecordbase import (
    SeqRecord,
    AlphaEnum,
    ALPHABETS
)
from libnano.seqstr import (
    reverseComplement,
    complement
)
import ssw
import math

def align_complement(fwd: str, rev: str):
    rc_rev = reverseComplement(rev)
    alignment = ssw.force_align(rc_rev, fwd)
    # print(alignment)
    print_force_align(reverse(rev), fwd, alignment)
# end def

class DSeq():
    def __init(self,
        fwd: str,
        rev: str = None,
        overhang: int = None,
        alphabet: int = AlphaEnum.DNA):

        assert(alphabet in ALPHABETS)
        self.alphabet: int = alphabet

        # self.overhang = 0 if overhang is None else overhang

        if rev is None:
            if overhang is not None:
                raise ValueError("overhang can't be defined for without a reverse strand")
            else:
                self.rev = reverseComplement(fwd)
                self.overhang = 0
        else:
            self.rev = rev

            # if overhang is None:
            #     self.

if __name__ == '__main__':
    print("dope")
    align_complement("GGATCCAAA", "TTTGGATC")
    align_complement("GGATCCAAA", "TTGGATC")
    align_complement("GGATCCAAA", "CTTGGATC")