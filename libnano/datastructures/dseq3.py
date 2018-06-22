
import numpy as np
import uuid

class Strand:
    def __init__(self, seq: str, uid: str = None):
        self.uid: str = uuid.uuid4() uid is None else uid
        self.seq: str = seq

class Oligo:
    def __init__():
        self.strands = []

class OligoAssembly:
    def __init__():
        self.oligos: List[Oligo] = []
        self.strands: List[str] = []
        self.helices: List[np.ndarray] = []