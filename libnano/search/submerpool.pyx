from seedmatch import SeedMatcher
from libnano.seqstr import thresholdRollingHammingList

# computed using seedfinder.pyx
PRECOMPUTED_SEEDS = {
    (6, 1): ['#-#-#'],
    #(6, 2): ['#--#'],
    (7, 1): ['##-#-'],
    (7, 2): ['##-#', '#--##'],
    (8, 1): ['###-#'],
    #(8, 2): [],
    (9, 1): ['##-##-#'],
    (9, 2): ['###-#-', '#--###'],
    (10, 1): ['######--', '##-#-###-'],
    #(10, 1): ['##-##-##'],
    (10, 2): ['####---', '-#--###'],
    (11, 1): ['###-##--'],
    (11, 2): ['####-#-', '#-#---###'],
    (12, 1): ['###-###--'],
    (12, 2): ['#####---', '#-#--###'],
    (13, 1): ['####-###-'],
    (13, 2): ['#####--#-', '#--#-#-###'],
    (14, 1): ['####-####-'],
    (14, 2): ['#####--##', '##--#-#-###'],
    (15, 1): ['###-###-###-'],
    (15, 2): ['######--#-', '##--#-#-#-##'], #['##-#--##-#']
    (16, 1): ['###-###-###-#'],
    (16, 2): ['#####-#-##-', '#-##--#####'],
    (17, 1): ['###-###-###-##'],
    (17, 2): ['#######-#--', '#-#-#--##--###'],
    (18, 1): ['##-##-##-##-##-#'],
    (18, 2): ['######-#-##-', '#-##--######'], #['###-#--###-#']
    (19, 1): ['###-###-###-###-'],
    (19, 2): ['####-###-#--##', '##-#--####-###']
}

cdef class SubmerPoolSearch:
    """ Uses `SeedMatcher`, `WordMatcher`, or `thresholdRollingHammingList`
    to screen for mismatches of different queries against the list of
    sequences using the `find` method.  This is a convenience class to take
    care of overhead in deciding what matcher to use
    """
    cdef list sequences
    cdef int submer_size
    cdef int mismatches
    cdef list matchers
    cdef list seeds

    def __init__(self, list sequences, int submer_size,
                        int mismatches=0, force_hamming=False):
        """
        Args:
            sequences (List[str]): list of sequences to screen
            submer_size (int): size of the submer to look for
            mismatches (Optional[int]): the number of mismatches to allow
                0 : exact match
                1, 2 : less than a 6 - 19 mer submer size will use a seed if
                    defined
                2 < : will perform a thresholded rolling hamming distance
            force_hamming (Optional[bool]): will perform thresholded rolling
                hamming distance no matter what
        """
        self.sequences = sequences
        self.submer_size = submer_size
        self.mismatches = mismatches

        self.matchers = None
        self.seeds = None

        if force_hamming:
            return
        if mismatches == 0:
            self.matchers = [WordMatcher(sequences, submer_size)]
        # end def
        else:
            seeds = PRECOMPUTED_SEEDS.get((submer_size, mismatches), None)
            if seeds is not None:
                self.matchers = [SeedMatcher(seed, sequences) for seed in seeds]
            self.seeds = seeds
    # end def

    def getMatchers(self):
        """
        Returns:
            The matcher objects to be used if the exist
        """
        return self.matchers

    def find(self, target):
        """
        Args:
            target (str): query string

        Returns:
            List[Tuple[int, int]]: list of tuples where the first value `i`, is
                the index into self.sequences list and the second `j` is the
                index into the self.sequences[i] where the hit occurs.
        """
        if self.matchers is None:
            return thresholdRollingHammingList( target,
                                                self.sequences,
                                                self.mismatches)
        else:
            hits = []
            for matcher in self.matchers:
                hits += matcher.match(target, self.mismatches)
            if len(self.matchers) > 1:
                # dedupe
                return list(set(hits))
            else:
                return hits
    #end def
# end class

cdef class WordMatcher:
    """
    Exact matches, just use python dict
    """
    cdef list sequences
    cdef int submer_size
    cdef dict table

    def __init__(self, list sequences, int submer_size):
        """
        Args:
            sequences (List[str]): list of sequences to screen
            submer_size (int): size of the submer to look for
        """
        cdef Py_ssize_t i, j
        cdef list submer_list
        self.sequences = sequences
        self.submer_size = submer_size
        self.table = {}
        table = self.table
        for i, seq in enumerate(sequences):
            for j in range(len(seq) - submer_size + 1):
                submer = seq[j:j+submer_size]
                submer_list = table.get(submer, None)
                if submer_list is None:
                    submer_list = []
                    table[submer] = submer_list
                submer_list.append((i, j))
            # end for
        # end for
    # end def

    def match(self, target, mismatches=0):
        """
        Args:
            target (str): query string
            mismatches (Optional[int]): does nothing here but exists for API
                compatibility will calls in SubmerPoolSearch

        Returns:
            List[Tuple[int, int]]: list of tuples where the first value `i`, is
                the index into self.sequences list and the second `j` is the
                index into the self.sequences[i] where the hit occurs.
        """
        return self.table.get(target, [])
    # end def
# end class
