/******************************************************************************
** seqscreen.c
**
** C functions and respective Python C API functions for filtering DNA or RNA
** sequences based on subsequence content or general sequnce attributes
** (i.e., GC clamps, homopolymer runs, GC content, etc.)
**
******************************************************************************/

int sf_containsRun(char* seq, int seq_length, int maxA, int maxT, int maxG,
              int maxC, int maxAT, int maxGC){
    int     Arun = 0, Trun = 0, Grun = 0, Crun = 0, ATrun = 0, GCrun = 0;
    int     i = 0;
    char    *seq_ptr = seq;

    for (i = 0; i < seq_length; i++, seq_ptr++) {
        if (*seq_ptr == 'A') {
            Arun++;
            ATrun++;
            Trun = Grun = Crun = GCrun = 0;
        } else if (*seq_ptr == 'T') {
            Trun++;
            ATrun++;
            Arun = Grun = Crun = GCrun = 0;
        } else if (*seq_ptr == 'G') {
            Grun++;
            GCrun++;
            Crun = Arun = Trun = ATrun = 0;
        } else if (*seq_ptr == 'C') {
            Crun++;
            GCrun++;
            Grun = Arun = Trun = ATrun = 0;
        }
        if ((GCrun > maxGC) || (ATrun > maxAT) || (Grun > maxG) || \
            (Crun > maxC) || (Trun > maxT) || (Arun > maxA))
            return 0;
    }
    return 1;
}

int sf_gcWindow(char* seq, int length, int gc_min_percent, int gc_max_percent,
              int window_size, int* window_start, int* gc_count) {
    int         i, gc_min, gc_max;
    char        *iptr, *mptr, *eptr;

    gc_min = (gc_min_percent * length) / 100;
    gc_max = (gc_max_percent * length) / 100;
    if (window_size > length) {
        return -1;
    }
    // init window
    *window_start = 0;
    *gc_count = 0;
    iptr = seq;
    mptr = seq;
    eptr = (seq + length);
    for (i = 0; i < window_size; i++) {
        *gc_count += (*mptr == 'G' || *mptr == 'C') ? 1:0;
        mptr++;
    }
    if (*gc_count < gc_min || *gc_count > gc_max) {
        return 0;
    }
    // scan seq
    while (mptr != eptr) {
        *gc_count -= (*iptr == 'G' || *iptr == 'C') ? 1:0;
        *gc_count += (*mptr == 'G' || *mptr == 'C') ? 1:0;
        ++*window_start;
        mptr++;
        iptr++;
        if (*gc_count < gc_min || *gc_count > gc_max) {
            return 0;
        }
    }
    return 1;
}

