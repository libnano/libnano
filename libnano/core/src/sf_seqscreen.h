/******************************************************************************
** sf_seqscreen.h
** 
** C functions and respective Python C API functions for filtering DNA or RNA
** sequences based on subsequence content or general sequnce attributes
** (i.e., GC clamps, homopolymer runs, GC content, etc.)
**
******************************************************************************/

// Return 1 if a sequence falls within the maximum allowable values for each
// type of run, otherwise return 0.
int sf_containsRun(char* seq, int seq_length, int maxA, int maxT, int maxG,
              int maxC, int maxAT, int maxGC);

int sf_gcWindow(char* seq, int length, int gc_min_percent, int gc_max_percent,
              int window_size, int* window_start, int* gc_count);    
