/******************************************************************************
** sf_seqmetric.h
**
** C functions for calculating basic DNA sequence metrics
**
******************************************************************************/

// Calculate the max A/T/G/C/AT/GC runs in `seq` and put return values in
// `ret_arr` ([maxA, maxC, maxG, maxT, maxAT, maxGC])
void sm_maxRuns(char* seq, int seq_length, int* ret_arr);

// Calculate the global GC content
float sm_gcContent(char* seq, int seq_length);
