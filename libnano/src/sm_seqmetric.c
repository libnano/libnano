/******************************************************************************
** sf_seqmetric.c
**
** C functions for calculating basic DNA sequence metrics
**
******************************************************************************/


void sm_maxRuns(char* seq, int seq_length, int* ret_arr){
    // ret_arr = [maxA, maxC, maxG, maxT, maxAT, maxGC]
    int     Arun = 0, Trun = 0, Grun = 0, Crun = 0, ATrun = 0, GCrun = 0;
    int     i = 0;
    char    *seq_ptr = seq;

    for (i = 0; i < seq_length; i++, seq_ptr++) {
        if (*seq_ptr == 'A') {
            Arun++;
            ATrun++;
            if (Arun > ret_arr[0]) {ret_arr[0] = Arun;}
            if (ATrun > ret_arr[4]) {ret_arr[4] = ATrun;}
            Trun = Grun = Crun = GCrun = 0;
        } else if (*seq_ptr == 'T') {
            Trun++;
            ATrun++;
            if (Trun > ret_arr[3]) {ret_arr[3] = Trun;}
            if (ATrun > ret_arr[4]) {ret_arr[4] = ATrun;}
            Arun = Grun = Crun = GCrun = 0;
        } else if (*seq_ptr == 'G') {
            Grun++;
            GCrun++;
            if (Grun > ret_arr[2]) {ret_arr[2] = Grun;}
            if (GCrun > ret_arr[5]) {ret_arr[5] = GCrun;}
            Crun = Arun = Trun = ATrun = 0;
        } else if (*seq_ptr == 'C') {
            Crun++;
            GCrun++;
            if (Crun > ret_arr[1]) {ret_arr[1] = Crun;}
            if (GCrun > ret_arr[5]) {ret_arr[5] = GCrun;}
            Grun = Arun = Trun = ATrun = 0;
        }
    }
}

float sm_gcContent(char* seq, int seq_length) {
    int i, gc_count = 0;
    char *seq_ptr = seq;

    for (i = 0; i < seq_length; i++, seq_ptr++) {
        gc_count += (*seq_ptr == 'G' || *seq_ptr == 'C') ? 1:0;
    }
    return (1.0 * gc_count) / seq_length;
}
