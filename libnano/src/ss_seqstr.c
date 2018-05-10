/******************************************************************************
** seqstr.c
**
** C functions for general manipulation of DNA, RNA, and amino acid sequences
** (i.e., computing complement / reverse / reverse complement etc).
**
******************************************************************************/

#include <stdio.h>
#include "ss_seqstr.h"

const char COMP_BASE_LUT[128] = {
 '?', '?', '?', '?', '?', '?', '?', '?', //7
 '?', '?','\n', '?', '?', '?', '?', '?', //15
 '?', '?', '?', '?', '?', '?', '?', '?', //23
 '?', '?', '?', '?', '?', '?', '?', '?', //31
 '?', '?', '?', '?', '?', '?', '?', '?', //39
 '?', '?', '?', '?', '?', '?', '?', '?', //47
 '?', '?', '?', '?', '?', '?', '?', '?', //55
 '?', '?', '?', '?', '?', '?', '?', '?', //63
 '?', 'T', 'V', 'G', 'H', '?', '?', 'C', //71
 'D', '?', '?', 'M', '?', 'K', 'N', '?', //79
 '?', '?', 'Y', 'S', 'A', 'A', 'B', 'W', //87
 '?', 'R', '?', '?', '?', '?', '?', '?', //95
 '?', 't', 'v', 'g', 'h', '?', '?', 'c', //103
 'd', '?', '?', 'm', '?', 'k', 'n', '?', //111
 '?', '?', 'y', 's', 'a', 'a', 'b', 'w', //119
 '?', 'r', '?', '?', '?', '?', '?', '?'  //127
};

int ss_revSeq(char* seq, int length){
    char        temp;
    char*       end_ptr;

    end_ptr = seq + (length - 1);
    while (end_ptr > seq) {
        temp = *end_ptr;
        *end_ptr-- = *seq;
        *seq++ = temp;
    }
    return 0;
}

int ss_compSeq(char* seq, int length) {
    int i;
    for (i=0; i < length; i++) {
        *seq = COMP_BASE_LUT[(unsigned char)*seq];
        seq++;
    }
    return 0;
}

int ss_revCompSeq(char* seq, int length) {
    char        temp;
    char*       end_ptr;

    end_ptr = seq + (length - 1);
    while (end_ptr > seq) {
        temp = *end_ptr;
        *end_ptr-- = COMP_BASE_LUT[(unsigned char)*seq];
        *seq++ = COMP_BASE_LUT[(unsigned char)temp];
    }

    if (length % 2) {
        *seq = COMP_BASE_LUT[(unsigned char)*seq];
    }

    return 0;
}

int ss_revCompSeqCopy(char *dest, const char *from, int length) {
    /* does not copy the 0 terminator */
    char *from_ptr = (char *)from + length - 1;
    const char *start_ptr = from - 1;
    char *to_ptr = dest;
    while (from_ptr > start_ptr) {
        *to_ptr++ = COMP_BASE_LUT[(unsigned char) *from_ptr--];
    }
    return 0;
}

int ss_hamming(const char *s1, const char *s2, int len) {
    int hamming_distance = 0;
    char *s1_ptr = (char *) s1;
    char *s2_ptr = (char *) s2;
    const char * s1_ptr_lim = s1_ptr + len;

    while (s1_ptr < s1_ptr_lim) {
        if (*s1_ptr++ != *s2_ptr++) {
            hamming_distance++;
        }
    }
    return hamming_distance;
}

int ss_rollingHamming(const char *s1, const char *s2, 
                        int len1, int len2, 
                        int overlap,
                        int* hamming_arr_ptr) {
    /*
    Assumes len1 < len2

    hamming_arr_ptr (int*): should be size len2 - len1 + 1 + 2*overlap

    Internal
    XXXXXXX          to           XXXXXXX
    YYYYYYYYYYYYY           YYYYYYYYYYYYY
    
    5' overlap
    XXXXXXX
      YYYYYYYYYYYYY

    3' overlap
            XXXXXXX
    YYYYYYYYYYYYY
    */
    const int num_positions = len2 - len1 + 1;  /* internal number of positions */
    int* hamming_arr_ptr_lim = hamming_arr_ptr + num_positions + overlap;
    int* ham_ptr;
    char *s1_ptr;
    char *s2_ptr;
    char *s2_start_ptr = (char *) s2;
    int i;
    // check internal
    for (ham_ptr = hamming_arr_ptr + overlap; 
        ham_ptr < hamming_arr_ptr_lim;
        ham_ptr++) {
        
        s1_ptr = (char *) s1;
        s2_ptr = s2_start_ptr++;
        const char *s2_ptr_lim = s2_ptr + len1;
        while (s2_ptr < s2_ptr_lim) {
            if (*s1_ptr++ != *s2_ptr++) {
                (*ham_ptr)++;
            }
        }
    }
    if (overlap > 0) {
        // check 5' overlap
        i = overlap;
        hamming_arr_ptr_lim = hamming_arr_ptr + overlap;
        for (ham_ptr = 0;
            ham_ptr < hamming_arr_ptr_lim;
            ham_ptr++, i--) {
            s1_ptr = (char *) s1 + i; // offset target
            s2_ptr = (char *) s2;
            const char *s2_ptr_lim = s2_ptr + len1 - i;
            while (s2_ptr < s2_ptr_lim) {
                if (*s1_ptr++ != *s2_ptr++) {
                    (*ham_ptr)++;
                }
            }
        }

        // check 3' overlap
        i = 1;
        hamming_arr_ptr_lim = hamming_arr_ptr + num_positions + 2*overlap;
        for (ham_ptr = hamming_arr_ptr + num_positions + overlap;
            ham_ptr < hamming_arr_ptr_lim;
            ham_ptr++, i++) {
            s1_ptr = (char *) s1; // offset target
            s2_ptr = (char *) s2 + num_positions - len1 + i;
            const char *s2_ptr_lim = s2_ptr + len1 - i;
            while (s2_ptr < s2_ptr_lim) {
                if (*s1_ptr++ != *s2_ptr++) {
                    (*ham_ptr)++;
                }
            }
        }
    }
    return 0;
}

int ss_minHamming(const char *s1, const char *s2, int len1, int len2) {

    const int num_positions = len2 - len1 + 1;
    char *s1_ptr;
    char *s2_ptr;
    char *s2_start_ptr = (char *) s2;
    int min_hd = 100000;
    int current_hd;
    int i=0;

    while (i++ < num_positions) {
        s1_ptr = (char *) s1;
        s2_ptr = s2_start_ptr++;
        current_hd = 0;
        const char *s2_ptr_lim = s2_ptr + len1;
        while (s2_ptr < s2_ptr_lim) {
            if (*s1_ptr++ != *s2_ptr++) {
                current_hd++;
            }
        }
        if (current_hd < min_hd) {
            min_hd = current_hd;
        }
    }
    return min_hd;
}

int ss_minHammingThreshold(const char *s1, const char *s2, 
                                int len1, int len2, int threshold) {

    const int num_positions = len2 - len1 + 1;
    char *s1_ptr;
    char *s2_ptr;
    char *s2_start_ptr = (char *) s2;
    int min_hd = 100000;
    int current_hd;
    int i=0;

    while (i++ < num_positions) {
        s1_ptr = (char *) s1;
        s2_ptr = s2_start_ptr++;
        current_hd = 0;
        const char *s2_ptr_lim = s2_ptr + len1;
        while (s2_ptr < s2_ptr_lim) {
            if (*s1_ptr++ != *s2_ptr++) {
                current_hd++;
            }
        }
        if (current_hd < min_hd) {
            min_hd = current_hd;
            if (threshold > min_hd) {
                goto mht_finish;
            }
        }
    }
    mht_finish:
        return min_hd;
}

#define MISMATCH_LENGTH 3
unsigned char ss_hammingAnd3pMismatchCheck(const char *s1, const char *s2, 
                 const int len, const int min_dist) {
    /* do hamming checks of s1 versus s2 in the following four (4) ways
     1. full s1, full s2
     2. MISMATCH_LENGTH 3p (right most) bases looking for a mismatch of s1, 
        and s2
     3. full s1, full reverse s2 (not the complement)
     4. MISMATCH_LENGTH 3p (right most) bases looking for a mismatch of s1, 
        and rev s2
     5. MISMATCH_LENGTH 3p (right most) bases looking for a mismatch of rev s1, 
        and s2.  This is needed for symmetry
    */
    int hamming_distance = 0;
    int hamming_distance_rev = 0;
    int mismatch_distance = 0;
    int mismatch_distance_revA = 0;
    int mismatch_distance_revB = 0;

    char *s1_ptr = (char *) s1;
    char *s2_ptr = (char *) s2;
    const char * s1_ptr_lim = s1_ptr + len;
    const char * s1_ptr_lim_3p = s1_ptr_lim - MISMATCH_LENGTH;
    const char * s1_ptr_lim_5p = s1_ptr + MISMATCH_LENGTH;

    // note:  due to operator precedence 
    // *++x increments and then derefereces, *x++ dereferences x and 
    // increments x after using *x

    // 1. get the hamming distance
    while (s1_ptr < s1_ptr_lim) {
        if (*s1_ptr++ != *s2_ptr++) {
            hamming_distance++;
        }
    }
    // 2. reverse the pointers to get the 3p mismatch
    while (s1_ptr > s1_ptr_lim_3p) {
        if (*--s1_ptr != *--s2_ptr) {
            mismatch_distance++;
        }
    }

    // 3. now do the reverse full hamming incrementing s2_ptr from the end
    s1_ptr = (char *) s1;
    s2_ptr = ((char *) s2) + len;
    while (s1_ptr < s1_ptr_lim) {
        if (*s1_ptr++ != *--s2_ptr) {
            hamming_distance_rev++;
        }
    }

    // 4. reverse the pointers to get the 3p mismatch of the reverse of s2
    while (s1_ptr > s1_ptr_lim_3p) {
        if (*--s1_ptr != *s2_ptr++) {
            mismatch_distance_revA++;
        }
    }

    // 5. reverse the pointers to get the 3p mismatch with the reverse of s1
    s1_ptr = (char *) s1;
    s2_ptr = ((char *) s2) + len;
    while (s1_ptr < s1_ptr_lim_5p) {
        if (*s1_ptr++ != *--s2_ptr) {
            mismatch_distance_revB++;
        }
    }

    
    // now check the

    // check this out.
    if ((hamming_distance > min_dist) && 
        (mismatch_distance) && 
        (hamming_distance_rev > min_dist) && 
        (mismatch_distance_revA) &&
        (mismatch_distance_revB)
        ) {
        return 1;
    } else {
        return 0;
    }
}

int ss_hamming2XCheck(const char *s1, const char *s2, int len, const int min_dist) {
    int hamming_distance = 0;
    int hamming_distance_rev = 0;
    
    char *s1_ptr = (char *) s1;
    char *s2_ptr = (char *) s2;
    const char * s1_ptr_lim = s1_ptr + len;

    while (s1_ptr < s1_ptr_lim) {
        if (*s1_ptr++ != *s2_ptr++) {
            hamming_distance++;
        }
    }

    // 3. now do the reverse full hamming incrementing s2_ptr from the end
    s1_ptr = (char *) s1;
    s2_ptr = ((char *) s2) + len;
    while (s1_ptr < s1_ptr_lim) {
        if (*s1_ptr++ != *--s2_ptr) {
            hamming_distance_rev++;
        }
    }

    // check this out.
    if ((hamming_distance > min_dist) && 
        (hamming_distance_rev > min_dist)
        ) {
        return 1;
    } else {
        return 0;
    }
}