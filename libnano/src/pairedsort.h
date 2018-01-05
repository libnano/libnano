/* Sort by an array of type_t while also sorting and array of slave_t
* this is done with the expectation that the array of type_t will be freed
* soon after this sorting is complete rather than building a struct to couple 
* them
* this is forked from ksort.h in klib
*/
#define PAIREDSORT_INIT(name, type_t, slave_t, __sort_lt) \
void mergesort_paired_##name(size_t n, type_t* keys, slave_t* idxs) { \
    type_t*sort_group_k[2], **a_k, **b_k; \
    slave_t *sort_group_i[2], *a_i, *b_i; \
    int curr, shift; \
\
    sort_group_k[0] = keys; \
    sort_group_k[1] = (type_t*) malloc(n*sizeof(type_t)); \
\
    sort_group_i[0] = idxs; \
    sort_group_i[1] = (slave_t *) malloc(n*sizeof(slave_t)); \
\
    for (curr = 0, shift = 0; (1ul << shift) < n; ++shift) { \
        a_k = sort_group_k[curr]; b_k = sort_group_k[1 - curr]; \
        a_i = sort_group_i[curr]; b_i = sort_group_i[1 - curr]; \
        if (shift == 0) { \
\
            type_t *p_k = b_k, **i_k, **eb_k = a_k + n; \
            slave_t *p_i = b_i, *i_i; \
\
            for (i_k = a_k, i_i = a_i; \
                    i_k < eb_k; \
                    i_k += 2, i_i +=2) { \
                if (i_k == eb_k - 1) { \
                    *p_k++ = *i_k; \
                    *p_i++ = *i_i; \
                } else { \
                    if (__sort_lt(*(i_k + 1), *i_k)) { \
                        *p_k++ = *(i_k + 1); \
                        *p_k++ = *i_k; \
                        *p_i++ = *(i_i + 1); \
                        *p_i++ = *i_i; \
                    } else { \
                        *p_k++ = *i_k; \
                        *p_k++ = *(i_k+1); \
                        *p_i++ = *i_i; \
                        *p_i++ = *(i_i+1); \
                    } \
                } \
            } \
        } else { \
            size_t i, step = 1ul << shift; \
            for (i = 0; i < n; i += step << 1) { \
                type_t *p_k, **j_k, **k_k, **ea_k, **eb_k; \
                slave_t *p_i, *j_i, *k_i; \
\
                if (n < i + step) { \
                    ea_k = a_k + n; eb_k = a_k; \
                } else { \
                    ea_k = a_k + i + step; \
                    eb_k = a_k + (n < i + (step << 1) ? n : i + (step << 1)); \
                } \
\
                j_k = a_k + i; \
                k_k = a_k + i + step; \
                p_k = b_k + i; \
\
                j_i = a_i + i; \
                k_i = a_i + i + step; \
                p_i = b_i + i; \
\
                while ( (j_k < ea_k) && (k_k < eb_k) ) { \
                    if (__sort_lt(*k_k, *j_k)) { \
                        *p_k++ = *k_k++; \
                        *p_i++ = *k_i++; \
                    } else { \
                        *p_k++ = *j_k++; \
                        *p_i++ = *j_i++; \
                    } \
                } \
                while (j_k < ea_k) { \
                    *p_k++ = *j_k++; \
                    *p_i++ = *j_i++; \
                } \
                while (k_k < eb_k) { \
                    *p_k++ = *k_k++; \
                    *p_i++ = *k_i++; \
                } \
            } \
        } \
        curr = 1 - curr; \
    } \
    if (curr == 1) { \
        type_t *p_k = sort_group_k[0], **i_k = sort_group_k[1], **eb_k = keys + n; \
        slave_t *p_i = sort_group_i[0], *i_i = sort_group_i[1]; \
\
        for (; p_k < eb_k; ++i_k, ++i_i) { \
            *p_k++ = *i_k; \
            *p_i++ = *i_i; \
        } \
    } \
    free(sort_group_k[1]); \
    free(sort_group_i[1]); \
}
