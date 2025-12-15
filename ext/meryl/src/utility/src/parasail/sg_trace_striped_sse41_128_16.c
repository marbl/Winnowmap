/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#if HAVE_INTEL

#include <stdint.h>
#include <stdlib.h>

#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <emmintrin.h>
#include <smmintrin.h>
#endif

#include "parasail.h"
#include "memory.h"
#include "internal_sse.h"

#define SG_TRACE
#define SG_SUFFIX _striped_sse41_128_16
#define SG_SUFFIX_PROF _striped_profile_sse41_128_16
#include "sg_helper.h"

#define SWAP(A,B) { __m128i* tmp = A; A = B; B = tmp; }



static inline void arr_store(
        __m128i *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    _mm_store_si128(array + (1LL*d*seglen+t), vH);
}

static inline __m128i arr_load(
        __m128i *array,
        int32_t t,
        int32_t seglen,
        int32_t d)
{
    return _mm_load_si128(array + (1LL*d*seglen+t));
}

#define FNAME parasail_sg_flags_trace_striped_sse41_128_16
#define PNAME parasail_sg_flags_trace_striped_profile_sse41_128_16

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
{
    /* declare local variables */
    parasail_profile_t *profile = NULL;
    parasail_result_t *result = NULL;

    /* validate inputs */
    PARASAIL_CHECK_NULL(s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);
    PARASAIL_CHECK_NULL(matrix);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(s1);
        PARASAIL_CHECK_GT0(s1Len);
    }

    /* initialize local variables */
    profile = parasail_profile_create_sse_128_16(s1, s1Len, matrix);
    if (!profile) return NULL;
    result = PNAME(profile, s2, s2Len, open, gap, s1_beg, s1_end, s2_beg, s2_end);

    parasail_profile_free(profile);

    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
{
    /* declare local variables */
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t s1Len = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const parasail_matrix_t *matrix = NULL;
    int32_t segWidth = 0;
    int32_t segLen = 0;
    int32_t offset = 0;
    int32_t position = 0;
    __m128i* restrict vProfile = NULL;
    __m128i* restrict pvHStore = NULL;
    __m128i* restrict pvHLoad = NULL;
    __m128i* restrict pvE = NULL;
    __m128i* restrict pvEaStore = NULL;
    __m128i* restrict pvEaLoad = NULL;
    __m128i* restrict pvHT = NULL;
    int16_t* restrict boundary = NULL;
    __m128i vGapO;
    __m128i vGapE;
    int16_t NEG_LIMIT = 0;
    int16_t POS_LIMIT = 0;
    int16_t score = 0;
    __m128i vNegLimit;
    __m128i vPosLimit;
    __m128i vSaturationCheckMin;
    __m128i vSaturationCheckMax;
    __m128i vMaxH;
    __m128i vPosMask;
    parasail_result_t *result = NULL;
    __m128i vTIns;
    __m128i vTDel;
    __m128i vTDiag;
    __m128i vTDiagE;
    __m128i vTInsE;
    __m128i vTDiagF;
    __m128i vTDelF;
    __m128i vTMask;
    __m128i vFTMask;

    /* validate inputs */
    PARASAIL_CHECK_NULL(profile);
    PARASAIL_CHECK_NULL(profile->profile16.score);
    PARASAIL_CHECK_NULL(profile->matrix);
    PARASAIL_CHECK_GT0(profile->s1Len);
    PARASAIL_CHECK_NULL(s2);
    PARASAIL_CHECK_GT0(s2Len);
    PARASAIL_CHECK_GE0(open);
    PARASAIL_CHECK_GE0(gap);

    /* initialize stack variables */
    i = 0;
    j = 0;
    k = 0;
    s1Len = profile->s1Len;
    end_query = s1Len-1;
    end_ref = s2Len-1;
    matrix = profile->matrix;
    segWidth = 8; /* number of values in vector unit */
    segLen = (s1Len + segWidth - 1) / segWidth;
    offset = (s1Len - 1) % segLen;
    position = (segWidth - 1) - (s1Len - 1) / segLen;
    vProfile = (__m128i*)profile->profile16.score;
    vGapO = _mm_set1_epi16(open);
    vGapE = _mm_set1_epi16(gap);
    NEG_LIMIT = (-open < matrix->min ? INT16_MIN + open : INT16_MIN - matrix->min) + 1;
    POS_LIMIT = INT16_MAX - matrix->max - 1;
    score = NEG_LIMIT;
    vNegLimit = _mm_set1_epi16(NEG_LIMIT);
    vPosLimit = _mm_set1_epi16(POS_LIMIT);
    vSaturationCheckMin = vPosLimit;
    vSaturationCheckMax = vNegLimit;
    vMaxH = vNegLimit;
    vPosMask = _mm_cmpeq_epi16(_mm_set1_epi16(position),
            _mm_set_epi16(0,1,2,3,4,5,6,7));
    vTIns  = _mm_set1_epi16(PARASAIL_INS);
    vTDel  = _mm_set1_epi16(PARASAIL_DEL);
    vTDiag = _mm_set1_epi16(PARASAIL_DIAG);
    vTDiagE= _mm_set1_epi16(PARASAIL_DIAG_E);
    vTInsE = _mm_set1_epi16(PARASAIL_INS_E);
    vTDiagF= _mm_set1_epi16(PARASAIL_DIAG_F);
    vTDelF = _mm_set1_epi16(PARASAIL_DEL_F);
    vTMask = _mm_set1_epi16(PARASAIL_ZERO_MASK);
    vFTMask= _mm_set1_epi16(PARASAIL_F_MASK);

    /* initialize result */
    result = parasail_result_new_trace(segLen, s2Len, 16, sizeof(__m128i));
    if (!result) return NULL;

    /* set known flags */
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_STRIPED
        | PARASAIL_FLAG_TRACE
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_8;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;

    /* initialize heap variables */
    pvHStore = parasail_memalign___m128i(16, segLen);
    pvHLoad  = parasail_memalign___m128i(16, segLen);
    pvE      = parasail_memalign___m128i(16, segLen);
    pvEaStore= parasail_memalign___m128i(16, segLen);
    pvEaLoad = parasail_memalign___m128i(16, segLen);
    pvHT     = parasail_memalign___m128i(16, segLen);
    boundary = parasail_memalign_int16_t(16, s2Len+1);

    /* validate heap variables */
    if (!pvHStore) return NULL;
    if (!pvHLoad) return NULL;
    if (!pvE) return NULL;
    if (!pvEaStore) return NULL;
    if (!pvEaLoad) return NULL;
    if (!pvHT) return NULL;
    if (!boundary) return NULL;

    /* initialize H and E */
    {
        int32_t index = 0;
        for (i=0; i<segLen; ++i) {
            int32_t segNum = 0;
            __m128i_16_t h;
            __m128i_16_t e;
            for (segNum=0; segNum<segWidth; ++segNum) {
                int64_t tmp = s1_beg ? 0 : (-open-gap*(segNum*segLen+i));
                h.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
                tmp = tmp - open;
                e.v[segNum] = tmp < INT16_MIN ? INT16_MIN : tmp;
            }
            _mm_store_si128(&pvHStore[index], h.m);
            _mm_store_si128(&pvE[index], e.m);
            _mm_store_si128(&pvEaStore[index], e.m);
            ++index;
        }
    }

    /* initialize uppder boundary */
    {
        boundary[0] = 0;
        for (i=1; i<=s2Len; ++i) {
            int64_t tmp = s2_beg ? 0 : (-open-gap*(i-1));
            boundary[i] = tmp < INT16_MIN ? INT16_MIN : tmp;
        }
    }

    for (i=0; i<segLen; ++i) {
        arr_store(result->trace->trace_table, vTDiagE, i, segLen, 0);
    }

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vEF_opn;
        __m128i vE;
        __m128i vE_ext;
        __m128i vF;
        __m128i vF_ext;
        __m128i vFa;
        __m128i vFa_ext;
        __m128i vH;
        __m128i vH_dag;
        const __m128i* vP = NULL;

        /* Initialize F value to -inf.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = vNegLimit;

        /* load final segment of pvHStore and shift left by 2 bytes */
        vH = _mm_load_si128(&pvHStore[segLen - 1]);
        vH = _mm_slli_si128(vH, 2);

        /* insert upper boundary condition */
        vH = _mm_insert_epi16(vH, boundary[j], 0);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        /* Swap the 2 H buffers. */
        SWAP(pvHLoad, pvHStore)
        SWAP(pvEaLoad, pvEaStore)

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm_adds_epi16(vH, _mm_load_si128(vP + i));
            vH = _mm_max_epi16(vH_dag, vE);
            vH = _mm_max_epi16(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);
            vSaturationCheckMax = _mm_max_epi16(vSaturationCheckMax, vH);
            vSaturationCheckMin = _mm_min_epi16(vSaturationCheckMin, vH);
            vSaturationCheckMin = _mm_min_epi16(vSaturationCheckMin, vE);
            vSaturationCheckMin = _mm_min_epi16(vSaturationCheckMin, vF);

            {
                __m128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                __m128i case1 = _mm_cmpeq_epi16(vH, vH_dag);
                __m128i case2 = _mm_cmpeq_epi16(vH, vF);
                __m128i vT = _mm_blendv_epi8(
                        _mm_blendv_epi8(vTIns, vTDel, case2),
                        vTDiag, case1);
                _mm_store_si128(pvHT + i, vT);
                vT = _mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i, segLen, j);
            }

            vEF_opn = _mm_subs_epi16(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm_subs_epi16(vE, vGapE);
            vE = _mm_max_epi16(vEF_opn, vE_ext);
            _mm_store_si128(pvE + i, vE);
            {
                __m128i vEa = _mm_load_si128(pvEaLoad + i);
                __m128i vEa_ext = _mm_subs_epi16(vEa, vGapE);
                vEa = _mm_max_epi16(vEF_opn, vEa_ext);
                _mm_store_si128(pvEaStore + i, vEa);
                if (j+1<s2Len) {
                    __m128i cond = _mm_cmpgt_epi16(vEF_opn, vEa_ext);
                    __m128i vT = _mm_blendv_epi8(vTInsE, vTDiagE, cond);
                    arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                }
            }

            /* Update vF value. */
            vF_ext = _mm_subs_epi16(vF, vGapE);
            vF = _mm_max_epi16(vEF_opn, vF_ext);
            if (i+1<segLen) {
                __m128i vTAll = arr_load(result->trace->trace_table, i+1, segLen, j);
                __m128i cond = _mm_cmpgt_epi16(vEF_opn, vF_ext);
                __m128i vT = _mm_blendv_epi8(vTDelF, vTDiagF, cond);
                vT = _mm_or_si128(vT, vTAll);
                arr_store(result->trace->trace_table, vT, i+1, segLen, j);
            }

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        vFa_ext = vF_ext;
        vFa = vF;
        for (k=0; k<segWidth; ++k) {
            int64_t tmp = s2_beg ? -open : (boundary[j+1]-open);
            int16_t tmp2 = tmp < INT16_MIN ? INT16_MIN : tmp;
            __m128i vHp = _mm_load_si128(&pvHLoad[segLen - 1]);
            vHp = _mm_slli_si128(vHp, 2);
            vHp = _mm_insert_epi16(vHp, boundary[j], 0);
            vEF_opn = _mm_slli_si128(vEF_opn, 2);
            vEF_opn = _mm_insert_epi16(vEF_opn, tmp2, 0);
            vF_ext = _mm_slli_si128(vF_ext, 2);
            vF_ext = _mm_insert_epi16(vF_ext, NEG_LIMIT, 0);
            vF = _mm_slli_si128(vF, 2);
            vF = _mm_insert_epi16(vF, tmp2, 0);
            vFa_ext = _mm_slli_si128(vFa_ext, 2);
            vFa_ext = _mm_insert_epi16(vFa_ext, NEG_LIMIT, 0);
            vFa = _mm_slli_si128(vFa, 2);
            vFa = _mm_insert_epi16(vFa, tmp2, 0);
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi16(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
                vSaturationCheckMin = _mm_min_epi16(vSaturationCheckMin, vH);
                vSaturationCheckMax = _mm_max_epi16(vSaturationCheckMax, vH);
                {
                    __m128i vTAll;
                    __m128i vT;
                    __m128i case1;
                    __m128i case2;
                    __m128i cond;
                    vHp = _mm_adds_epi16(vHp, _mm_load_si128(vP + i));
                    case1 = _mm_cmpeq_epi16(vH, vHp);
                    case2 = _mm_cmpeq_epi16(vH, vF);
                    cond = _mm_andnot_si128(case1,case2);
                    vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    vT = _mm_load_si128(pvHT + i);
                    vT = _mm_blendv_epi8(vT, vTDel, cond);
                    _mm_store_si128(pvHT + i, vT);
                    vTAll = _mm_and_si128(vTAll, vTMask);
                    vTAll = _mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                /* Update vF value. */
                {
                    __m128i vTAll = arr_load(result->trace->trace_table, i, segLen, j);
                    __m128i cond = _mm_cmpgt_epi16(vEF_opn, vFa_ext);
                    __m128i vT = _mm_blendv_epi8(vTDelF, vTDiagF, cond);
                    vTAll = _mm_and_si128(vTAll, vFTMask);
                    vTAll = _mm_or_si128(vTAll, vT);
                    arr_store(result->trace->trace_table, vTAll, i, segLen, j);
                }
                vEF_opn = _mm_subs_epi16(vH, vGapO);
                vF_ext = _mm_subs_epi16(vF, vGapE);
                {
                    __m128i vEa = _mm_load_si128(pvEaLoad + i);
                    __m128i vEa_ext = _mm_subs_epi16(vEa, vGapE);
                    vEa = _mm_max_epi16(vEF_opn, vEa_ext);
                    _mm_store_si128(pvEaStore + i, vEa);
                    if (j+1<s2Len) {
                        __m128i cond = _mm_cmpgt_epi16(vEF_opn, vEa_ext);
                        __m128i vT = _mm_blendv_epi8(vTInsE, vTDiagE, cond);
                        arr_store(result->trace->trace_table, vT, i, segLen, j+1);
                    }
                }
                if (! _mm_movemask_epi8(
                            _mm_or_si128(
                                _mm_cmpgt_epi16(vF_ext, vEF_opn),
                                _mm_cmpeq_epi16(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm_max_epi16(vEF_opn, vF_ext);*/
                vF = vF_ext;
                vFa_ext = _mm_subs_epi16(vFa, vGapE);
                vFa = _mm_max_epi16(vEF_opn, vFa_ext);
                vHp = _mm_load_si128(pvHLoad + i);
            }
        }
end:
        {
            /* extract vector containing last value from the column */
            __m128i vCompare;
            vH = _mm_load_si128(pvHStore + offset);
            vCompare = _mm_and_si128(vPosMask, _mm_cmpgt_epi16(vH, vMaxH));
            vMaxH = _mm_max_epi16(vH, vMaxH);
            if (_mm_movemask_epi8(vCompare)) {
                end_ref = j;
            }
        }
    }

    /* max last value from all columns */
    if (s2_end)
    {
        for (k=0; k<position; ++k) {
            vMaxH = _mm_slli_si128(vMaxH, 2);
        }
        score = (int16_t) _mm_extract_epi16(vMaxH, 7);
        end_query = s1Len-1;
    }

    /* max of last column */
    if (s1_end)
    {
        /* Trace the alignment ending position on read. */
        int16_t *t = (int16_t*)pvHStore;
        int32_t column_len = segLen * segWidth;
        for (i = 0; i<column_len; ++i, ++t) {
            int32_t temp = i / segWidth + i % segWidth * segLen;
            if (temp >= s1Len) continue;
            if (*t > score) {
                score = *t;
                end_query = temp;
                end_ref = s2Len-1;
            }
            else if (*t == score && end_ref == s2Len-1 && temp < end_query) {
                end_query = temp;
            }
        }
    }

    if (!s1_end && !s2_end) {
        /* extract last value from the last column */
        {
            __m128i vH = _mm_load_si128(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 2);
            }
            score = (int16_t) _mm_extract_epi16 (vH, 7);
            end_ref = s2Len - 1;
            end_query = s1Len - 1;
        }
    }

    if (_mm_movemask_epi8(_mm_or_si128(
            _mm_cmplt_epi16(vSaturationCheckMin, vNegLimit),
            _mm_cmpgt_epi16(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(boundary);
    parasail_free(pvHT);
    parasail_free(pvEaLoad);
    parasail_free(pvEaStore);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}

SG_IMPL_ALL
SG_IMPL_PROF_ALL

#endif  //  HAVE_INTEL
