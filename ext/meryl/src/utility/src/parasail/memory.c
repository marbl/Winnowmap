/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "memory.h"

/* validate inputs */
#define VALIDATE_A_B                                       \
if (a <= 0 || b <= 0) {                                    \
    fprintf(stderr, "%s: inputs must be > 0\n", __func__); \
    return NULL;                                           \
}

void* parasail_memalign(size_t alignment, size_t size)
{
    void *ptr = NULL;
#if defined(HAVE__ALIGNED_MALLOC)
    ptr = _aligned_malloc(size, alignment);
#elif defined(HAVE_POSIX_MEMALIGN)
    int retcode = posix_memalign(&ptr, alignment, size);
    if (0 != retcode) {
        fprintf(stderr, "%s: posix_memalign failed: %s\n",
                __func__, strerror(retcode));
        return NULL;
    }
#elif defined(HAVE_ALIGNED_ALLOC)
    ptr = aligned_alloc(alignment, size);
#elif defined(HAVE_MEMALIGN)
    ptr = memalign(alignment, size);
#else
#error "No suitable memory alignment routine found."
#endif
    if (NULL == ptr) {
        fprintf(stderr, "%s: failed\n", __func__);
    }
    return ptr;
}

void parasail_free(void *ptr)
{
#if defined(HAVE__ALIGNED_MALLOC)
     _aligned_free(ptr);
#else
    free(ptr);
#endif
}

void parasail_free_unaligned(void *ptr)
{
    free(ptr);
}

int * parasail_memalign_int(size_t alignment, size_t size)
{
    return (int *) parasail_memalign(alignment, size*sizeof(int));
}

int8_t * parasail_memalign_int8_t(size_t alignment, size_t size)
{
    return (int8_t *) parasail_memalign(alignment, size*sizeof(int8_t));
}

int16_t * parasail_memalign_int16_t(size_t alignment, size_t size)
{
    return (int16_t *) parasail_memalign(alignment, size*sizeof(int16_t));
}

int32_t * parasail_memalign_int32_t(size_t alignment, size_t size)
{
    return (int32_t *) parasail_memalign(alignment, size*sizeof(int32_t));
}

int64_t * parasail_memalign_int64_t(size_t alignment, size_t size)
{
    return (int64_t *) parasail_memalign(alignment, size*sizeof(int64_t));
}

void parasail_memset(void *b, int c, size_t len)
{
    (void)memset(b, c, len);
}

void parasail_memset_int(int *b, int c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int8_t(int8_t *b, int8_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int16_t(int16_t *b, int16_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int32_t(int32_t *b, int32_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

void parasail_memset_int64_t(int64_t *b, int64_t c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        b[i] = c;
    }
}

parasail_result_t* parasail_result_new(void)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    PARASAIL_NEW(result, parasail_result_t);

    result->score = 0;
    result->end_query = 0;
    result->end_ref = 0;
    result->flag = 0;
    result->extra = NULL;

    return result;
}

parasail_result_t* parasail_result_new_stats(void)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* allocate struct to hold memory */
    result = parasail_result_new();
    if (!result) return NULL;

    /* allocate only tables */
    PARASAIL_NEW(result->stats, parasail_result_extra_stats_t);

    return result;
}

parasail_result_t* parasail_result_new_table1(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    /* validate inputs */
    VALIDATE_A_B
    
    /* allocate struct to hold memory */
    result = parasail_result_new();
    if (!result) return NULL;

    /* allocate only score table */
    PARASAIL_NEW(result->tables, parasail_result_extra_tables_t);
    PARASAIL_CALLOC(result->tables->score_table, int, a*b);

    return result;
}

parasail_result_t* parasail_result_new_rowcol1(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    VALIDATE_A_B
    
    /* allocate struct to hold memory */
    result = parasail_result_new();
    if (!result) return NULL;

    /* allocate only score col and row */
    PARASAIL_NEW(result->rowcols, parasail_result_extra_rowcols_t);
    PARASAIL_CALLOC(result->rowcols->score_row, int, b);
    PARASAIL_CALLOC(result->rowcols->score_col, int, a);

    return result;
}

parasail_result_t* parasail_result_new_table3(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    VALIDATE_A_B
    
    /* allocate struct to hold memory */
    result = parasail_result_new();
    if (!result) return NULL;
    
    /* allocate only tables */
    PARASAIL_NEW(result->stats, parasail_result_extra_stats_t);
    PARASAIL_NEW(result->stats->tables, parasail_result_extra_stats_tables_t);
    PARASAIL_CALLOC(result->stats->tables->score_table, int, a*b);
    PARASAIL_CALLOC(result->stats->tables->matches_table, int, a*b);
    PARASAIL_CALLOC(result->stats->tables->similar_table, int, a*b);
    PARASAIL_CALLOC(result->stats->tables->length_table, int, a*b);

    return result;
}

parasail_result_t* parasail_result_new_trace(const int a, const int b, const size_t alignment, const size_t size)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    VALIDATE_A_B

    /* allocate struct to hold memory */
    result = parasail_result_new();
    if (!result) return NULL;

    PARASAIL_NEW(result->trace, parasail_result_extra_trace_t);
    result->trace->trace_table = parasail_memalign(alignment, size*a*b);
    if (!result->trace->trace_table) return NULL;
    result->trace->trace_ins_table = NULL;
    result->trace->trace_del_table = NULL;

    return result;
}

parasail_result_t* parasail_result_new_rowcol3(const int a, const int b)
{
    /* declare all variables */
    parasail_result_t *result = NULL;

    VALIDATE_A_B
    
    /* allocate struct to hold memory */
    result = parasail_result_new();
    if (!result) return NULL;
    
    PARASAIL_NEW(result->stats, parasail_result_extra_stats_t);
    PARASAIL_NEW(result->stats->rowcols, parasail_result_extra_stats_rowcols_t);
    PARASAIL_CALLOC(result->stats->rowcols->score_row, int, b);
    PARASAIL_CALLOC(result->stats->rowcols->matches_row, int, b);
    PARASAIL_CALLOC(result->stats->rowcols->similar_row, int, b);
    PARASAIL_CALLOC(result->stats->rowcols->length_row, int, b);

    PARASAIL_CALLOC(result->stats->rowcols->score_col, int, a);
    PARASAIL_CALLOC(result->stats->rowcols->matches_col, int, a);
    PARASAIL_CALLOC(result->stats->rowcols->similar_col, int, a);
    PARASAIL_CALLOC(result->stats->rowcols->length_col, int, a);

    return result;
}

void parasail_result_free(parasail_result_t *result)
{
    /* validate inputs */
    if (!result) {
        fprintf(stderr, "%s: attempted free of NULL result pointer\n", __func__);
        return;
    }
    
    if (result->flag & PARASAIL_FLAG_STATS) {
        if (result->flag & PARASAIL_FLAG_TABLE) {
            free(result->stats->tables->score_table);
            free(result->stats->tables->matches_table);
            free(result->stats->tables->similar_table);
            free(result->stats->tables->length_table);
            free(result->stats->tables);
        }
        if (result->flag & PARASAIL_FLAG_ROWCOL) {
            free(result->stats->rowcols->score_row);
            free(result->stats->rowcols->matches_row);
            free(result->stats->rowcols->similar_row);
            free(result->stats->rowcols->length_row);
            free(result->stats->rowcols->score_col);
            free(result->stats->rowcols->matches_col);
            free(result->stats->rowcols->similar_col);
            free(result->stats->rowcols->length_col);
            free(result->stats->rowcols);
        }
        free(result->stats);
    }
    else {
        if (result->flag & PARASAIL_FLAG_TABLE) {
            free(result->tables->score_table);
            free(result->tables);
        }
        if (result->flag & PARASAIL_FLAG_ROWCOL) {
            free(result->rowcols->score_row);
            free(result->rowcols->score_col);
            free(result->rowcols);
        }
    }
    if (result->flag & PARASAIL_FLAG_TRACE) {
        parasail_free(result->trace->trace_table);
        if (NULL != result->trace->trace_ins_table)
            parasail_free(result->trace->trace_ins_table);
        if (NULL != result->trace->trace_del_table)
            parasail_free(result->trace->trace_del_table);
        free(result->trace);
    }

    free(result);
}

void parasail_version(int *major, int *minor, int *patch)
{
    *major = PARASAIL_VERSION_MAJOR;
    *minor = PARASAIL_VERSION_MINOR;
    *patch = PARASAIL_VERSION_PATCH;
}

static parasail_matrix_t* parasail_matrix_create_internal(
        const char *alphabet, const int match, const int mismatch, int case_sensitive)
{
    parasail_matrix_t *retval = NULL;
    int *matrix = NULL;
    int *mapper = NULL;
    char *alphabet_copy = NULL;
    size_t size = 0;
    size_t size1 = 0;
    size_t i = 0;
    size_t j = 0;
    size_t c = 0;

    PARASAIL_CHECK_NULL(alphabet);

    size = strlen(alphabet);
    if (size >= INT_MAX) {
        fprintf(stderr, "%s: alphabet is too large\n", __func__);
        return NULL;
    }
    size1 = size + 1;

    PARASAIL_CALLOC(matrix, int, size1*size1);

    for (i=0; i<size; ++i) {
        for (j=0; j<size; ++j) {
            if (i == j) {
                matrix[c++] = match;
            }
            else {
                matrix[c++] = mismatch;
            }
        }
        matrix[c++] = 0;
    }
    for (j=0; j<size1; ++j) {
        matrix[c++] = 0;
    }

    PARASAIL_CALLOC(mapper, int, 256);
    parasail_memset_int(mapper, (int)size, 256);
    if (case_sensitive) {
        for (i=0; i<size; ++i) {
            mapper[(unsigned char)alphabet[i]] = (int)i;
        }
    }
    else {
        for (i=0; i<size; ++i) {
            mapper[toupper((unsigned char)alphabet[i])] = (int)i;
            mapper[tolower((unsigned char)alphabet[i])] = (int)i;
        }
    }

    PARASAIL_CALLOC(alphabet_copy, char, size1+1);
    (void)memcpy(alphabet_copy, alphabet, sizeof(char) * size);
    alphabet_copy[size] = '*';
    alphabet_copy[size1] = '\0';

    PARASAIL_NEW(retval, parasail_matrix_t);
    retval->name = "";
    retval->matrix = matrix;
    retval->mapper = mapper;
    retval->size = (int)size1;
    retval->max = match > mismatch ? match : mismatch;
    retval->min = match > mismatch ? mismatch : match;
    retval->user_matrix = matrix;
    retval->type = PARASAIL_MATRIX_TYPE_SQUARE;
    retval->length = retval->size;
    retval->alphabet = alphabet_copy;
    retval->query = NULL;
    return retval;
}

parasail_matrix_t* parasail_matrix_create(
        const char *alphabet, const int match, const int mismatch)
{
    return parasail_matrix_create_internal(alphabet, match, mismatch, 0);
}

parasail_matrix_t* parasail_matrix_create_case_sensitive(
        const char *alphabet, const int match, const int mismatch)
{
    return parasail_matrix_create_internal(alphabet, match, mismatch, 1);
}

static parasail_matrix_t* parasail_matrix_pssm_create_internal(
        const char *alphabet, const int *values, const int length, int case_sensitive)
{
    parasail_matrix_t *retval = NULL;
    int *matrix = NULL;
    int *mapper = NULL;
    char *alphabet_copy = NULL;
    size_t size = 0;
    size_t size1 = 0;
    size_t i = 0;
    size_t j = 0;
    size_t c = 0;
    size_t v = 0;
    int min = INT_MAX;
    int max = INT_MIN;

    PARASAIL_CHECK_NULL(alphabet);
    PARASAIL_CHECK_NULL(values);
    if (length <= 0) {
        fprintf(stderr, "%s: length must be > 0\n", __func__);
        return NULL;
    }

    size = strlen(alphabet);
    if (size >= INT_MAX) {
        fprintf(stderr, "%s: alphabet is too large\n", __func__);
        return NULL;
    }
    size1 = size + 1;

    PARASAIL_CALLOC(matrix, int, size1*length);

    /* find min and max values */
    for (i=0; i<size*length; ++i) {
        min = values[i] < min ? values[i] : min;
        max = values[i] > max ? values[i] : max;
    }

    for (i=0; i<(size_t)length; ++i) {
        for (j=0; j<size; ++j) {
            matrix[c++] = values[v++];
        }
        matrix[c++] = min;
    }

    PARASAIL_CALLOC(mapper, int, 256);
    parasail_memset_int(mapper, (int)size, 256);
    if (case_sensitive) {
        for (i=0; i<size; ++i) {
            mapper[(unsigned char)alphabet[i]] = (int)i;
        }
    }
    else {
        for (i=0; i<size; ++i) {
            mapper[toupper((unsigned char)alphabet[i])] = (int)i;
            mapper[tolower((unsigned char)alphabet[i])] = (int)i;
        }
    }

    PARASAIL_CALLOC(alphabet_copy, char, size1+1);
    (void)memcpy(alphabet_copy, alphabet, sizeof(char) * size);
    alphabet_copy[size] = '*';
    alphabet_copy[size1] = '\0';

    PARASAIL_NEW(retval, parasail_matrix_t);
    retval->name = "";
    retval->matrix = matrix;
    retval->mapper = mapper;
    retval->size = (int)size1;
    retval->max = max;
    retval->min = min;
    retval->user_matrix = matrix;
    retval->type = PARASAIL_MATRIX_TYPE_PSSM;
    retval->length = length;
    retval->alphabet = alphabet_copy;
    retval->query = NULL;
    return retval;
}

parasail_matrix_t* parasail_matrix_pssm_create(
        const char *alphabet, const int *values, const int length)
{
    return parasail_matrix_pssm_create_internal(alphabet, values, length, 0);
}

parasail_matrix_t* parasail_matrix_pssm_create_case_sensitive(
        const char *alphabet, const int *values, const int length)
{
    return parasail_matrix_pssm_create_internal(alphabet, values, length, 1);
}

parasail_matrix_t* parasail_matrix_copy(const parasail_matrix_t *matrix)
{
    parasail_matrix_t *retval = NULL;

    PARASAIL_CHECK_NULL(matrix);

    PARASAIL_NEW(retval, parasail_matrix_t);
    retval->name = matrix->name;
    retval->size = matrix->size;
    retval->max = matrix->max;
    retval->min = matrix->min;
    retval->type = matrix->type;
    retval->length = matrix->length;

    {
        size_t matrix_size = 0;
        size_t alphabet_size = (matrix->size+1);
        int *new_mapper = NULL;
        int *new_matrix = NULL;
        char *new_alphabet = NULL;
        char *new_query = NULL;

        if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
            matrix_size = matrix->size*matrix->size;
        }
        else if (matrix->type == PARASAIL_MATRIX_TYPE_PSSM) {
            matrix_size = matrix->size*matrix->length;
        }

        PARASAIL_CALLOC(new_mapper, int, 256);
        (void)memcpy(new_mapper, matrix->mapper, sizeof(int)*256);

        PARASAIL_CALLOC(new_matrix, int, matrix_size);
        (void)memcpy(new_matrix, matrix->matrix, sizeof(int)*matrix_size);

        PARASAIL_CALLOC(new_alphabet, char, alphabet_size);
        (void)memcpy(new_alphabet, matrix->alphabet, sizeof(char)*alphabet_size);

        if (matrix->query) {
            size_t query_size = strlen(matrix->query);
            PARASAIL_CALLOC(new_query, char, query_size+1);
            (void)memcpy(new_query, matrix->query, sizeof(char)*(query_size+1));
            new_query[query_size] = '\0';
        }

        retval->mapper = new_mapper;
        retval->matrix = new_matrix;
        retval->user_matrix = new_matrix;
        retval->alphabet = new_alphabet;
        retval->query = new_query;
    }

    return retval;
}

parasail_matrix_t* parasail_matrix_convert_square_to_pssm(
        const parasail_matrix_t *matrix,
        const char *s1,
        int s1Len)
{
    parasail_matrix_t *retval = NULL;
    size_t matrix_size = matrix->size*s1Len;
    size_t alphabet_size = matrix->size+1;
    int *new_mapper = NULL;
    int *new_matrix = NULL;
    char *new_alphabet = NULL;
    char *new_query = NULL;
    int i;

    PARASAIL_CHECK_NULL(matrix);

    if (matrix->type != PARASAIL_MATRIX_TYPE_SQUARE) {
        fprintf(stderr, "%s: attempted to convert non-square matrix to pssm\n", __func__);
        return NULL;
    }

    PARASAIL_NEW(retval, parasail_matrix_t);

    PARASAIL_CALLOC(new_mapper, int, 256);
    (void)memcpy(new_mapper, matrix->mapper, sizeof(int)*256);

    PARASAIL_CALLOC(new_matrix, int, matrix_size);

    PARASAIL_CALLOC(new_alphabet, char, alphabet_size);
    (void)memcpy(new_alphabet, matrix->alphabet, sizeof(char)*alphabet_size);

    PARASAIL_CALLOC(new_query, char, s1Len+1);
    (void)memcpy(new_query, s1, sizeof(char)*(s1Len+1));

    for (i=0; i<s1Len; ++i) {
        (void)memcpy(
                &new_matrix[matrix->size*i],
                &matrix->matrix[matrix->size*matrix->mapper[(unsigned char)s1[i]]],
                sizeof(int)*matrix->size);
    }

    retval->name = matrix->name;
    retval->matrix = new_matrix;
    retval->mapper = new_mapper;
    retval->size = matrix->size;
    retval->max = matrix->max;
    retval->min = matrix->min;
    retval->user_matrix = new_matrix;
    retval->type = PARASAIL_MATRIX_TYPE_PSSM;
    retval->length = s1Len;
    retval->alphabet = new_alphabet;
    retval->query = new_query;

    return retval;
}

void parasail_matrix_set_value(parasail_matrix_t *matrix, int row, int col, int value)
{
    PARASAIL_CHECK_NULL_NORETVAL(matrix);

    if (NULL == matrix->user_matrix) {
        fprintf(stderr, "%s: attempted to set value of built-in matrix '%s'\n",
                __func__, matrix->name);
        return;
    }

    matrix->user_matrix[row*matrix->size + col] = value;
    if (value > matrix->max) {
        matrix->max = value;
    }
    if (value < matrix->min) {
        matrix->min = value;
    }
}

void parasail_matrix_free(parasail_matrix_t *matrix)
{
    /* validate inputs */
    PARASAIL_CHECK_NULL_NORETVAL(matrix);

    if (NULL != matrix->user_matrix) {
        free((void*)matrix->matrix);
        free((void*)matrix->mapper);
        free((void*)matrix->alphabet);
        if (matrix->query) free((void*)matrix->query);
        free(matrix);
    }
    else {
        fprintf(stderr, "%s: attempted to free built-in matrix '%s'\n",
                __func__, matrix->name);
    }
}

parasail_profile_t* parasail_profile_new(
        const char * s1, const int s1Len, const parasail_matrix_t *matrix)
{
    /* declare all variables */
    parasail_profile_t *profile = NULL;

    PARASAIL_CHECK_NULL(matrix);
    if (matrix->type == PARASAIL_MATRIX_TYPE_SQUARE) {
        PARASAIL_CHECK_NULL(s1);
    }

    PARASAIL_NEW(profile, parasail_profile_t);

    profile->s1 = s1;
    profile->s1Len = s1Len;
    profile->matrix = matrix;
    profile->profile8.score = NULL;
    profile->profile8.matches = NULL;
    profile->profile8.similar = NULL;
    profile->profile16.score = NULL;
    profile->profile16.matches = NULL;
    profile->profile16.similar = NULL;
    profile->profile32.score = NULL;
    profile->profile32.matches = NULL;
    profile->profile32.similar = NULL;
    profile->profile64.score = NULL;
    profile->profile64.matches = NULL;
    profile->profile64.similar = NULL;
    profile->free = NULL;
    profile->stop = INT32_MAX;

    if (matrix->type == PARASAIL_MATRIX_TYPE_PSSM) {
        profile->s1Len = matrix->length;
    }

    return profile;
}

void parasail_profile_free(parasail_profile_t *profile)
{
    if (!profile) {
        fprintf(stderr, "%s: attempted free of NULL profile pointer\n", __func__);
        return;
    }

    if (NULL != profile->profile8.score) {
        profile->free(profile->profile8.score);
    }
    if (NULL != profile->profile8.matches) {
        profile->free(profile->profile8.matches);
    }
    if (NULL != profile->profile8.similar) {
        profile->free(profile->profile8.similar);
    }

    if (NULL != profile->profile16.score) {
        profile->free(profile->profile16.score);
    }
    if (NULL != profile->profile16.matches) {
        profile->free(profile->profile16.matches);
    }
    if (NULL != profile->profile16.similar) {
        profile->free(profile->profile16.similar);
    }

    if (NULL != profile->profile32.score) {
        profile->free(profile->profile32.score);
    }
    if (NULL != profile->profile32.matches) {
        profile->free(profile->profile32.matches);
    }
    if (NULL != profile->profile32.similar) {
        profile->free(profile->profile32.similar);
    }

    if (NULL != profile->profile64.score) {
        profile->free(profile->profile64.score);
    }
    if (NULL != profile->profile64.matches) {
        profile->free(profile->profile64.matches);
    }
    if (NULL != profile->profile64.similar) {
        profile->free(profile->profile64.similar);
    }

    free(profile);
}

char* parasail_reverse(const char *s, size_t length)
{
    char *r = NULL;
    size_t i = 0;
    size_t j = 0;

    PARASAIL_CALLOC(r, char, length+1);
    r[length] = '\0';
    for (i=0,j=length-1; i<length; ++i,--j) {
        r[i] = s[j];
    }

    return r;
}

uint32_t* parasail_reverse_uint32_t(const uint32_t *s, size_t length)
{
    uint32_t *r = NULL;
    size_t i = 0;
    size_t j = 0;

    PARASAIL_CALLOC(r, uint32_t, length);
    for (i=0,j=length-1; i<length; ++i,--j) {
        r[i] = s[j];
    }

    return r;
}

int parasail_result_is_nw(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_NW;
}

int parasail_result_is_sg(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_SG;
}

int parasail_result_is_sw(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_SW;
}

int parasail_result_is_saturated(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_SATURATED;
}

int parasail_result_is_banded(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_BANDED;
}

int parasail_result_is_scan(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_SCAN;
}

int parasail_result_is_striped(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_STRIPED;
}

int parasail_result_is_diag(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_DIAG;
}

int parasail_result_is_blocked(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_BLOCKED;
}

int parasail_result_is_stats(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_STATS;
}

int parasail_result_is_stats_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return (result->flag & PARASAIL_FLAG_STATS) && (result->flag & PARASAIL_FLAG_TABLE);
}

int parasail_result_is_stats_rowcol(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return (result->flag & PARASAIL_FLAG_STATS) && (result->flag & PARASAIL_FLAG_ROWCOL);
}

int parasail_result_is_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_TABLE;
}

int parasail_result_is_rowcol(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_ROWCOL;
}

int parasail_result_is_trace(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->flag & PARASAIL_FLAG_TRACE;
}

int parasail_result_get_score(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->score;
}

int parasail_result_get_end_query(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->end_query;
}

int parasail_result_get_end_ref(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    return result->end_ref;
}

int parasail_result_get_matches(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats(result));
    return result->stats->matches;
}

int parasail_result_get_similar(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats(result));
    return result->stats->similar;
}

int parasail_result_get_length(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats(result));
    return result->stats->length;
}

int* parasail_result_get_score_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_table(result) || parasail_result_is_stats_table(result));
    if (parasail_result_is_stats_table(result)) {
        return result->stats->tables->score_table;
    }
    if (parasail_result_is_table(result)) {
        return result->tables->score_table;
    }
    return NULL; /* should not reach */
}

int* parasail_result_get_matches_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_table(result));
    return result->stats->tables->matches_table;
}

int* parasail_result_get_similar_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_table(result));
    return result->stats->tables->similar_table;
}

int* parasail_result_get_length_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_table(result));
    return result->stats->tables->length_table;
}

int* parasail_result_get_score_row(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result) || parasail_result_is_rowcol(result));
    if (parasail_result_is_stats_rowcol(result)) {
        return result->stats->rowcols->score_row;
    }
    if (parasail_result_is_rowcol(result)) {
        return result->rowcols->score_row;
    }
    return NULL; /* should not reach */
}

int* parasail_result_get_matches_row(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->matches_row;
}

int* parasail_result_get_similar_row(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->similar_row;
}

int* parasail_result_get_length_row(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->length_row;
}

int* parasail_result_get_score_col(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result) || parasail_result_is_rowcol(result));
    if (parasail_result_is_stats_rowcol(result)) {
        return result->stats->rowcols->score_col;
    }
    if (parasail_result_is_rowcol(result)) {
        return result->rowcols->score_col;
    }
    return NULL; /* should not reach */
}

int* parasail_result_get_matches_col(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->matches_col;
}

int* parasail_result_get_similar_col(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->similar_col;
}

int* parasail_result_get_length_col(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_stats_rowcol(result));
    return result->stats->rowcols->length_col;
}

int* parasail_result_get_trace_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_trace(result));
    return result->trace->trace_table;
}

int* parasail_result_get_trace_ins_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_trace(result));
    return result->trace->trace_ins_table;
}

int* parasail_result_get_trace_del_table(const parasail_result_t * const restrict result)
{
    PARASAIL_CHECK_NULL(result);
    PARASAIL_ASSERT(parasail_result_is_trace(result));
    return result->trace->trace_del_table;
}

