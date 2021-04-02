/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 *
 * Helpers for implementing SG cases.
 */
#ifndef SG_HELPERS_H_
#define SG_HELPERS_H_

#ifdef PARASAIL_TABLE
#define SG_TABCOL _table
#else
#ifdef PARASAIL_ROWCOL
#define SG_TABCOL _rowcol
#else
#define SG_TABCOL
#endif
#endif

#ifdef SG_STATS
#define SG_STATSTRACE _stats
#else
#ifdef SG_TRACE
#define SG_STATSTRACE _trace
#else
#define SG_STATSTRACE
#endif
#endif

#ifndef SG_SUFFIX
#define SG_SUFFIX
#endif

#ifndef SG_SUFFIX_PROF
#define SG_SUFFIX_PROF
#endif

#define SG_FNAME1 parasail_sg_qb   
#define SG_FNAME2 parasail_sg_qe   
#define SG_FNAME3 parasail_sg_qx   
#define SG_FNAME4 parasail_sg_db   
#define SG_FNAME5 parasail_sg_de   
#define SG_FNAME6 parasail_sg_dx   
#define SG_FNAME7 parasail_sg_qb_de
#define SG_FNAME8 parasail_sg_qe_db
#define SG_FNAME9 parasail_sg      

#define SG_CASE1 1,0,0,0
#define SG_CASE2 0,1,0,0
#define SG_CASE3 1,1,0,0
#define SG_CASE4 0,0,1,0
#define SG_CASE5 0,0,0,1
#define SG_CASE6 0,0,1,1
#define SG_CASE7 1,0,0,1
#define SG_CASE8 0,1,1,0
#define SG_CASE9 1,1,1,1

#define CONCAT(A,B,C,D) A ## B ## C ## D

#define SG_IMPL(NAME, STATS, TABCOL, SUFFIX, CASE) \
parasail_result_t* CONCAT(NAME,STATS,TABCOL,SUFFIX) ( \
        const char * const restrict s1, const int s1Len, \
        const char * const restrict s2, const int s2Len, \
        const int open, const int gap, const parasail_matrix_t *matrix) \
{ \
    return FNAME(s1, s1Len, s2, s2Len, open, gap, matrix, CASE); \
}

#define SG_IMPL_ALL \
SG_IMPL(SG_FNAME1, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE1) \
SG_IMPL(SG_FNAME2, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE2) \
SG_IMPL(SG_FNAME3, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE3) \
SG_IMPL(SG_FNAME4, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE4) \
SG_IMPL(SG_FNAME5, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE5) \
SG_IMPL(SG_FNAME6, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE6) \
SG_IMPL(SG_FNAME7, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE7) \
SG_IMPL(SG_FNAME8, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE8) \
SG_IMPL(SG_FNAME9, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX, SG_CASE9)

#define SG_IMPL_PROF(NAME, STATS, TABCOL, SUFFIX, CASE) \
parasail_result_t* CONCAT(NAME,STATS,TABCOL,SUFFIX) ( \
		const parasail_profile_t * const restrict profile, \
		const char * const restrict s2, const int s2Len, \
		const int open, const int gap) \
{ \
    return PNAME(profile, s2, s2Len, open, gap, CASE); \
}

#define SG_IMPL_PROF_ALL \
SG_IMPL_PROF(SG_FNAME1, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE1) \
SG_IMPL_PROF(SG_FNAME2, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE2) \
SG_IMPL_PROF(SG_FNAME3, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE3) \
SG_IMPL_PROF(SG_FNAME4, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE4) \
SG_IMPL_PROF(SG_FNAME5, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE5) \
SG_IMPL_PROF(SG_FNAME6, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE6) \
SG_IMPL_PROF(SG_FNAME7, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE7) \
SG_IMPL_PROF(SG_FNAME8, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE8) \
SG_IMPL_PROF(SG_FNAME9, SG_STATSTRACE, SG_TABCOL, SG_SUFFIX_PROF, SG_CASE9)

#endif /* SG_HELPERS_H_ */
