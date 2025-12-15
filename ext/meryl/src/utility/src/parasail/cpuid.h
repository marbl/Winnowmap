/**
 * @file
 *
 * @author jeffrey.daily@gmail.com
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _PARASAIL_CPUID_H_
#define _PARASAIL_CPUID_H_

#include "parasail.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int parasail_can_use_avx512vbmi(void);
extern int parasail_can_use_avx512bw(void);
extern int parasail_can_use_avx512f(void);
extern int parasail_can_use_avx2(void);
extern int parasail_can_use_sse41(void);
extern int parasail_can_use_sse2(void);
extern int parasail_can_use_altivec(void);
extern int parasail_can_use_neon(void);

#ifdef __cplusplus
}
#endif

#endif /* _PARASAIL_CPUID_H_ */

