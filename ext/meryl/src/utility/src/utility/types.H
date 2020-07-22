
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
 *
 *  This software is based on:
 *    'Canu' v2.0              (https://github.com/marbl/canu)
 *  which is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef TYPES_H
#define TYPES_H

//  Do NOT enable _GLIBCXX_PARALLEL.  Performance is atrocious, and it will
//  not sort in-place, so memory will blow up.  More comments in AS_global.C.
#undef  _GLIBCXX_PARALLEL

//  ISO C99 says that to get INT32_MAX et al, these must be defined.
//  (7.18.2, 7.18.4, 7.8.1)
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

//  Tell gcc (and others, maybe) about unused parameters.  This is important for gcc (especially
//  newer ones) that complain about unused parameters.  Thanks to ideasman42 at
//  http://stackoverflow.com/questions/3599160/unused-parameter-warnings-in-c-code.
#ifdef __GNUC__
#define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#define UNUSED(x) UNUSED_ ## x
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include <float.h>
#include <math.h>

#include <assert.h>
#include <errno.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <omp.h>

#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS == 32)
#error I do not support 32-bit off_t.
#endif

typedef int8_t  int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef __int128            int128;
typedef unsigned __int128  uint128;

#define  uint128NUMBER(A,B)   ((((uint128)A) << 64) | ((uint128)B))
#define  uint64NUMBER(X) X ## LLU
#define  uint32NUMBER(X) X ## LU

#define  uint64ZERO      uint64NUMBER(0x0000000000000000)
#define  uint64ONE       uint64NUMBER(0x0000000000000001)
#define  uint64MAX       uint64NUMBER(0xffffffffffffffff)
#define  uint64MASK(X)   (((~uint64ZERO) >> (64 - (X))) & (-(uint64)((X) != 0)))

#define  uint32ZERO      uint32NUMBER(0x00000000)
#define  uint32ONE       uint32NUMBER(0x00000001)
#define  uint32MAX       uint32NUMBER(0xffffffff)
#define  uint32MASK(X)   (((~uint32ZERO) >> (32 - (X))) & (-(uint32)((X) != 0)))

#define  uint16ZERO      (0x0000)
#define  uint16ONE       (0x0001)
#define  uint16MAX       (0xffff)
#define  uint16MASK(X)   (((~uint16ZERO) >> (16 - (X))) & (-(uint16)((X) != 0)))

#define  uint8ZERO       (0x00)
#define  uint8ONE        (0x01)
#define  uint8MAX        (0xff)
#define  uint8MASK(X)    (((~uint8ZERO) >> (8 - (X))) & (-(uint8)((X) != 0)))


inline  int32  strtoint32 (char *str)  {  return( (int32)strtol  (str, NULL, 10));  }
inline uint32  strtouint32(char *str)  {  return((uint32)strtoul (str, NULL, 10));  }
inline  int64  strtoint64 (char *str)  {  return( (int64)strtoll (str, NULL, 10));  }
inline uint64  strtouint64(char *str)  {  return((uint64)strtoull(str, NULL, 10));  }
inline double  strtodouble(char *str)  {  return((double)strtod  (str, NULL));      }

inline void    strtonumber( int32 &num, char *str, char **rem)  {  num =  (int32)strtol  (str, rem, 10);  }
inline void    strtonumber(uint32 &num, char *str, char **rem)  {  num = (uint32)strtoul (str, rem, 10);  }
inline void    strtonumber( int64 &num, char *str, char **rem)  {  num =  (int64)strtoll (str, rem, 10);  }
inline void    strtonumber(uint64 &num, char *str, char **rem)  {  num = (uint64)strtoull(str, rem, 10);  }
inline void    strtonumber(double &num, char *str, char **rem)  {  num = (double)strtod  (str, rem);      }

inline  int32  doubletoint32(double d)   { return((int32) ((1.0 + 16.0 * DBL_EPSILON) * d));  };
inline  int64  doubletoint64(double d)   { return((int64) ((1.0 + 16.0 * DBL_EPSILON) * d));  };

template<typename uintType> char const *toHex(uintType v);
template<typename uintType> char const *toDec(uintType v);

// These macros are use to eliminate inter-platform differnces between
// calculated results
//#define DBL_TO_INT(X)   ((int)((1.0+16.0*DBL_EPSILON)*(X)))
//#define ROUNDPOS(X)     (DBL_TO_INT((X)+0.5) )
//#define ROUND(X)        (((X)>0.0) ? ROUNDPOS(X) : -ROUNDPOS(-(X)) )
//#define ZERO_PLUS       ( 16.0*DBL_EPSILON)
//#define ZERO_MINUS      (-16.0*DBL_EPSILON)
//#define ONE_PLUS        (1.0+ZERO_PLUS)
//#define ONE_MINUS       (1.0+ZERO_MINUS)
//#define INT_EQ_DBL(I,D) (fabs((double)(I)-(D)) < 16.0*DBL_EPSILON  )
//#define DBL_EQ_DBL(A,B) (fabs((A)-(B))<16.0*DBL_EPSILON)

// Pointers
#define F_PTR    "0x%016p"

// Characters
#define F_C           "%c"
#define F_CP           "c"
#define F_CI         "%*c"

// Strings
#define F_STR         "%s"
#define F_STRP         "s"
#define F_STRI       "%*s"

// Integers
#define F_S16    "%" PRId16
#define F_S16P       PRId16
#define F_S16I  "%*" PRId16
#define F_U16    "%" PRIu16
#define F_U16P       PRIu16
#define F_U16I  "%*" PRIu16
#define F_S32    "%" PRId32
#define F_S32P       PRId32
#define F_S32I  "%*" PRId32
#define F_U32    "%" PRIu32
#define F_U32P       PRIu32
#define F_U32I  "%*" PRIu32
#define F_S64    "%" PRId64
#define F_S64P       PRId64
#define F_S64I  "%*" PRId64
#define F_U64    "%" PRIu64
#define F_U64P       PRIu64
#define F_U64I  "%*" PRIu64
#define F_X64 "%016" PRIx64
#define F_X64P       PRIx64
#define F_X64I  "%*" PRIx64

// Floating points
#define F_F32         "%f"
#define F_F32P         "f"
#define F_F32I       "%*f"
#define F_F64        "%lf"
#define F_F64P        "lf"
#define F_F64I      "%*lf"

// Standard typedefs
#define F_SIZE_T     "%zu"
#define F_SIZE_TP     "zu"
#define F_SIZE_TI   "%*zu"

#define F_OFF_T     F_S64
#define F_OFF_TP    F_S64P
#define F_OFF_TI    F_S64I


#endif  //  TYPES_H
