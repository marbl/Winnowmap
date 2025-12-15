
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

#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstdint>
#include <cinttypes>
#include <climits>
#include <unistd.h>
#include <cstring>
#include <cctype>

#include <cfloat>
#include <cmath>

#include <cassert>
#include <cerrno>

#include <omp.h>

#include <limits>
#include <set>
#include <vector>

#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS == 32)
#error I do not support 32-bit off_t.
#endif

//  Make the basic int types a bit more friendly.

__extension__
typedef unsigned __int128  uint128;
__extension__
typedef          __int128   int128;

typedef uint64_t           uint64;
typedef  int64_t            int64;

typedef uint32_t           uint32;
typedef  int32_t            int32;

typedef uint16_t           uint16;
typedef  int16_t            int16;

typedef  uint8_t           uint8;
typedef   int8_t            int8;

//  There's no way to assign a constant value to the 128-bit integers
//  directly, but with a helper function we can assign it using two 64-bit
//  integers.  This only really makes sense for the unsigned flavor, e.g.,
//  when used for bit packed quantities.

constexpr inline uint128   build_uint128(uint64 a, uint64 b)    { return(((uint128)a << 64) | ((uint128)b)); }
constexpr inline  int128   build_int128 ( int64 a,  int64 b)    { return( ((int128)a << 64) |  ((int128)b)); }

//  Some handy constants.
//
//  numeric_limits<> on the 128-bit types is undefined, so we're forced to do
//  it the hard way.

constexpr uint128   uint128zero = 0;
constexpr uint128   uint128one  = 1;
constexpr uint128   uint128min  = 0;
constexpr uint128   uint128max  = (uint128)(0xffffffffffffffffllu) << 64 | (uint128)(0xffffffffffffffffllu);

constexpr  int128    int128zero = 0;
constexpr  int128    int128one  = 1;
constexpr  int128    int128min  = (uint128)(0x8000000000000000llu) << 64 | (uint128)(0x0000000000000000llu);
constexpr  int128    int128max  = (uint128)(0x7fffffffffffffffllu) << 64 | (uint128)(0xffffffffffffffffllu);

constexpr uint64    uint64zero  = 0;
constexpr uint64    uint64one   = 1;
constexpr uint64    uint64min   = std::numeric_limits<uint64>::min();
constexpr uint64    uint64max   = std::numeric_limits<uint64>::max();

constexpr  int64     int64zero  = 0;
constexpr  int64     int64one   = 1;
constexpr  int64     int64min   = std::numeric_limits<int64>::min();
constexpr  int64     int64max   = std::numeric_limits<int64>::max();

constexpr uint32    uint32zero  = 0;
constexpr uint32    uint32one   = 1;
constexpr uint32    uint32min   = std::numeric_limits<uint32>::min();
constexpr uint32    uint32max   = std::numeric_limits<uint32>::max();

constexpr  int32     int32zero  = 0;
constexpr  int32     int32one   = 1;
constexpr  int32     int32min   = std::numeric_limits<int32>::min();
constexpr  int32     int32max   = std::numeric_limits<int32>::max();

constexpr uint16    uint16zero  = 0;
constexpr uint16    uint16one   = 1;
constexpr uint16    uint16min   = std::numeric_limits<uint16>::min();
constexpr uint16    uint16max   = std::numeric_limits<uint16>::max();

constexpr  int16     int16zero  = 0;
constexpr  int16     int16one   = 1;
constexpr  int16     int16min   = std::numeric_limits<int16>::min();
constexpr  int16     int16max   = std::numeric_limits<int16>::max();

constexpr uint8     uint8zero   = 0;
constexpr uint8     uint8one    = 1;
constexpr uint8     uint8min    = std::numeric_limits<uint8>::min();
constexpr uint8     uint8max    = std::numeric_limits<uint8>::max();

constexpr  int8      int8zero   = 0;
constexpr  int8      int8one    = 1;
constexpr  int8      int8min    = std::numeric_limits<int8>::min();
constexpr  int8      int8max    = std::numeric_limits<int8>::max();

//  Conversion from floating point to integer.  lrint() rounds the
//  floating-point argument to an integer value, using the current rounding
//  mode.  This mode can be set with std::fesetround().

#if !defined(__clang__)
inline  int64  doubletoint64(double d)   { return(std::llrint(d)); }
inline  int32  doubletoint32(double d)   { return(std:: lrint(d)); }
#else
inline  int64  doubletoint64(double d)   { return(llrint(d)); }
inline  int32  doubletoint32(double d)   { return(lrint(d)); }
#endif

//  Decoding stings into numbers (and a boolean).
//   - The first set simply convert the string to a number and return that
//     number.
//   - The second set converts the string to a number and returns a pointer
//     to the letter in the string just after the number.
//
//  There probably should be a strtobool() of the second form, but I'm not
//  really sure what to do with the 'invalid' case that is currently treated
//  as 'false'.

uint128 strtoullll(char const *nptr, char **endptr);   //  The original strtoul() et al take char**
 int128 strtollll (char const *nptr, char **endptr);   //  as the second arg.

inline uint128 strtouint128(char const *str)  {  return((uint128)strtoullll(str, nullptr));      }
inline  int128 strtoint128 (char const *str)  {  return( (int128)strtollll (str, nullptr));      }
inline uint64  strtouint64 (char const *str)  {  return((uint64) strtoull  (str, nullptr, 10));  }
inline  int64  strtoint64  (char const *str)  {  return( (int64) strtoll   (str, nullptr, 10));  }
inline uint32  strtouint32 (char const *str)  {  return((uint32) strtoul   (str, nullptr, 10));  }
inline  int32  strtoint32  (char const *str)  {  return( (int32) strtol    (str, nullptr, 10));  }
inline uint16  strtouint16 (char const *str)  {  return((uint16) strtoul   (str, nullptr, 10));  }   //  WARNING: these convert to
inline  int16  strtoint16  (char const *str)  {  return( (int16) strtol    (str, nullptr, 10));  }   //  a 32-bit integer then cast
inline uint8   strtouint8  (char const *str)  {  return((uint8)  strtoul   (str, nullptr, 10));  }   //  to 16- or 8-bit integers.
inline  int8   strtoint8   (char const *str)  {  return( (int8)  strtol    (str, nullptr, 10));  }
inline  float  strtofloat  (char const *str)  {  return( (float) strtof    (str, nullptr));      }
inline double  strtodouble (char const *str)  {  return((double) strtod    (str, nullptr));      }

inline char const *strtonumber(char const *str, uint128 &num)  {  char *rem=nullptr;  num = (uint128)strtoullll(str, &rem);      return(rem);  }
inline char const *strtonumber(char const *str,  int128 &num)  {  char *rem=nullptr;  num =  (int128)strtollll (str, &rem);      return(rem);  }
inline char const *strtonumber(char const *str, uint64  &num)  {  char *rem=nullptr;  num = (uint64) strtoull  (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str,  int64  &num)  {  char *rem=nullptr;  num =  (int64) strtoll   (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str, uint32  &num)  {  char *rem=nullptr;  num = (uint32) strtoul   (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str,  int32  &num)  {  char *rem=nullptr;  num =  (int32) strtol    (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str, uint16  &num)  {  char *rem=nullptr;  num = (uint16) strtoul   (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str,  int16  &num)  {  char *rem=nullptr;  num =  (int16) strtol    (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str, uint8   &num)  {  char *rem=nullptr;  num = (uint8)  strtoul   (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str,  int8   &num)  {  char *rem=nullptr;  num =  (int8)  strtol    (str, &rem, 10);  return(rem);  }
inline char const *strtonumber(char const *str,  float  &num)  {  char *rem=nullptr;  num = (double) strtof    (str, &rem);      return(rem);  }
inline char const *strtonumber(char const *str, double  &num)  {  char *rem=nullptr;  num = (double) strtod    (str, &rem);      return(rem);  }

inline bool    strtobool(char const *str)  {
  if ((str == nullptr) ||
      (str[0] == 0))
    return(false);

  if (((str[0] == 'y') && (str[1] == 0)) ||
      ((str[0] == 'Y') && (str[1] == 0)) ||
      ((str[0] == 't') && (str[1] == 0)) ||
      ((str[0] == 'T') && (str[1] == 0)) ||
      ((str[0] == '1') && (str[1] == 0)) ||
      ((str[0] == '+') && (str[1] == 0)))
    return(true);

  if ((strcasecmp(str, "yes")  == 0) ||
      (strcasecmp(str, "true") == 0))
    return(true);

  return(false);
}


//  A fancy string-to-integer conversion that handles positive and negative
//  numbers, different bases, SI units (and the binary form) and different
//  types.
//
//  Errors are reutned in a vector, which is convenient but only for command
//  line processing.
//
template<typename X>
X
decodeInteger(char const *str, uint64 bgn, uint64 end,
              X &result,
              std::vector<char const *> &err);


//  Test if a character or string is of the desired encoding.

#define inRange(A,c,B) ( ((A) <= c) && (c <= (B)) )

inline bool   isNUL(char c)        { return c == 0x00; }  //
inline bool   isBEL(char c)        { return c == 0x07; }  //  '\a'
inline bool   isBS(char c)         { return c == 0x08; }  //  '\b'
inline bool   isTab(char c)        { return c == 0x09; }  //  '\t'
inline bool   isLF(char c)         { return c == 0x0a; }  //  '\n'
inline bool   isVT(char c)         { return c == 0x0b; }  //  '\v', maybe
inline bool   isFF(char c)         { return c == 0x0c; }  //  '\f'
inline bool   isCR(char c)         { return c == 0x0d; }  //  '\r'
inline bool   isSpace(char c)      { return c == 0x20; }  //  ' '

inline bool   isEndOfLine(char c)  { return isLF(c)       //  Most common.
                                         || isCR(c); }
inline bool   isWhiteSpace(char c) { return isSpace(c)    //  Most common.
                                         || isTab(c)      //
                                         || isLF(c)       //
                                         || isCR(c)       //
                                         || isFF(c)       //
                                         || isVT(c); }    //  Least common.
inline bool   isVisible(char c)    { return inRange('!', c, '~'); }

inline bool   isLetter(char c)     { return inRange('a', c, 'z')
                                         || inRange('A', c, 'Z');  }
inline bool   isBinDigit(char c)   { return inRange('0', c, '1');  }
inline bool   isOctDigit(char c)   { return inRange('0', c, '7');  }
inline bool   isDecDigit(char c)   { return inRange('0', c, '9');  }
inline bool   isHexDigit(char c)   { return inRange('0', c, '9')
                                         || inRange('a', c, 'f')
                                         || inRange('A', c, 'F');  }

inline bool   isComment(char c)    { return (c == '!') || (c == '#') || (c == 0); }
inline bool   isDelimiter(char c)  { return (c == ':') || (c == '=');             }

#undef inRange

bool          isBinNumber (char const *s);
bool          isOctNumber (char const *s);
bool          isDecNumber (char const *s, char dot='.');
bool inline   isDecInteger(char const *s)   { return(isDecNumber(s,  0));  }
bool inline   isDecFloat  (char const *s)   { return(isDecNumber(s, '.')); }
bool          isHexNumber (char const *s);

//  Disallow the usual character tests becuse of their goofy return values.

#undef  isalnum
#undef  isalpha
#undef  iscntrl
#undef  isdigit
#undef  isgraph
#undef  islower
#undef  isprint
#undef  ispunct
#undef  isspace
#undef  isupper
#undef  isxdigit
#undef  isnumber

int inline isalnum (char c) = delete;
int inline isalpha (char c) = delete;
int inline iscntrl (char c) = delete;
int inline isdigit (char c) = delete;
int inline isgraph (char c) = delete;
int inline islower (char c) = delete;
int inline isprint (char c) = delete;
int inline ispunct (char c) = delete;
int inline isspace (char c) = delete;
int inline isupper (char c) = delete;
int inline isxdigit(char c) = delete;
int inline isnumber(char c) = delete;

//  Convert an ascii binary, octal, decimal or hexadecimal letter to an
//  integer.  No type checking is performed; you've already called
//  isHexNumber() et al, right?  Inputs not in the type return undefined
//  values -- asciiBinToInteger('2'), asciiOctToInteger('9'),
//  asciiHexToInteger('g') are undefined.
//
//  The pieces of asciiHexToInteger() are as follows:
//    (d & 0xf)        //  Decodes '0'-'9' as 0-9, 'a' - 'f' as 1-6
//    (d >> 6)         //  Decodes digits as 0, letters as 1.
//   ((d >> 6) << 3)   //  Decodes digits as 0, letters as 8.
//      the above was for an older implementation, but the
//      same idea applies the current implementation.

inline uint8  asciiBinToInteger(char d)   { return d - '0'; }   //  Pretty trivial.
inline uint8  asciiOctToInteger(char d)   { return d - '0'; }
inline uint8  asciiDecToInteger(char d)   { return d - '0'; }
inline uint8  asciiHexToInteger(char d)   {
  uint8  n = (uint8)d      & 0xf;  //  [0-9] -> 0-9; [AaQq] -> 1, [BbRr] -> 2, ..., [FfVv] -> 6, [GgWw] -> 7, ..., [Oo_<del>] -> f
  uint8  s = (uint8)d >> 6 & 0x1;  //  [@-\xff] -> 1

  return (n + 9 * s) & 0xf;
}

//  Convert an integer to a printable letter.  If it's not a printable
//  letter, returns '.'.

inline
char
integerToLetter(uint32 i) {
  return(((' ' <= i) && (i <= '~')) ? i : '.');
}

//  Convert a string representing a set of numbers to
//   - the first and last values (for form '#' or '#-#')
//   - a vector of the low and high values
//   - a set of the values
//
//  The string should be comprised of multiple comma separated ranges:
//   - #     a single number
//   - #-#   a range of numbers
//   - #/#   a one-out-of-N specification
//
//  The first form returns a pointer to the letter after the decoded values.
//
//  If a single number is encountered in the first or second forms, both
//  'bgn' and 'end' are set to that value.
//
//  If 'numberType' is a 128-bit integer, only 64-bit integers can be
//  converted.

template<typename numberType> char const *decodeRange(char const *range, numberType &bgn, numberType &end);
template<typename numberType> void        decodeRange(char const *range, std::vector<numberType> &bgn, std::vector<numberType> &end);
template<typename numberType> void        decodeRange(char const *range, std::set<numberType> &values);

//  Convert an unsigned integer representing bits or bytes to
//  a floating point number representing GB or MB.

inline double bitsToGB(uint64 bits)   { return(bits / 8 / 1024.0 / 1024.0 / 1024.0); }
inline double bitsToMB(uint64 bits)   { return(bits / 8 / 1024.0 / 1024.0);          }

//  Convert an unsigned integer to one with 3 significant digit number, and
//  also return the correct SI base.

uint64      scaledNumber(uint64 n, uint32 div=1024);   //  Return n between 0 and div,
char        scaledUnit  (uint64 n, uint32 div=1024);   //  and the SI unit of that
const char *scaledName  (uint64 n, uint32 div=1024);   //  scaling.

//  Convert an unsigned integer to a character string in the desired base.
//
//    char *toXXX(v, str)
//      Expects a pre-allocated character buffer 'str' with enough space for
//      the output string and a NUL terminating byte.  It returns a pointer
//      to the NUL byte.  A 128-bit integer in:
//        binary      needs 129 bytes
//        octal       needs  44 bytes
//        decimal     needs  40 bytes (it's 340,282,366,920,938,463,463,374,607,431,768,211,455)
//        hexadecimal needs  33 bytes
//
//    char const *toXX(v)
//      Returns a pointer to one of 32 private string buffers.  This is
//      thread safe, as long as you don't use it more than 32 times at once.
//
//  Both forms take an optional 'width' (in bits) to display.  The actual
//  width used is the minimum of this width and the number of bits in the
//  type.  toDec() accepts the width, but doesn't use it.

char toBinDigit(uint8 value);
char toOctDigit(uint8 value);
char toDecDigit(uint8 value);
char toHexDigit(uint8 value);

template<typename uintType> char       *toBin(uintType value, char *out, uint32 width=128);
template<typename uintType> char       *toOct(uintType value, char *out, uint32 width=128);
template<typename  intType> char       *toDec( intType value, char *out, uint32 width=128);
template<typename uintType> char       *toHex(uintType value, char *out, uint32 width=128);

template<typename uintType> void        toBin(uintType value, char *&str, uint64 &strLen, uint64 &strMax, uint32 width=128);
template<typename uintType> void        toOct(uintType value, char *&str, uint64 &strLen, uint64 &strMax, uint32 width=128);
template<typename  intType> void        toDec( intType value, char *&str, uint64 &strLen, uint64 &strMax, uint32 width=128);
template<typename uintType> void        toHex(uintType value, char *&str, uint64 &strLen, uint64 &strMax, uint32 width=128);

template<typename uintType> char const *toBin(uintType value, uint32 width=128);
template<typename uintType> char const *toOct(uintType value, uint32 width=128);
template<typename  intType> char const *toDec( intType value, uint32 width=128);
template<typename uintType> char const *toHex(uintType value, uint32 width=128);

//  Format specifications for printf()

#define F_PTR    "0x%016p"   // Pointers
#define F_C           "%c"   // Characters
#define F_CP           "c"
#define F_CI         "%*c"
#define F_STR         "%s"   // Strings
#define F_STRP         "s"
#define F_STRI       "%*s"
#define F_S16    "%" PRId16  // Integers
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
#define F_F32         "%f"   // Floating points
#define F_F32P         "f"
#define F_F32I       "%*f"
#define F_F64        "%lf"
#define F_F64P        "lf"
#define F_F64I      "%*lf"
#define F_SIZE_T     "%zu"   // Standard typedefs
#define F_SIZE_TP     "zu"
#define F_SIZE_TI   "%*zu"
#define F_OFF_T     F_S64
#define F_OFF_TP    F_S64P
#define F_OFF_TI    F_S64I


#endif  //  TYPES_H
