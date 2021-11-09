
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

#include "types.H"
#include "strings.H"


////////////////////////////////////////////////////////////
//
//  Sadly, there is no equivalent of strtoull() for 128-bit integers, so I
//  provide my own.  Only base 10 is supported.  Overflow isn't handled.
//
//  The obvious implementation of strtollll() -- that being to sum all the
//  digits and then negate the sum for negative values -- doesn't handle
//  int128min technically correct.  It ends up overflowing the (positive)
//  int128 by one.  This results in int128min.  Fortunately, int128min ==
//  -int128min, and the negation done by 'neg' doesn't do anything.  As
//  implemented below, though, we instead add or subtract each new digit,
//  giving us overflow-free results (but a little slower).
//

uint128
strtoullll(char const *nptr, char **endptr) {
  uint128     res = 0;
  char const *ptr = nptr;

  if (isEmptyString(ptr))
    return(res);

  while ((*ptr != 0) && (isWhiteSpace(*ptr) == true))
    ptr++;

  while ((*ptr != 0) && (isDecDigit(*ptr) == true)) {
    res *= 10;
    res += asciiDecToInteger(*ptr);

    ptr++;
  }

  if (endptr)
    *endptr = (char *)ptr;
  return(res);
}

int128
strtollll(char const *nptr, char **endptr) {
  int128      res = 0;
  bool        neg = false;
  char const *ptr = nptr;

  if (isEmptyString(ptr))
    return(res);

  while ((*ptr != 0) && (isWhiteSpace(*ptr) == true))
    ptr++;

  switch (*ptr) {
    case '-':  ptr++;  neg = true;  break;
    case '+':  ptr++;               break;
    default:                        break;
  }

  while ((*ptr != 0) && (isDecDigit(*ptr) == true)) {
    res *= 10;

    if (neg == false)
      res += asciiDecToInteger(*ptr);
    else
      res -= asciiDecToInteger(*ptr);

    ptr++;
  }

  if (endptr)
    *endptr = (char *)ptr;
  return(res);
}



////////////////////////////////////////////////////////////
//
//  Test if a string is a number in the desired encoding.
//

bool
isBinNumber(char const *s) {
  if (isEmptyString(s) == true)
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if (isBinDigit(s[ii]) == false)
      return(false);

  return(true);
}


bool
isOctNumber(char const *s) {
  if (isEmptyString(s) == true)
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if (isOctDigit(s[ii]) == false)
      return(false);

  return(true);
}


bool
isDecNumber(char const *s, char dot) {
  if (isEmptyString(s) == true)
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if ((isDecDigit(s[ii]) == false) && (s[ii] != dot))
      return(false);

  return(true);
}


bool
isHexNumber(char const *s) {
  if (isEmptyString(s) == true)
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if (isHexDigit(s[ii]) == false)
      return(false);

  return(true);
}



////////////////////////////////////////////////////////////
//
//  Convert a string of numbers to a pair of numbers, a vector of ranges, or
//  a set of values.
//

template<typename numberType>
char const *
decodeRange(char const *range, numberType &bgn, numberType &end) {
  char const    *ap = range;

  ap = strtonumber(ap, bgn);       //  Grab the first number.

  end = bgn;                       //  Set the second to that.

  if ((*ap == '-') ||              //  If this is a range,
      (*ap == '/'))                //  or a one-of-many selection,
    ap = strtonumber(ap+1, end);   //  grab the second number

  if (*ap == ',')                  //  If the next letter continues
    return(ap + 1);                //  move past that and return.

  if (*ap == 0)                    //  If the next letter is the end
    return(nullptr);               //  of the string, return nullptr.

  //  Otherwise, we can't decode this range.

  fprintf(stderr, "ERROR: invalid range '%s'\n", range);
  exit(1);

  return(nullptr);
}


template<typename numberType>
void
decodeRange(char const *range, std::vector<numberType> &bgn, std::vector<numberType> &end) {
  char const  *ap = range;
  numberType   av = 0;
  numberType   bv = 0;

  while (isEmptyString(ap) == false) {
    ap = decodeRange(ap, av, bv);

    bgn.push_back(av);
    end.push_back(bv);
  }
}


template<typename numberType>
void
decodeRange(char const *range, std::set<numberType> &values) {
  char const  *ap = range;
  numberType   av = 0;
  numberType   bv = 0;

  while (isEmptyString(ap) == false) {
    ap = decodeRange(ap, av, bv);

    for (numberType xx=av; xx<=bv; xx++)
      values.insert(xx);
  }
}


template  char const *decodeRange<uint128>(char const *range, uint128 &bgn, uint128 &end);
template  char const *decodeRange <int128>(char const *range,  int128 &bgn,  int128 &end);
template  char const *decodeRange<uint64> (char const *range, uint64  &bgn, uint64  &end);
template  char const *decodeRange <int64> (char const *range,  int64  &bgn,  int64  &end);
template  char const *decodeRange<uint32> (char const *range, uint32  &bgn, uint32  &end);
template  char const *decodeRange <int32> (char const *range,  int32  &bgn,  int32  &end);
template  char const *decodeRange<uint16> (char const *range, uint16  &bgn, uint16  &end);
template  char const *decodeRange <int16> (char const *range,  int16  &bgn,  int16  &end);
template  char const *decodeRange<uint8>  (char const *range, uint8   &bgn, uint8   &end);
template  char const *decodeRange <int8>  (char const *range,  int8   &bgn,  int8   &end);
template  char const *decodeRange<double> (char const *range, double  &bgn, double  &end);

template  void  decodeRange<uint128>(char const *range, std::vector<uint128> &bgn, std::vector<uint128> &end);
template  void  decodeRange <int128>(char const *range, std::vector <int128> &bgn, std::vector <int128> &end);
template  void  decodeRange<uint64> (char const *range, std::vector<uint64>  &bgn, std::vector<uint64>  &end);
template  void  decodeRange <int64> (char const *range, std::vector <int64>  &bgn, std::vector <int64>  &end);
template  void  decodeRange<uint32> (char const *range, std::vector<uint32>  &bgn, std::vector<uint32>  &end);
template  void  decodeRange <int32> (char const *range, std::vector <int32>  &bgn, std::vector <int32>  &end);
template  void  decodeRange<uint16> (char const *range, std::vector<uint16>  &bgn, std::vector<uint16>  &end);
template  void  decodeRange <int16> (char const *range, std::vector <int16>  &bgn, std::vector <int16>  &end);
template  void  decodeRange<uint8>  (char const *range, std::vector<uint8>   &bgn, std::vector<uint8>   &end);
template  void  decodeRange <int8>  (char const *range, std::vector <int8>   &bgn, std::vector <int8>   &end);
template  void  decodeRange <double>(char const *range, std::vector <double> &bgn, std::vector <double> &end);

template  void  decodeRange<uint128>(char const *range, std::set<uint128> &values);
template  void  decodeRange <int128>(char const *range, std::set <int128> &values);
template  void  decodeRange<uint64> (char const *range, std::set<uint64>  &values);
template  void  decodeRange <int64> (char const *range, std::set <int64>  &values);
template  void  decodeRange<uint32> (char const *range, std::set<uint32>  &values);
template  void  decodeRange <int32> (char const *range, std::set <int32>  &values);
template  void  decodeRange<uint16> (char const *range, std::set<uint16>  &values);
template  void  decodeRange <int16> (char const *range, std::set <int16>  &values);
template  void  decodeRange<uint8>  (char const *range, std::set<uint8>   &values);
template  void  decodeRange <int8>  (char const *range, std::set <int8>   &values);
template  void  decodeRange <double>(char const *range, std::set <double> &values);



////////////////////////////////////////////////////////////
//
//  Convert an unsigned integer to one with 3 significant digit number, and
//  also return the correct SI base.
//
//  This does NOT round correctly.  We'd need to track the remainder
//  and increment 'n' if the remainder is more than half 'div'.
//

uint64
scaledNumber(uint64 n, uint32 div) {

  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;
  if (n > 9999)   n /= div;

  return(n);
}

char
scaledUnit(uint64 n, uint32 div) {
  char u = ' ';

  if (n > 9999)  {  n /= div; u = 'k';  }   //  kilo
  if (n > 9999)  {  n /= div; u = 'M';  }   //  mega
  if (n > 9999)  {  n /= div; u = 'G';  }   //  giga
  if (n > 9999)  {  n /= div; u = 'T';  }   //  tera
  if (n > 9999)  {  n /= div; u = 'P';  }   //  peta
  if (n > 9999)  {  n /= div; u = 'E';  }   //  exa
  if (n > 9999)  {  n /= div; u = 'Z';  }   //  zetta
  if (n > 9999)  {  n /= div; u = 'Y';  }   //  yotta

  return(u);
}

const char *
scaledName(uint64 n, uint32 div) {
  const char *s = "";

  if (n > 9999)  {  n /= div; s = " thousand";     }
  if (n > 9999)  {  n /= div; s = " million";      }
  if (n > 9999)  {  n /= div; s = " billion";      }
  if (n > 9999)  {  n /= div; s = " trillion";     }
  if (n > 9999)  {  n /= div; s = " quadrillion";  }
  if (n > 9999)  {  n /= div; s = " quintillion";  }
  if (n > 9999)  {  n /= div; s = " sextillion";   }
  if (n > 9999)  {  n /= div; s = " septillion";   }

  return(s);
}



////////////////////////////////////////////////////////////
//
//  Convert an unsigned integer to a character string in the desired base.
//
//  All follow the same pattern except for the constants, and except for
//  toDec() which also differs in the 'shift' operation.
//
//  The helper function getNextString() is the only part that needs
//  to worry about thread safety.  Everything else operates on that
//  returned buffer space.
//
//  Instead of allocating 32 strings of max length, we could allocate 4096
//  bytes and dole pieces out of the appropriate max length as needed.
//

char     alpha[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
char    *strAlloc  =   nullptr;
char    *str[32]   = { nullptr };
uint32   pos       =   0;

static
char *                   //  Helper function to return the next available buffer.
getNextString(void) {    //  This is the only part that needs to care about threads.
  char *ret = nullptr;   //  Everything else operates on the buffer returned.

#pragma omp critical (toHEXlock)
  {
    if (strAlloc == nullptr) {
      strAlloc = new char [32 * 129];

      for (uint32 ii=0; ii<32; ii++)
        str[ii] = strAlloc + ii * 129;
    }

    ret = str[pos++];

    if (pos >= 32)
      pos = 0;
  }

  return(ret);
}



template<typename uintType>
char *
toBin(uintType v, char *ret, uint32 w) {
  uint32  L = w;
  uint32  W = sizeof(uintType) * 8;
  uint32  l = std::min(L, W);
  uint32  p = l;
  uint32  s = 0;

  while (p > 0) {
    p -= 1;
    ret[p] = alpha[ (v >> s) & 0x1 ];
    s += 1;
  }

  ret[w] = 0;

  return(ret + w);
}

template<typename uintType>
char const *
toBin(uintType v, uint32 w) {
  char *ret = getNextString();
  toBin(v, ret, w);
  return(ret);
}

template char       *toBin<uint128>(uint128 v, char *out, uint32 width);
template char       *toBin<uint64> (uint64  v, char *out, uint32 width);
template char       *toBin<uint32> (uint32  v, char *out, uint32 width);
template char       *toBin<uint16> (uint16  v, char *out, uint32 width);
template char       *toBin<uint8>  (uint8   v, char *out, uint32 width);

template char const *toBin<uint128>(uint128 v, uint32 width);
template char const *toBin<uint64> (uint64  v, uint32 width);
template char const *toBin<uint32> (uint32  v, uint32 width);
template char const *toBin<uint16> (uint16  v, uint32 width);
template char const *toBin<uint8>  (uint8   v, uint32 width);



template<typename uintType>
char *
toOct(uintType v, char *ret, uint32 w) {
  uint32  L = (w + 2) / 3;
  uint32  W = sizeof(uintType) * 8 / 3 + 1;
  uint32  l = std::min(L, W);
  uint32  p = l;
  uint32  s = 0;

  while (p > 0) {
    p -= 1;
    ret[p] = alpha[ (v >> s) & 0x7 ];
    s += 3;
  }

  ret[l] = 0;

  return(ret + l);
}

template<typename uintType>
char const *
toOct(uintType v, uint32 w) {
  char *ret = getNextString();
  toOct(v, ret, w);
  return(ret);
}

template char       *toOct<uint128>(uint128 v, char *out, uint32 width);
template char       *toOct<uint64> (uint64  v, char *out, uint32 width);
template char       *toOct<uint32> (uint32  v, char *out, uint32 width);
template char       *toOct<uint16> (uint16  v, char *out, uint32 width);
template char       *toOct<uint8>  (uint8   v, char *out, uint32 width);

template char const *toOct<uint128>(uint128 v, uint32 width);
template char const *toOct<uint64> (uint64  v, uint32 width);
template char const *toOct<uint32> (uint32  v, uint32 width);
template char const *toOct<uint16> (uint16  v, uint32 width);
template char const *toOct<uint8>  (uint8   v, uint32 width);



template<typename uintType>
char *
toDec(uintType v, char *ret, uint32 w) {
  uint32   p = 1;
  uint32   e = 0;

  for (uintType t=v/10; t > 0; t /= 10)   //  Count how long the output string will
    p++;                                  //  be; we build backwards, right-to-left.

  ret[p]   =  0;                          //  Terminate the string, and remember the
  e = p;                                  //  end position so we can return it to the user.

  ret[--p] = alpha[ v % 10 ];             //  Convert the last digit; this handles v=0 too.

  for (uintType t=v/10; t > 0; t /= 10)   //  Convert the next low order digit to
    ret[--p] = alpha[ t % 10 ];           //  an ASCII letter, repeat.

  assert(p == 0);

  return(ret + e);
}

template<typename uintType>
char const *
toDec(uintType v, uint32 w) {
  char *ret = getNextString();
  toDec(v, ret, w);
  return(ret);
}

template char       *toDec<uint128>(uint128 v, char *out, uint32 width);
template char       *toDec<uint64> (uint64  v, char *out, uint32 width);
template char       *toDec<uint32> (uint32  v, char *out, uint32 width);
template char       *toDec<uint16> (uint16  v, char *out, uint32 width);
template char       *toDec<uint8>  (uint8   v, char *out, uint32 width);

template char const *toDec<uint128>(uint128 v, uint32 width);
template char const *toDec<uint64> (uint64  v, uint32 width);
template char const *toDec<uint32> (uint32  v, uint32 width);
template char const *toDec<uint16> (uint16  v, uint32 width);
template char const *toDec<uint8>  (uint8   v, uint32 width);



template<typename uintType>
char *
toHex(uintType v, char *ret, uint32 w) {
  uint32  L = sizeof(uintType) * 8 / 4;   //  The maximum possible width
  uint32  W = (w + 3) / 4;                //  The user suggested width
  uint32  l = std::min(L, W);
  uint32  p = l;
  uint32  s = 0;

  while (p > 0) {
    p -= 1;
    ret[p] = alpha[ (v >> s) & 0xf ];
    s += 4;
  }

  ret[l] = 0;

  return(ret + l);
}

template<typename uintType>
char const *
toHex(uintType v, uint32 w) {
  char *ret = getNextString();
  toHex(v, ret, w);
  return(ret);
}

template char       *toHex<uint128>(uint128 v, char *out, uint32 width);
template char       *toHex<uint64> (uint64  v, char *out, uint32 width);
template char       *toHex<uint32> (uint32  v, char *out, uint32 width);
template char       *toHex<uint16> (uint16  v, char *out, uint32 width);
template char       *toHex<uint8>  (uint8   v, char *out, uint32 width);

template char const *toHex<uint128>(uint128 v, uint32 width=32);
template char const *toHex<uint64> (uint64  v, uint32 width=16);
template char const *toHex<uint32> (uint32  v, uint32 width=8);
template char const *toHex<uint16> (uint16  v, uint32 width=4);
template char const *toHex<uint8>  (uint8   v, uint32 width=2);

