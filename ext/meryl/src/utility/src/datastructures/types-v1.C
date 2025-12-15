
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

#include <vector>
#include <tuple>

using merylutil::isEmptyString;

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
//  Convert a string to an integer type.
//
//  Decodes an integer constant in str[bgn...end-1] to a supplied 'result'
//  which is also returned from the function.
//
//  Whitespace at the begin or end of the string is ignored.
//
//  The constant can begin with minus sign in which case the value is negated
//  before being returned.
//
//  The constant can end in an optional single letter base indentifier:
//    'b' - base  2, binary
//    'o' - base  8, octal
//    'd' - base 10, decimal
//    'h' - base 16, hexadecimal
//  If no base is specified, decimal is assumed.
//
//  Some SI suffixes are also accepted instead of the base letter.  The
//  number is assumed to be decimal in this case.  A further 'i' can be added
//  to indicate a binary multiplier (2^10, etc).
//     k - 1e3    ki - 2^10   kilo
//     m - 1e6    mi - 2^20   mega
//     g - 1e9    gi - 2^30   giga
//     t - 1e12   ti - 2^40   tera
//     p - 1e15   pi - 2^50   peta
//     e - 1e18   ei - 2^60   exa
//  (see also https://physics.nist.gov/cuu/Units/binary.html)
//
//  If letters not in the base are encountered in the string, an error
//  message is pushed onto the err vector.
//
template<typename X>
X
decodeInteger(char const *str, uint64 bgn, uint64 end,
              X &result,
              std::vector<char const *> &err) {
  uint8  decode[256];
  bool   negate = false;
  uint64 scale  = 1000;

  //  Clear the result so we can return early.

  result = 0;

  //  Initialize the ASCII-to-value decoding table to all invalid values.

  for (uint32 qq=0; qq<256; qq++)
    decode[qq] = 0xff;

  //  We don't want to / can't use the normal decoding functions
  //  (strtouint64, atoi, strtoul, etc) as those will simply stop on the
  //  first non-decimal letter, and we want to declare that an error.
  //  So we rolled our own that also handles any base.
  //
  auto convertNumber = [&] (char const *str, uint64 b, uint64 e, uint8 shift) -> std::pair<uint64, bool> {
                         uint64  val = 0;
                         bool    inv = false;

                         for (uint64 ii=b; ii<e; ii++) {
                           uint8 num = decode[str[ii]];

                           if (num == 0xff)
                             inv = true;

                           val *= shift;
                           val += num;
                         }

                         return(std::make_pair(val, inv));
                       };

  //  Find the end?

  if (end == 0)
    for (end=bgn; str[end]; end++)
      ;

  //  Ignore spaces at the start and end.

  while ((bgn < end) && (isWhiteSpace(str[bgn]) == true))
    bgn++;
  while ((bgn < end) && (isWhiteSpace(str[end-1]) == true))
    end--;

  //  Remember and skip any negative sign.

  if (str[bgn] == '-') {
    negate = true;
    bgn++;
  }

  //  Remember and strip off the binary SI indicator.

  if (str[end-1] == 'i') {
    scale = 1024;
    end--;
  }

  //  Return if there is no number to decode.

  if (bgn == end)
    return(result);

  //  Find the base, decode the integer.  Or make errors.

  uint64      value = 0;
  bool        invalidNumber = false;
  char const *expectedType  = nullptr;

  if      (str[end-1] == 'b') {
    expectedType = "binary";

    decode['0'] = 0x00;   decode['1'] = 0x01;

    std::tie(value, invalidNumber) = convertNumber(str, bgn, end-1, 2);
  }

  else if (str[end-1] == 'o') {
    expectedType = "octal";

    decode['0'] = 0x00;   decode['1'] = 0x01;   decode['2'] = 0x02;   decode['3'] = 0x03;
    decode['4'] = 0x04;   decode['5'] = 0x05;   decode['6'] = 0x06;   decode['7'] = 0x07;

    std::tie(value, invalidNumber) = convertNumber(str, bgn, end-1, 8);
  }

  else if (str[end-1] == 'h') {
    expectedType = "hexadecimal";

    decode['0'] = 0x00;   decode['1'] = 0x01;   decode['2'] = 0x02;   decode['3'] = 0x03;   decode['4'] = 0x04;
    decode['5'] = 0x05;   decode['6'] = 0x06;   decode['7'] = 0x07;   decode['8'] = 0x08;   decode['9'] = 0x09;
    decode['a'] = 0x0a;   decode['b'] = 0x0b;   decode['c'] = 0x0c;   decode['d'] = 0x0d;   decode['e'] = 0x0e;   decode['f'] = 0x0f;
    decode['A'] = 0x0a;   decode['B'] = 0x0b;   decode['C'] = 0x0c;   decode['D'] = 0x0d;   decode['E'] = 0x0e;   decode['F'] = 0x0f;

    std::tie(value, invalidNumber) = convertNumber(str, bgn, end-1, 16);
  }

  else if ((str[end-1] == 'k') ||
           (str[end-1] == 'm') ||
           (str[end-1] == 'g') ||
           (str[end-1] == 't') ||
           (str[end-1] == 'p') ||
           (str[end-1] == 'e') ||
           (str[end-1] == 'd')) {
    expectedType = "decimal";

    decode['0'] = 0x00;   decode['1'] = 0x01;   decode['2'] = 0x02;   decode['3'] = 0x03;   decode['4'] = 0x04;
    decode['5'] = 0x05;   decode['6'] = 0x06;   decode['7'] = 0x07;   decode['8'] = 0x08;   decode['9'] = 0x09;

    std::tie(value, invalidNumber) = convertNumber(str, bgn, end-1, 10);

    switch (str[end-1]) {
      case 'e':   value *= scale;  [[fallthrough]];   //  Scale the value by 1000 or 1024 for
      case 'p':   value *= scale;  [[fallthrough]];   //  each multiplier increment.
      case 't':   value *= scale;  [[fallthrough]];   //
      case 'g':   value *= scale;  [[fallthrough]];   //  A 'g' will multiply by 'scale' three times.
      case 'm':   value *= scale;  [[fallthrough]];   //
      case 'k':   value *= scale;  [[fallthrough]];   //  (thank -Wimplicit-fallthrough for this)
      default:
        break;
    }
  }

  else if (('0' <= str[end-1]) &&
           (str[end-1] <= '9')) {
    expectedType = "decimal";

    decode['0'] = 0x00;   decode['1'] = 0x01;   decode['2'] = 0x02;   decode['3'] = 0x03;   decode['4'] = 0x04;
    decode['5'] = 0x05;   decode['6'] = 0x06;   decode['7'] = 0x07;   decode['8'] = 0x08;   decode['9'] = 0x09;

    std::tie(value, invalidNumber) = convertNumber(str, bgn, end, 10);
  }

  else {
    char *a    = new char [1024];
    char *b    = new char [1024];
    char *c    = new char [1024];
    int32 bpos = 0;

    snprintf(a, 1024, "Can't decode '%s': didn't find type of number ('b', 'o', 'd' or 'h') at end.", str);
    snprintf(b, 1024, "              ");
    c[0] = 0;

    bpos = strlen(b);

    for (uint64 ii=0; ii<bgn; ii++)
      b[bpos++] = ' ';

    for (uint64 ii=bgn; ii<end; ii++)
      b[bpos++] = '^';

    b[bpos] = 0;

    err.push_back(a);
    err.push_back(b);
    err.push_back(c);
  }

  if (invalidNumber == true) {
    char *a    = new char [1024];
    char *b    = new char [1024];
    char *c    = new char [1024];
    int32 bpos = 0;

    snprintf(a, 1024, "Can't decode '%s': %s number has invalid letters.", str, expectedType);
    snprintf(b, 1024, "              ");
    c[0] = 0;

    bpos = strlen(b);

    for (uint64 ii=0; ii<bgn; ii++)
      b[bpos++] = ' ';

    for (uint64 ii=bgn; ii<end; ii++)
      b[bpos++] = '^';

    b[bpos] = 0;

    err.push_back(a);
    err.push_back(b);
    err.push_back(c);
  }

  //  Copy the (negated) value to the result and return.

  if (invalidNumber)
    value = 0;

  if (negate)
    result = -value;
  else
    result =  value;

  return(result);
}

template  int8  decodeInteger(char const *str, uint64 bgn, uint64 end,  int8  &result, std::vector<char const *> &err);
template uint8  decodeInteger(char const *str, uint64 bgn, uint64 end, uint8  &result, std::vector<char const *> &err);
template  int16 decodeInteger(char const *str, uint64 bgn, uint64 end,  int16 &result, std::vector<char const *> &err);
template uint16 decodeInteger(char const *str, uint64 bgn, uint64 end, uint16 &result, std::vector<char const *> &err);
template  int32 decodeInteger(char const *str, uint64 bgn, uint64 end,  int32 &result, std::vector<char const *> &err);
template uint32 decodeInteger(char const *str, uint64 bgn, uint64 end, uint32 &result, std::vector<char const *> &err);
template  int64 decodeInteger(char const *str, uint64 bgn, uint64 end,  int64 &result, std::vector<char const *> &err);
template uint64 decodeInteger(char const *str, uint64 bgn, uint64 end, uint64 &result, std::vector<char const *> &err);


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

  //  We CANNOT unambiguously support empty start ranges unless
  //  unsigned integers are assumed.
  //    -123  -> negative 123?
  //    -123  -> range 0-123?
  //
  //if ((*ap == '-') ||              //  If this is a range but missing
  //    (*ap == '/')) {              //  the first number, set to min,
  //  bgn = std::numeric_limits<numberType>::min();
  //  if (*++ap != 0)                //  grab the second number and return.
  //    ap = strtonumber(ap+1, end);
  //  else
  //    end = std::numeric_limits<numberType>::max();
  //}

  ap = strtonumber(ap, bgn);       //  Grab the first number.

  end = bgn;                       //  Set the second to that.

  if ((*ap == '-') ||              //  If this is a range, or a
      (*ap == '/')) {              //  one-of-many selection, grab
    if (*++ap != 0)                //  the second number or set to
      ap = strtonumber(ap, end);   //  max if there isn't one.
    else
      end = std::numeric_limits<numberType>::max();
  }

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

char toBinDigit(uint8 value) { return alpha[value & 0x01]; }
char toOctDigit(uint8 value) { return alpha[value & 0x07]; }
char toDecDigit(uint8 value) { return alpha[value %   10]; }
char toHexDigit(uint8 value) { return alpha[value & 0x0f]; }


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
void
toBin(uintType v, char *&str, uint64 &strLen, uint64 &strMax, uint32 w) {
  merylutil::increaseArray(str, strLen + 129, strMax, 1024);
  strLen = toBin(v, str + strLen, w) - str;
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

template void        toBin<uint128>(uint128 v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toBin<uint64> (uint64  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toBin<uint32> (uint32  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toBin<uint16> (uint16  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toBin<uint8>  (uint8   v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);

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
void
toOct(uintType v, char *&str, uint64 &strLen, uint64 &strMax, uint32 w) {
  merylutil::increaseArray(str, strLen + 129, strMax, 1024);
  strLen = toOct(v, str + strLen, w) - str;
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

template void        toOct<uint128>(uint128 v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toOct<uint64> (uint64  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toOct<uint32> (uint32  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toOct<uint16> (uint16  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toOct<uint8>  (uint8   v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);

template char const *toOct<uint128>(uint128 v, uint32 width);
template char const *toOct<uint64> (uint64  v, uint32 width);
template char const *toOct<uint32> (uint32  v, uint32 width);
template char const *toOct<uint16> (uint16  v, uint32 width);
template char const *toOct<uint8>  (uint8   v, uint32 width);



template<typename intType>
char *
toDec(intType v, char *ret, uint32 w) {
  uint32   p = 1;
  uint32   e = 0;

  if (v < 0) {                            //  Add a negative sign, and make
    ret[0] = '-';                         //  the number positive.
    p++;
    v = -v;
  }

  for (intType t=v/10; t > 0; t /= 10)    //  Count how long the output string will
    p++;                                  //  be; we build backwards, right-to-left.

  ret[p]   =  0;                          //  Terminate the string, and remember the
  e = p;                                  //  end position so we can return it to the user.

  ret[--p] = alpha[ v % 10 ];             //  Convert the last digit; this handles v=0 too.

  for (intType t=v/10; t > 0; t /= 10)    //  Convert the next low order digit to
    ret[--p] = alpha[ t % 10 ];           //  an ASCII letter, repeat.

  if (ret[0] == '-')                      //  If a negative number, the last
    assert(p == 1);                       //  digit we added will be in [1],
  else                                    //  otherwise, it is in [0].
    assert(p == 0);

  return(ret + e);
}

template<typename intType>
void
toDec(intType v, char *&str, uint64 &strLen, uint64 &strMax, uint32 w) {
  merylutil::increaseArray(str, strLen + 129, strMax, 1024);
  strLen = toDec(v, str + strLen, w) - str;
}

template<typename intType>
char const *
toDec(intType v, uint32 w) {
  char *ret = getNextString();
  toDec(v, ret, w);
  return(ret);
}

template char       *toDec<uint128>(uint128 v, char *out, uint32 width);
template char       *toDec< int128>( int128 v, char *out, uint32 width);
template char       *toDec<uint64> (uint64  v, char *out, uint32 width);
template char       *toDec< int64> ( int64  v, char *out, uint32 width);
template char       *toDec<uint32> (uint32  v, char *out, uint32 width);
template char       *toDec< int32> ( int32  v, char *out, uint32 width);
template char       *toDec<uint16> (uint16  v, char *out, uint32 width);
template char       *toDec< int16> ( int16  v, char *out, uint32 width);
template char       *toDec<uint8>  (uint8   v, char *out, uint32 width);
template char       *toDec< int8>  ( int8   v, char *out, uint32 width);

template void        toDec<uint128>(uint128 v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec< int128>( int128 v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec<uint64> (uint64  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec< int64> ( int64  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec<uint32> (uint32  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec< int32> ( int32  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec<uint16> (uint16  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec< int16> ( int16  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec<uint8>  (uint8   v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toDec< int8>  ( int8   v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);

template char const *toDec<uint128>(uint128 v, uint32 width);
template char const *toDec< int128>( int128 v, uint32 width);
template char const *toDec<uint64> (uint64  v, uint32 width);
template char const *toDec< int64> ( int64  v, uint32 width);
template char const *toDec<uint32> (uint32  v, uint32 width);
template char const *toDec< int32> ( int32  v, uint32 width);
template char const *toDec<uint16> (uint16  v, uint32 width);
template char const *toDec< int16> ( int16  v, uint32 width);
template char const *toDec<uint8>  (uint8   v, uint32 width);
template char const *toDec< int8>  ( int8   v, uint32 width);



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
void
toHex(uintType v, char *&str, uint64 &strLen, uint64 &strMax, uint32 w) {
  merylutil::increaseArray(str, strLen + 129, strMax, 1024);
  strLen = toHex(v, str + strLen, w) - str;
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

template void        toHex<uint128>(uint128 v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toHex<uint64> (uint64  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toHex<uint32> (uint32  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toHex<uint16> (uint16  v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);
template void        toHex<uint8>  (uint8   v, char *&str, uint64 &strLen, uint64 &strMax, uint32 width);

template char const *toHex<uint128>(uint128 v, uint32 width=32);
template char const *toHex<uint64> (uint64  v, uint32 width=16);
template char const *toHex<uint32> (uint32  v, uint32 width=8);
template char const *toHex<uint16> (uint16  v, uint32 width=4);
template char const *toHex<uint8>  (uint8   v, uint32 width=2);

