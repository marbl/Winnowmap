
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

#include "strings.H"

namespace merylutil::inline strings::inline v1 {

////////////////////////////////////////////////////////////
//
//  Strip whitespace from the end of a line.
//
void
chomp(char *S) {
  char *t = S;

  while (*t != 0)
    t++;

  t--;

  while ((t >= S) && (isWhiteSpace(*t) == true))
    *t-- = 0;
}



////////////////////////////////////////////////////////////
//
//  In-place removes whitespace from either end of a string.
//
//  The non-whitespace portion of the string is shifted to the front of the
//  storage and unused space is filled with NUL bytes.  Any bytes after the
//  first NUL byte are not modified.
//
//  Using '@' to represent NUL bytes:
//
//    input:   s = '   hello world    @'
//    output:  s = 'hello world@@@@@@@@'
//
//    input:   s = '   hello@world@    '
//    output:  s = 'hello@@@@world@    '
//
//    input:   s = 'hello@kq4Bq4tab9q5B'
//    output:  s = 'hello@kq4Bq4tab9q5B'
//
//  The pointer to s is returned from the function.
//
char *
trimString(char *s, bool bgn, bool end) {
  int64 i=0, b=0, e=0;

  while (s[e] != 0)                          //  Find the end of the string.
    e++;                                     //  We later expect 'e' to be the
  e--;                                       //  last letter of the string.

  if (end)                                   //  Trim spaces at the end.
    while ((e >= 0) && (s[e] == ' '))        //
      s[e--] = 0;                            //

  if (bgn)                                   //  Trim spaces from the start.
    while ((b < e) && (s[b] == ' '))         //
      b++;                                   //

  while (b <= e)                             //  Shift the string to actually
    s[i++] = s[b++];                         //  remove the spaces at the start.

  while (i <= e)                             //  Terminate the new string and
    s[i++] = 0;                              //  erase left over crud.

  return s;
}





static uint64   dsLen[32] = {0};
static uint64   dsMax[32] = {0};
static char    *dsStr[32]    = { nullptr };
static uint32   pos        = 31;


char const *
displayString(char const *s) {
  pos = (pos + 1) % 32;

  return displayString(s, dsStr[pos], dsLen[pos], dsMax[pos]);
}


char const *
displayString(char const *s, char *&d, uint64 &dLen, uint64 &dMax) {

  resizeArray(d, 0, dMax, 1024, _raAct::doNothing);   //  Allocate default space.

  if (s == nullptr) {                                 //  If no string input,
    d[0] = '<';                                       //  state so.
    d[1] = 'n'; d[2] = 'u'; d[3] = 'l'; d[4] = 'l'; d[5] = 'p'; d[6] = 't'; d[7] = 'r';
    d[8] = '>'; d[9] =  0;
    return d;
  }

  dLen = 1;                                           //  Compute (a little bit more than)
  for (char const *c=s; *c; c++) {                    //  the actual length of string
    if      (*c < 0x20)  dLen += 4;                   //  we'll be returning.
    else if (*c < 0x7f)  dLen += 1;
    else                 dLen += 4;
  }
  resizeArray(d, 0, dMax, dLen, _raAct::doNothing);   //  Then make sure we have enough space.

  dLen = 0;                                           //  Now just process each letter into
  for (char const *c=s; *c; c++) {                    //  a nice symbol.
    if (*c < 0x20) {
      switch (*c) {
        case 0x00:  d[dLen++] = '<';  d[dLen++] = 'N';  d[dLen++] = 'U';  d[dLen++] = 'L';  d[dLen++] = '>';  break;
        case 0x01:  d[dLen++] = '<';  d[dLen++] = 'S';  d[dLen++] = 'O';  d[dLen++] = 'H';  d[dLen++] = '>';  break;
        case 0x02:  d[dLen++] = '<';  d[dLen++] = 'S';  d[dLen++] = 'T';  d[dLen++] = 'X';  d[dLen++] = '>';  break;
        case 0x03:  d[dLen++] = '<';  d[dLen++] = 'E';  d[dLen++] = 'T';  d[dLen++] = 'X';  d[dLen++] = '>';  break;
        case 0x04:  d[dLen++] = '<';  d[dLen++] = 'E';  d[dLen++] = 'O';  d[dLen++] = 'T';  d[dLen++] = '>';  break;
        case 0x05:  d[dLen++] = '<';  d[dLen++] = 'E';  d[dLen++] = 'N';  d[dLen++] = 'Q';  d[dLen++] = '>';  break;
        case 0x06:  d[dLen++] = '<';  d[dLen++] = 'A';  d[dLen++] = 'C';  d[dLen++] = 'K';  d[dLen++] = '>';  break;
        case 0x07:  d[dLen++] = '<';  d[dLen++] = 'B';  d[dLen++] = 'E';  d[dLen++] = 'L';  d[dLen++] = '>';  break;
        case 0x08:  d[dLen++] = '<';  d[dLen++] = 'B';  d[dLen++] = 'S';                    d[dLen++] = '>';  break;
        case 0x09:  d[dLen++] = '<';  d[dLen++] = 'H';  d[dLen++] = 'T';                    d[dLen++] = '>';  break;
        case 0x0a:  d[dLen++] = '<';  d[dLen++] = 'L';  d[dLen++] = 'F';                    d[dLen++] = '>';  break;
        case 0x0b:  d[dLen++] = '<';  d[dLen++] = 'V';  d[dLen++] = 'T';                    d[dLen++] = '>';  break;
        case 0x0c:  d[dLen++] = '<';  d[dLen++] = 'F';  d[dLen++] = 'F';                    d[dLen++] = '>';  break;
        case 0x0d:  d[dLen++] = '<';  d[dLen++] = 'C';  d[dLen++] = 'R';                    d[dLen++] = '>';  break;
        case 0x0e:  d[dLen++] = '<';  d[dLen++] = 'S';  d[dLen++] = 'O';                    d[dLen++] = '>';  break;
        case 0x0f:  d[dLen++] = '<';  d[dLen++] = 'S';  d[dLen++] = 'I';                    d[dLen++] = '>';  break;

        case 0x10:  d[dLen++] = '<';  d[dLen++] = 'D';  d[dLen++] = 'L';  d[dLen++] = 'E';  d[dLen++] = '>';  break;
        case 0x11:  d[dLen++] = '<';  d[dLen++] = 'D';  d[dLen++] = 'C';  d[dLen++] = '1';  d[dLen++] = '>';  break;
        case 0x12:  d[dLen++] = '<';  d[dLen++] = 'D';  d[dLen++] = 'C';  d[dLen++] = '2';  d[dLen++] = '>';  break;
        case 0x13:  d[dLen++] = '<';  d[dLen++] = 'D';  d[dLen++] = 'C';  d[dLen++] = '3';  d[dLen++] = '>';  break;
        case 0x14:  d[dLen++] = '<';  d[dLen++] = 'D';  d[dLen++] = 'C';  d[dLen++] = '4';  d[dLen++] = '>';  break;
        case 0x15:  d[dLen++] = '<';  d[dLen++] = 'N';  d[dLen++] = 'A';  d[dLen++] = 'K';  d[dLen++] = '>';  break;
        case 0x16:  d[dLen++] = '<';  d[dLen++] = 'S';  d[dLen++] = 'Y';  d[dLen++] = 'N';  d[dLen++] = '>';  break;
        case 0x17:  d[dLen++] = '<';  d[dLen++] = 'E';  d[dLen++] = 'T';  d[dLen++] = 'B';  d[dLen++] = '>';  break;
        case 0x18:  d[dLen++] = '<';  d[dLen++] = 'C';  d[dLen++] = 'A';  d[dLen++] = 'N';  d[dLen++] = '>';  break;
        case 0x19:  d[dLen++] = '<';  d[dLen++] = 'E';  d[dLen++] = 'M';                    d[dLen++] = '>';  break;
        case 0x1a:  d[dLen++] = '<';  d[dLen++] = 'S';  d[dLen++] = 'U';  d[dLen++] = 'B';  d[dLen++] = '>';  break;
        case 0x1b:  d[dLen++] = '<';  d[dLen++] = 'E';  d[dLen++] = 'S';  d[dLen++] = 'C';  d[dLen++] = '>';  break;
        case 0x1c:  d[dLen++] = '<';  d[dLen++] = 'F';  d[dLen++] = 'S';                    d[dLen++] = '>';  break;
        case 0x1d:  d[dLen++] = '<';  d[dLen++] = 'G';  d[dLen++] = 'S';                    d[dLen++] = '>';  break;
        case 0x1e:  d[dLen++] = '<';  d[dLen++] = 'R';  d[dLen++] = 'S';                    d[dLen++] = '>';  break;
        case 0x1f:  d[dLen++] = '<';  d[dLen++] = 'U';  d[dLen++] = 'S';                    d[dLen++] = '>';  break;
      }
    }
    else if (*c < 0x7f) {
      d[dLen++] = *c;
    }
    else if (*c == 0x7f) {
      d[dLen++] = '<';  d[dLen++] = 'D';  d[dLen++] = 'E';  d[dLen++] = 'L';  d[dLen++] = '>';  break;
    }
    else {
      d[dLen++] = '<';  d[dLen++] = 'x';  d[dLen++] = toHexDigit(*c >> 4);
                                          d[dLen++] = toHexDigit(*c);         d[dLen++] = '>';  break;
    }
  }

  d[dLen++] = 0;

  return d;
}


}  //  merylutil::strings::v1
