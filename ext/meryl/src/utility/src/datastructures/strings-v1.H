
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

#ifndef MERYLUTIL_DS_STRINGS_V1_H
#define MERYLUTIL_DS_STRINGS_V1_H

#include "types.H"

#include <set>
#include <vector>

namespace merylutil::inline strings::inline v1 {

//  Remove whitespace from the end of a line.
void   chomp(char *S);
char  *trimString(char *s, bool bgn=true, bool end=true);

//  Return true if string s is unallocated or empty.
inline
bool
isEmptyString(char const *s) {
  return((s == nullptr) || (s[0] == 0));
}

//  Display a string, replacing binary control signals (LF, CR, FF, and all
//  the non-printable stuff) with codes '<LF>', '<CR>', '<FF'>, '<xHH>'.
//
//  The first form uses an internal collection of 32 buffers.
//  The second form returns a pointer to d.
//
char const *displayString(char const *s);
char const *displayString(char const *s, char *&d, uint64 &dLen, uint64 &dMax);

inline
char const *displayLetter(char c) {
  char s[2] = {c, 0};
  return displayString(s);
}

}  //  merylutil::strings::v1

#endif  //  MERYLUTIL_DS_STRINGS_V1_H
