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

#ifndef MERYLUTIL_DS_REGEX_V1_H
#define MERYLUTIL_DS_REGEX_V1_H

#include <regex.h>

#include "types.H"
#include "arrays.H"
#include "system.H"

namespace merylutil::inline regex::inline v1 {

class RegEx {
public:
  RegEx();
  ~RegEx();

  bool        compile(char const *pattern, std::vector<char const *> *errors = nullptr);
  bool        compile(char const *pattern, std::vector<char const *> &errors) {
    return compile(pattern, &errors);
  }

  //ol        compile_concat(char const *pattern, ...);

  bool        match(char const *line, std::vector<char const *> *errors = nullptr);
  bool        match(char const *line, std::vector<char const *> &errors) {
    return match(line, &errors);
  }

  uint32      matchesLen(void)      { return _rmLen; }

  uint32      getMatchBgn(uint32 m) { return _rm[m].rm_so; }
  uint32      getMatchEnd(uint32 m) { return _rm[m].rm_eo; }
  uint32      getMatchLen(uint32 m) { return _rm[m].rm_eo - _rm[m].rm_so; }

  //  getMatch(0) is the whole string that matched.
  //  getMatch(1) and on are the individual subexpressions that matched.
  char const *getMatch(uint32 m)    { return _mm[m]; }

private:
  bool        reportError(int r, char const *fmt, char const *arg, std::vector<char const *> *errors);

private:
  regex_t     *_rx       = nullptr;   //  The compiled regex.

  uint32       _rmLen    = 0;         //  Output from executing the regex on
  regmatch_t  *_rm       = nullptr;   //  some input string.

  uint32       _mmMax    = 0;         //  _mmString is a copy of ALL subexpressions
  char        *_mmString = nullptr;   //  that matched, NUL terminated; _mm has pointers
  char       **_mm       = nullptr;   // into that string for consumption by the client.

  uint32       _esMax    = 0;         //  An error string, simplifies error reorting
  char        *_es       = nullptr;   //  in compile() and match().
};

}  //  merylutil::regex::v1

#endif  //  MERYLUTIL_DS_REGEX_V1_H
