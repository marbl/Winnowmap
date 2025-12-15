
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


namespace merylutil::inline regex::inline v1 {

RegEx::RegEx() {
}

RegEx::~RegEx() {
  if (_rx)
    regfree(_rx);

  delete    _rx;
  delete [] _rm;

  delete [] _mmString;
  delete [] _mm;

  delete [] _es;
}


bool
RegEx::compile(char const *pattern, std::vector<char const *> *errors) {

  if (_rx)                                     //  If an existing regex_t,
    regfree(_rx);                              //  release it so we can compile
  else                                         //  a new expression, otherwise,
    _rx = new regex_t;                         //  allocate an empty one.

  int r = regcomp(_rx, pattern, REG_EXTENDED | REG_ICASE);

  if (r)
    return reportError(r, "pattern '%s' failed to compile: %s", pattern, errors);

  delete [] _rm;
  delete [] _mm;

  _rmLen = _rx->re_nsub + 1;                   //  We could recycle these,
  _rm    = new regmatch_t [_rmLen];            //  but it's a little simpler
  _mm    = new char *     [_rmLen];            //  to just get new ones.

  return true;
}


#if 0
//  Untested.  Was supposed to allow compilation of patters composed from
//  various strings:
//     char const *hexdigit = "[0-9A-Fa-f]";
//     compile_concat("number=0x", hexdigit, hexdigit, "...")
//
bool
RegEx::compile_concat(char const *pattern, ...) {
  char    *pat, *p;
  uint32   patLen = strlen(pattern) + 1;       //  For the NUL byte.
  va_list  ap;

  va_start(ap, pattern);
  for (char const *a = va_arg(ap, char const *); a; a = va_arg(ap, char const *))
    patLen += strlen(a);
  va_end(ap);

  pat = p = new char [patLen];

  for (char const *a=pattern; *a; )
    *p++ = *a++;
  *p = 0;

  va_start(ap, pattern);
  for (char const *a = va_arg(ap, char const *); a; a = va_arg(ap, char const *)) {
    while (*a)
      *p++ = *a++;
    *p = 0;
  }
  va_end(ap);

  return compile(pat);
}
#endif


bool
RegEx::match(char const *line, std::vector<char const *> *errors) {

  if (_rx == nullptr)
    return false;

  for (uint32 ii=0; ii<_rmLen; ii++)           //  Initialize all submatches
    _rm[ii].rm_so = _rm[ii].rm_eo = -1;        //  to empty intervals.

  int r = regexec(_rx, line, _rmLen, _rm, 0);  //

  if (r == REG_NOMATCH)                        //  If no match, return
    return false;                              //  failure.

  if (r)
    return reportError(r, "regex failed on string '%s': %s", line, errors);

  uint32 ll = 0;                               //  Compute total length of all
  for (uint32 ii=0; ii<_rmLen; ii++)           //  submatches and then allocate
    ll += _rm[ii].rm_eo - _rm[ii].rm_so + 1;   //  a string to store them all.

  resizeArray(_mmString, 0, _mmMax, ll, _raAct::doNothing);

  for (uint32 pp=0, ii=0, xx=0; ii<_rmLen; ii++) {
    _mm[ii] = _mmString + pp;                  //  Set pointer to start of submatch string.

    xx = _rm[ii].rm_so;                        //  Copy submatch string
    while (xx<_rm[ii].rm_eo)                   //  to storage.  Submatches with no match
      _mmString[pp++] = line[xx++];            //  are just the NUL byte.

    _mmString[pp++] = 0;                       //  Terminate submatch string in storage.
  }

  return true;
}


bool
RegEx::reportError(int r, char const *fmt, char const *arg, std::vector<char const *> *errors) {
  resizeArray(_es, 0, _esMax,                //  Allocate the error string to
              regerror(r, _rx, nullptr, 0),  //  the size reported by regerror().
              _raAct::doNothing);            //  (it's freed in the desctructor)

  regerror(r, _rx, _es, _esMax);             //  Then copy the error to our string.

  return fatalError(true, errors, fmt, arg, _es);
}


}  //  merylutil::regex::v1
