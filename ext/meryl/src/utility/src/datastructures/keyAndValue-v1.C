
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
#include "arrays.H"

namespace merylutil::inline strings::inline v1 {

bool
KeyAndValue::find(const char *line) {
  char  *ptr = nullptr;

  //  Reset our state, but return fail if there is no line.

  _key = nullptr;
  _val = nullptr;

  if (isEmptyString(line) == true)
    return(false);

  //  Copy the string so we can do bad things to it.

  duplicateArray(_line, _lineLen, _lineMax, line, (uint32)strlen(line) + 1);

  //  Zip ahead until the first non-space letter.
  //
  //  If the letter is a comment or the delimiter, we're done; there is no key.

  ptr = _line;

  while (isWhiteSpace(*ptr) == true)          //  Spaces before the key.
    ptr++;

  if ((*ptr == 0) ||
      (isComment(*ptr) == true) ||
      (isDelimiter(*ptr) == true))
    return(false);

  _key = ptr;

  //  Keep zipping ahead until the end of the line.
  //    Detect the first comment mark that is preceeded by a space.
  //      Change it to NUL to terminate the string and return.
  //
  //    Detect the first key=value delimiter.
  //      Change it to a space so we can iterate over it.
  //      Change all consecutive delimiters to a space too.

  char *equals    = nullptr;
  char *eol       = nullptr;
  bool  lastspace = false;

  while (1) {
    eol = ptr;

    if ((lastspace == true) && (isComment(*ptr) == true)) {
      *ptr = 0;
      break;
    }

    lastspace = isWhiteSpace(*ptr);

    if ((isDelimiter(*ptr) == true) && (equals == nullptr)) {
      equals = ptr;    //  Remember the first letter in the delimiter

      while ((*ptr != 0) && (isDelimiter(*ptr) == true))
        *ptr++ = ' ';

      ptr--;           //  Back up to the last delimiter letter.

      eol    = ptr;    //  Update eol since we possibly moved ptr ahead.
    }

    if (*ptr == 0)
      break;

    ptr++;
  }

  //  If no delimiter, we're done.  There cannot be a key/value pair.

  if (equals == nullptr)
    return(false);

  //  Cleanup 1:  Find the last letter in the key make the key stop there.

  while (isWhiteSpace(*equals) == true)
    equals--;

  equals++;      //  Move from the last letter of the key.
  *equals = 0;   //  Terminate the key string.
  equals++;      //  Move to the next letter, either space or the value.

  //  Cleanup 2: Find the first letter of the value.
  //  If we're at eol now, return true with an empty value string.

  while (isWhiteSpace(*equals) == true)
    equals++;

  _val = equals;

  if (equals == eol)
    return(true);

  //  Cleanup 3: Find the last letter of the value and make the value stop
  //  there.

  assert(*eol == 0);

  eol--;

  while (isWhiteSpace(*eol) == true) {
    *eol = 0;
    eol--;
  }

  return(true);
}

}  //  merylutil::strings::v1
