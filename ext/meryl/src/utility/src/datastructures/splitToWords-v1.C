
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

void
splitToWords::split(const char *line, splitType type) {

  //  Initialize to no words and no characters.
  //  Then return if the input line is empty.

  _wordsLen = 0;
  _charsLen = 0;

  if (isEmptyString(line) == true)
    return;

  //  Initialize the separator array based on splitType, if needed.

  if ((type == splitWords) ||
      (type == splitWhitespace)) {
    clearsc();
    setsc(' ');
    setsc('\t');
    setsc('\n');
    setsc('\r');
  }

  if (type == splitPaths) {
    clearsc();
    setsc('/');
  }

  if (type == splitTabs) {
    clearsc();
    setsc('\t');
  }

  if (type == splitLines) {
    clearsc();
    setsc('\n');
    setsc('\r');
  }

  //  Count the number of words and chars in the input line, then make
  //  sure there is space for us to store them.

  while (line[_charsLen] != 0)
    if (issc(line[_charsLen++]))
      _wordsLen++;

  resizeArray(_words, 0, _wordsMax, _wordsLen + 1);
  resizeArray(_chars, 0, _charsMax, _charsLen + 1);

  //  Clear all the words pointers, and copy the input line to our storage.
  //  This greatly simplifies the loop, as we don't need to worry about
  //  terminating the final word.

  memset(_words, 0,    sizeof(char *) * (_wordsLen + 1));
  memcpy(_chars, line, sizeof(char)   * (_charsLen + 1));

  //  Scan the line copy, converting word separators to NUL bytes.
  //  counting and saving the start of each word in _words.

  _wordsLen = 0;

  for (uint32 st=1, ii=0; ii < _charsLen; ii++) {
    if (issc(line[ii])) {                     //  If the character is a word
      _chars[ii] = 0;                         //  separator, convert to NUL,
      st         = true;                      //  and flag the next character
    }                                         //  as the start of a new word.

    else if (st) {                            //  Otherwise, if this is the
      _words[_wordsLen++] = _chars + ii;      //  start of a word, make
      st                  = false;            //  a new word.
    }
  }
}


void
splitToWords::clear(void) {
  _wordsLen = 0;
  _charsLen = 0;
}


void
splitToWords::erase(void) {

  delete [] _words;
  delete [] _chars;

  _wordsLen = 0;
  _wordsMax = 0;
  _words    = nullptr;

  _charsLen = 0;
  _charsMax = 0;
  _chars    = nullptr;
}


}  //  merylutil::strings::v1
