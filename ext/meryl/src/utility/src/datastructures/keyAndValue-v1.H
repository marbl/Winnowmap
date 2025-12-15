
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

#ifndef MERYLUTIL_DS_KEYANDVALUE_V1_H
#define MERYLUTIL_DS_KEYANDVALUE_V1_H

#include "types.H"

namespace merylutil::inline strings::inline v1 {

//  Convert a line into a key-value pair.
//
//  The line should be of the form:
//   - 'key'           find() returns true with value == nullptr
//   - 'key = value'   find() returns true
//   - 'key : value'   find() returns true
//
//  In all cases the line is modified by removing whitespace and inserting
//  NUL characters at the end of the key and value.
//
//  Comments are accepted with either '!' or '#' at the start
//  of the line, or preceeded by a white-space letter.
//
class KeyAndValue {
public:
  KeyAndValue(const char *line = nullptr)   {  find(line);       };
  ~KeyAndValue()                            {  delete [] _line;  };

public:
  bool     find(const char *line);

public:
  char    *key(void)           { return(_key); };
  char    *value(void)         { return(_val); };

private:
  uint32  _lineMax = 0;
  uint32  _lineLen = 0;
  char   *_line    = nullptr;

  char   *_key     = nullptr;
  char   *_val     = nullptr;
};

}  //  merylutil::strings::v1

#endif  //  MERYLUTIL_DS_KEYANDVALUE_V1_H
