
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

#ifndef MERYLUTIL_DS_STRINGLIST_V1_H
#define MERYLUTIL_DS_STRINGLIST_V1_H

#include <vector>

#include "types.H"

//  stringList - efficiently store a collection of strings.
//
//  The intended purpose is to load a list of words read, one per line, from
//  a file.  Extended to support breaking at tabs or generic white space.
//  However, empty records are ignored - in tab mode, 'a b\t\tc\td' will
//  return strings 'a b', 'c', 'd'; the two tabs between the first two words
//  is treated as one separator.

namespace merylutil::inline strings::inline v1 {

class stringList {
public:
  stringList()                                              {                      };
  stringList(char const *filename, splitType st=splitWords) {  load(filename, st); };

  ~stringList() {
    for (uint64 ii=0; ii<_data.size(); ii++)
      delete [] _data[ii];
  }

  void            load(char const *filename, splitType st=splitLines);
  void            add(char const *line) {  _pointers.push_back(line);  }

  uint64          size(void)            {  return(_pointers.size());   }

  char const *    operator[](uint64 ii) {
    if (ii < _pointers.size())
      return(_pointers[ii]);
    else
      return(nullptr);
  }

  std::vector<char const *>  &getVector(void) { return(_pointers); }

private:
  std::vector<char const *>  _pointers;
  std::vector<char const *>  _data;
};

}  //  merylutil::strings::v1

#endif  //  MERYLUTIL_DS_STRINGLIST_V1_H
