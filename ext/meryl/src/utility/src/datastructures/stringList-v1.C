
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
#include "files.H"

namespace merylutil::inline strings::inline v1 {

static
bool
loadMore(char *&chunk, uint32 &cMax, uint32 &cLen, uint32 &cPos,
         compressedFileReader *file) {

  char *nhunk = new char [cMax];                //  Allocate another chunk of data.

  for (uint32 ii=0; ii<cLen - cPos; ii++)       //  Move unused bytes from the end
    nhunk[ii] = chunk[cPos + ii];               //  of the chunk to the new chunk.

  cLen -= cPos;                                 //  Reset cLen and cPos to show
  cPos -= cPos;                                 //  the moved data.

  chunk = nhunk;                                //  Save the new chunk.

  if (cMax - cLen < cMax / 2)                   //  Grab more space if the moved
    resizeArray(chunk, cLen, cMax, 2 * cMax);   //  data is big.

  cLen += loadFromFile(chunk + cLen,            //  Load more data.
                       "stringListData",
                       1, cMax - cLen, file->file(), false);

  if (cLen == 0) {                              //  Return true if there is
    delete [] chunk;                            //  data in the chunk,
    chunk = nullptr;                            //  or delete the chunk we
    return(false);                              //  just allocated.
  } else {
    return(true);
  }
}

void
stringList::load(char const *filename, splitType st) {
  uint32                 cMax  = 16;
  uint32                 cLen  = 0;
  uint32                 cPos  = 0;
  compressedFileReader  *file  = new compressedFileReader(filename);
  char                  *chunk = nullptr;

  bool   splitWS  = (st == splitWords || st == splitWhitespace);
  bool   splitTAB = (st == splitTabs);
  bool   splitEOL = (st == splitLines);

  while (loadMore(chunk, cMax, cLen, cPos, file)) {

    _data.push_back(chunk);

    //  Skip whitespace at the start of the chunk.  Usually only occurs on
    //  the first chunk load, but can occur if we're lucky enough to find
    //  the end-of-string before the last chunk ended.

    while ((cPos < cLen) &&
           (isWhiteSpace(chunk[cPos]) == true))
      cPos++;

    //  Search ahead until the next separator - either generic whitespace
    //  or newline/linefeed.
    //    cPos is the start of the next string to load
    //    sEnd is the end of the string
    //
    for (uint32 sEnd=cPos; (sEnd < cLen); sEnd++) {
      if ((splitWS  && isWhiteSpace(chunk[sEnd])) ||
          (splitTAB && isTab       (chunk[sEnd])) ||
          (splitEOL && isEndOfLine (chunk[sEnd]))) {

        chunk[sEnd] = 0;                              //  Terminate string and
        _pointers.push_back(chunk + cPos);            //  save the pointer to it.

        cPos = sEnd + 1;                              //  Advance to the next letter.

        while ((cPos < cLen) &&                       //  Skip whitespace to find
               (isWhiteSpace(chunk[cPos]) == true))   //  the next string start.
          cPos++;

        //  Update loop iterator to the last letter we processed, cPos, and
        //  let the next iteration process the letter after the first.
        sEnd = cPos;
      }
    }

    //  We finished processing all strings in the block.  Leave cPos as is,
    //  instead of updating to sEnd, so that we'll save (and process again)
    //  all bytes from the start of the current string.  The next iteration
    //  of the while() will loadMore() for us.
  }

  delete file;
}

}  //  merylutil::strings::v1
