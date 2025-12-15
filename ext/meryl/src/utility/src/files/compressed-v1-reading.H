
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

#ifndef MERYLUTIL_FILES_COMPRESSED_V1_READING_H
#define MERYLUTIL_FILES_COMPRESSED_V1_READING_H

#include "types.H"
#include "compressed-v1.H"

namespace merylutil::inline files::inline v1 {

class compressedFileReader {
public:
  compressedFileReader(char const *filename, int32 threads=0, cftType type=cftType::cftNONE);
  ~compressedFileReader();

  void  reopen(int32 threads=0);
  void  close(void);

  FILE *operator*(void)     {  return _file;  }
  FILE *file(void)          {  return _file;  }

  char *filename(void)      {  return _filename;  }

  bool  isCompressed(void)  {  return (_pipe == true);                       }
  bool  isNormal(void)      {  return (_pipe == false) && (_stdi == false);  }

  //  Seekable if it is not stdin and not a pipe.
  //  Reopenable if it is not stdin.
  bool  isSeekable(void)    {  return (_stdi == false) && (_pipe == false);   }
  bool  isReopenable(void)  {  return (_stdi == false);                       }

  //  Return the file line-by-line.
  //    while (F->readLine())
  //      puts(F->line());
  //
  bool    readLine(void) {
    if (merylutil::readLine(_line, _lineLen, _lineMax, _file) == true) {
      _lineNum++;
      return(true);
    } else {
      return(false);
    }
  }
  char   *line   (void)   { return(_line);    }
  uint32  lineLen(void)   { return(_lineLen); }
  uint32  lineNum(void)   { return(_lineNum); }


private:
  FILE     *_file     = nullptr;
  char     *_filename = nullptr;

  int32     _nThreads = 2;
  cftType   _type     = cftType::cftNONE;

  bool      _pipe     = false;
  bool      _stdi     = false;

  uint32    _lineMax  = 0;
  uint32    _lineLen  = 0;
  uint64    _lineNum  = 0;
  char     *_line     = nullptr;     //  For readLine()
};

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_FILES_COMPRESSED_V1_READING_H
