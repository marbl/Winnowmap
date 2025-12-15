
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

#ifndef MERYLUTIL_FILES_COMPRESSED_V1_WRITING_H
#define MERYLUTIL_FILES_COMPRESSED_V1_WRITING_H

#include "types.H"
#include "compressed-v1.H"

namespace merylutil::inline files::inline v1 {

class compressedFileWriter {
public:
  compressedFileWriter(char const *filename, uint32 cLevel=1, uint32 nThreads=2);
  ~compressedFileWriter();

  void  close(void);

  FILE *operator*(void)     {  return(_file);          };
  FILE *file(void)          {  return(_file);          };

  char *filename(void)      {  return(_filename);      };

  bool  isCompressed(void)  {  return(_pipe == true);  };

private:
  FILE  *_file     = nullptr;
  char  *_filename = nullptr;

  bool   _pipe     = false;
  bool   _stdi     = false;
};

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_FILES_COMPRESSED_V1_WRITING_H
