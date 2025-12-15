
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

#ifndef MERYLUTIL_FILES_COMPRESSED_V1_H
#define MERYLUTIL_FILES_COMPRESSED_V1_H

#include "types.H"

namespace merylutil::inline files::inline v1 {

enum cftType {
  cftNONE  = 0,
  cftGZ    = 1,
  cftLZIP  = 2,
  cftBZ2   = 3,
  cftXZ    = 4,
  cftZSTD  = 5,
  cftSTDIN = 6
};

cftType  compressedFileType(char const *filename);

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_FILES_COMPRESSED_V1_H
