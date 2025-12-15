
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

#include "arrays.H"
#include "files.H"

namespace merylutil::inline files::inline v1 {


cftType
compressedFileType(char const *filename) {

  if ((filename == NULL) || (filename[0] == 0) || (strcmp(filename, "-") == 0))
    return(cftSTDIN);

  int32  len = strlen(filename);

  if      ((len > 3) && (strcasecmp(filename + len - 3, ".gz") == 0))
    return(cftGZ);

  else if ((len > 3) && (strcasecmp(filename + len - 3, ".lz") == 0))
    return(cftLZIP);

  else if ((len > 4) && (strcasecmp(filename + len - 4, ".bz2") == 0))
    return(cftBZ2);

  else if ((len > 3) && (strcasecmp(filename + len - 3, ".xz") == 0))
    return(cftXZ);

  else if ((len > 4) && (strcasecmp(filename + len - 4, ".zst") == 0))
    return(cftZSTD);

  else
    return(cftNONE);
}


}  //  merylutil::files::v1


