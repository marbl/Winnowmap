
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

#include "writing-v1.H"

//
//  Write a (potentially ginormous) array of data to disk.
//
//  The write is broken into a series of writes of at most 32 MB at once.
//  Comments on earlier versions indicated this was needed for OSF1 (V5.1).
//  In testing, certainly FreeBSD 11 still isn't happy writing 16 GB of data
//  at once; it seems to truncate to 32-bit somewhere.
//

namespace merylutil::inline files::inline v1 {

void
writeToFile(void const  *objects,
            char const  *description,
            uint64       objectSize,
            uint64       nObjects,
            FILE        *file) {
  uint64  nWritten  = 0;
  uint64  blockSize = (uint64)32 * 1024 * 1024 / objectSize;

  while (nWritten < nObjects) {
    uint64  toWrite = std::min(blockSize, nObjects - nWritten);

    errno = 0;
    uint64 written = fwrite(((char *)objects) + nWritten * objectSize, objectSize, toWrite, file);
    nWritten += written;

    if (written != toWrite)
      fprintf(stderr, "writeToFile()-- After writing %lu out of %lu '%s' objects (%lu bytes each): %s\n",
              nWritten, nObjects, description, objectSize, strerror(errno)), exit(1);
  }

  if (nWritten != nObjects)
    fprintf(stderr, "writeToFile()-- After writing %lu out of %lu '%s' objects (%lu bytes each): Short write\n",
            nWritten, nObjects, description, objectSize), exit(1);
}

}  //  namespace merylutil::files::v1
