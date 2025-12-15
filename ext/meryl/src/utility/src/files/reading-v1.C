
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

#include "reading-v1.H"

//
//  Read a (potentially ginormous) array of data from disk.
//
//  The read is broken into a series of 32 MB chunks.  See comments in
//  writing-v1.C for details.
//

namespace merylutil::inline files::inline v1 {

uint64
loadFromFile(void        *objects,
             char const  *description,
             uint64       objectSize,
             uint64       nObjects,
             FILE        *file,
             bool         exact) {
  uint64  nLoaded   = 0;
  uint64  blockSize = (uint64)32 * 1024 * 1024 / objectSize;

  if (file == nullptr)
    return(nLoaded);

  while (nLoaded < nObjects) {
    uint64  toLoad = std::min(blockSize, nObjects - nLoaded);

    errno = 0;
    uint64 loaded = fread(((char *)objects) + nLoaded * objectSize, objectSize, toLoad, file);
    nLoaded += loaded;

    //  If we've loaded all requested objects, stop loading.
    if (nLoaded == nObjects)
      break;

    //  If end-of-file and expecting exact, fail.
    if (feof(file) && (exact == true))
      fprintf(stderr, "loadFromFile()-- After loading %lu out of %lu '%s' objects (%lu bytes each): End of file\n",
              nLoaded, nObjects, description, objectSize), exit(1);

    //  If end-of-file, stop loading.
    if (feof(file))
      break;

    //  If an error (that is not resumable), fail.
    if (ferror(file) && (errno != EINTR))
      fprintf(stderr, "loadFromFile()-- After loading %lu out of %lu '%s' objects (%lu bytes each): %s\n",
              nLoaded, nObjects, description, objectSize, strerror(errno)), exit(1);
  }

  //  Successful load.  Return how many objects were loaded.
  return(nLoaded);
}

}  //  namespace merylutil::files::v1
