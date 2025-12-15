
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

#ifndef MERYLUTIL_FILES_READING_V1_H
#define MERYLUTIL_FILES_READING_V1_H

#include "types.H"
#include "accessing-v1.H"



////////////////////////////////////////
//
//  Read a single line into the supplied array.
//   - The array can be unallocated (L == nullptr).
//   - Returns false if input file F is nullptr.
//   - Returns false if the input file F is EOF when the function is called.
//   - Otherwise returns true.
//
//  Trims traling whitespace from the line, including any new-lines and/or
//  carriage-returns.  Because of this, it is possible for the string returned to
//  be completely empty.  Do not use L[0] == 0 to test for EOF!
//
//  The v0 implementation is _probably_ slower; it is certainly not recently tested.
//
//  Implementations are in files/readLine-v0.C and files/readLine-v1.C.
//

namespace merylutil::inline files::v0 {
bool
readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F);
}

namespace merylutil::inline files::inline v1 {
bool
readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F);
}



////////////////////////////////////////
//
//  Read objects from a disk file.  Two forms are provided:
//
//    loadFromFile(POINTER, description, NUMBER_OF_OBJECTS, file);
//    loadFileFile(OBJECT,  description,                    file);
//
//  The first will load some number of objects into an array, and expects a
//  pointer.  The second will load one object into an instance, and expects
//  an actual object.
//
//  The second form has two implementations, one for OBJECTS and one for
//  POINTERS.  The OBJECT variant loads into the object directly, while the
//  POINTER varient loads into the pointed to object.
//
//  If 'exact' is true, the load will fail unless exactly the correct number
//  of objects is loaded.
//
//  For writing, see writeToFile() in files-writing.H.
//

namespace merylutil::inline files::inline v1 {

uint64
loadFromFile(void        *objects,
             char const  *description,
             uint64       objectSize,
             uint64       nObjects,
             FILE        *file,
             bool         exact=true);

template<typename OBJ>
uint64
loadFromFile(OBJ *array, char const *description, uint64 nObjects, FILE *file, bool exact=true) {
  return(loadFromFile(array, description, sizeof(OBJ), nObjects, file, exact));
}

template<typename OBJ>
uint64
loadFromFile(OBJ &object, char const *description, FILE *file, bool exact=true) {
  return(loadFromFile(&object, description, sizeof(OBJ), 1, file, exact));
}

template<typename OBJ>
uint64
loadFromFile(OBJ *object, char const *description, FILE *file, bool exact=true) {
  return(loadFromFile(object, description, sizeof(OBJ), 1, file, exact));
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Read an array of objects from file '{prefix}{separator}{suffix}'.
//  If {suffix} is nullptr, then both {separator} and {suffix} are ignored.
//

namespace merylutil::inline files::inline v1 {

template<typename OBJ>
void
loadFile(char const *prefix, char separator, char const *suffix, OBJ *objects, uint64 numberToLoad) {
  FILE    *file   = openInputFile(prefix, separator, suffix);
  uint64   length = sizeOfFile(file);

  if (numberToLoad == 0)
    numberToLoad = length / sizeof(OBJ);

  if (length / sizeof(OBJ) < numberToLoad) {
    if (suffix)
      fprintf(stderr, "loadFile()-- File '%s%c%s' contains " F_U64 " objects, but asked to load " F_U64 ".\n",
              prefix, separator, suffix, length / sizeof(OBJ), numberToLoad);
    else
      fprintf(stderr, "loadFile()-- File '%s' contains " F_U64 " objects, but asked to load " F_U64 ".\n",
              prefix, length / sizeof(OBJ), numberToLoad);
    exit(1);
  }

  loadFromFile(objects, "loadFile", numberToLoad, file);

  closeFile(file, prefix, separator, suffix);
}

template<typename OBJ>
void
loadFile(char const *prefix, OBJ *objects, uint64 numberToLoad) {
  loadFile(prefix, '.', nullptr, objects, numberToLoad);
}

}  //  namespace merylutil::files::v1

#endif  //  MERYLUTIL_FILES_READING_V1_H
