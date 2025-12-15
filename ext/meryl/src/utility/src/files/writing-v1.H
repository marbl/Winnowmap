
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

#ifndef MERYLUTIL_FILES_WRITING_V1_H
#define MERYLUTIL_FILES_WRITING_V1_H

#include "types.H"
#include "accessing-v1.H"



////////////////////////////////////////
//
//  Write objects to a disk file.  Two forms are provided:
//
//    writeToFile(POINTER, description, NUMBER_OF_OBJECTS, file);
//    writeToFile(OBJECT,  description,                    file);
//
//  The first will write some number of objects from an array, and expects a
//  pointer.  The second will write one object from an instance, and expects
//  an actual object.
//
//  The second form has two implementations, one for OBJECTS and one for
//  POINTERS.  The OBJECT variant writes the object directly, while the
//  POINTER varient writes from the pointed to object.  If you try to write a
//  pointer itself, you'll get unexpected results (but you shouldn't be
//  writing pointers to files anyway!).
//
//  For reading, see loadFromFile() in files-reading.H.
//

namespace merylutil::inline files::inline v1 {

void
writeToFile(void const  *objects,
            char const  *description,
            uint64       objectSize,
            uint64       nObjects,
            FILE        *file);

template<typename OBJ>
void
writeToFile(OBJ const *array, char const *description, uint64 nObjects, FILE *file) {
  writeToFile(array, description, sizeof(OBJ), nObjects, file);
}

template<typename OBJ>
void
writeToFile(OBJ &object, char const *description, FILE *file) {
  writeToFile(&object, description, sizeof(OBJ), 1, file);
}

template<typename OBJ>
void
writeToFile(OBJ const *object, char const *description, FILE *file) {
  writeToFile(object, description, sizeof(OBJ), 1, file);
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Write an array of objects to file '{prefix}{separator}{suffix}'.
//  If {suffix} is nullptr, then both {separator} and {suffix} are ignored.
//

namespace merylutil::inline files::inline v1 {

template<typename OBJ>
void    saveFile(char const *prefix, char separator, char const *suffix, OBJ *objects, uint64 numberToWrite) {
  FILE    *file = openOutputFile(prefix, separator, suffix);

  writeToFile(objects, "saveFile", numberToWrite, file);

  closeFile(file, prefix, separator, suffix);
}

template<typename OBJ>
void    saveFile(char const *prefix, OBJ *objects, uint64 numberToWrite) {
  saveFile(prefix, '.', nullptr, objects, numberToWrite);
}

}  //  namespace merylutil::files::v1


#endif  //  MERYLUTIL_FILES_WRITING_V1_H
