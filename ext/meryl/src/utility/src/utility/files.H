
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

#ifndef FILES_H
#define FILES_H

#include "types.H"
#include <vector>


//  Provides a safe and reliable mechanism for reading / writing
//  binary data.
//
//  Split writes/reads into smaller pieces, check the result of each
//  piece.  Really needed by OSF1 (V5.1), useful on other platforms to
//  be a little more friendly (big writes are usually not
//  interruptable).

#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

void    AS_UTL_findBaseFileName(char *basename, char const *filename);

bool    AS_UTL_readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F);

void    AS_UTL_mkdir(char const *dirname);
void    AS_UTL_rmdir(char const *dirname);

void    AS_UTL_symlink(char const *pathToFile, char const *pathToLink);

void    AS_UTL_unlink(char const *prefix, char separator='.', char const *suffix=NULL);

void    AS_UTL_rename(char const *oldname, char const *newname);
void    AS_UTL_rename(char const *oldprefix, char oldseparator, char const *oldsuffix,
                      char const *newprefix, char newseparator, char const *newsuffix);

bool    AS_UTL_makeReadOnly(char const *prefix, char separator='.', char const *suffix=NULL);
bool    AS_UTL_makeWritable(char const *prefix, char separator='.', char const *suffix=NULL);

//  Tests on existence of paths, files and directories.
bool    pathExists     (char const *prefix, char separator='.', char const *suffix=NULL);                        //  true if file or dir 'path' exists
bool    fileExists     (char const *prefix, char separator='.', char const *suffix=NULL, bool writable=false);   //  true if file        'path' exists (and is writable)
bool    directoryExists(char const *prefix, char separator='.', char const *suffix=NULL);                        //  true if         dir 'path' exists

off_t   AS_UTL_sizeOfFile(char const *path);
off_t   AS_UTL_sizeOfFile(FILE *file);

uint64  AS_UTL_timeOfFile(char const *path);
uint64  AS_UTL_timeOfFile(FILE *file);

off_t   AS_UTL_ftell(FILE *stream);
void    AS_UTL_fseek(FILE *stream, off_t offset, int whence);

//  Search CANU_INSTALL_PATH, MERYL_INSTALL_PATH and PATH for a data file.
char const *findSharedFile(char const *relpath, char const *filename);

//  Read a file-of-files into a vector
void    AS_UTL_loadFileList(char const *fileName, std::vector<char const *> &FILE);

FILE   *AS_UTL_openInputFile (char const *prefix, char separator='.', char const *suffix=NULL, bool doOpen=true);
FILE   *AS_UTL_openOutputFile(char const *prefix, char separator='.', char const *suffix=NULL, bool doOpen=true);

void    AS_UTL_closeFile(FILE *&F, char const *filename, char separator, char const *suffix, bool critical=true);
void    AS_UTL_closeFile(FILE *&F, char const *filename=NULL, bool critical=true);

void    AS_UTL_createEmptyFile(char const *prefix, char separator='.', char const *suffix=NULL);

template<typename OBJ>
void    AS_UTL_loadFile(char const *prefix, char separator, char const *suffix, OBJ *objects, uint64 numberToLoad) {
  FILE    *file   = AS_UTL_openInputFile(prefix, separator, suffix);
  uint64   length = AS_UTL_sizeOfFile(file);

  if (numberToLoad == 0)
    numberToLoad = length / sizeof(OBJ);

  if (length / sizeof(OBJ) < numberToLoad) {
    if (suffix)
      fprintf(stderr, "AS_UTL_loadFile()-- File '%s%c%s' contains " F_U64 " objects, but asked to load " F_U64 ".\n",
              prefix, separator, suffix, length / sizeof(OBJ), numberToLoad);
    else
      fprintf(stderr, "AS_UTL_loadFile()-- File '%s' contains " F_U64 " objects, but asked to load " F_U64 ".\n",
              prefix, length / sizeof(OBJ), numberToLoad);
    exit(1);
  }

  loadFromFile(objects, "loadFile", numberToLoad, file);

  AS_UTL_closeFile(file, prefix, separator, suffix);
}

template<typename OBJ>
void    AS_UTL_loadFile(char const *prefix, OBJ *objects, uint64 numberToLoad) {
  AS_UTL_loadFile(prefix, '.', NULL, objects, numberToLoad);
}

template<typename OBJ>
void    AS_UTL_saveFile(char const *prefix, char separator, char const *suffix, OBJ *objects, uint64 numberToWrite) {
  FILE    *file = AS_UTL_openOutputFile(prefix, separator, suffix);

  writeToFile(objects, "saveFile", numberToWrite, file);

  AS_UTL_closeFile(file, prefix, separator, suffix);
}

template<typename OBJ>
void    AS_UTL_saveFile(char const *prefix, OBJ *objects, uint64 numberToWrite) {
  AS_UTL_saveFile(prefix, '.', NULL, objects, numberToWrite);
}







//  Read or write an object from or to a disk file.  Two forms are provided:
//
//    loadFromFile(POINTER, description, NUMBER_OF_OBJECTS, file);
//    loadFileFile(OBJECT,  description,                    file);
//
//  The first will load some number of objects into an array, and expects a
//  pointer.  The second will load one object into an instance, and expects
//  an actual object.
//
//  The second form hsa two implementations, one for OBJECTS and one for
//  POINTERS.  The OBJECT variant loads into the object directly, while the
//  POINTER varient loads into the pointed to object.  If you try to write a
//  pointer itself, you'll get unexpected results (but you shouldn't be
//  writing pointers to files anyway!).
//
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



#include "files-fasta-fastq.H"
#include "files-compressed.H"
#include "files-buffered.H"
#include "files-buffered-implementation.H"
#include "files-memoryMapped.H"


#endif  //  FILES_H
