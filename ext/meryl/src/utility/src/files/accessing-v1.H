
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

#ifndef MERYLUTIL_FILES_ACCESSING_V1_H
#define MERYLUTIL_FILES_ACCESSING_V1_H

#include "types.H"


////////////////////////////////////////
//
//  A collection of wrappers for creating, renaming, finding et cetera files
//  and directories.
//
//  See reading-v1.H and writing-v1.H for operations that read and write data.
//

namespace merylutil::inline files::inline v1 {

//  OPERATIONS ON FILE NAMES

char const *
constructPathName(char const *prefix, char separator, char const *suffix);

void    findBaseFileName(char *basename, char const *filename);
void    stripLastSuffix(char *basename, char const *filename);

bool    rename(char const *oldname, char const *newname, bool fatal=true);
bool    rename(char const *oldprefix, char oldseparator, char const *oldsuffix,
               char const *newprefix, char newseparator, char const *newsuffix, bool fatal=true);

bool    symlink(char const *pathToFile, char const *pathToLink, bool fatal=true);

bool    pathExists     (char const *prefix, char separator='.', char const *suffix=nullptr);                        //  true if file or dir 'path' exists
bool    fileExists     (char const *prefix, char separator='.', char const *suffix=nullptr, bool writable=false);   //  true if file        'path' exists (and is writable)
bool    directoryExists(char const *prefix, char separator='.', char const *suffix=nullptr);                        //  true if         dir 'path' exists

//  OPERATIONS ON DIRECTORIES

bool    mkdir (char const *dirname, bool fatal=true);
bool    rmdir (char const *dirname, bool fatal=true);

bool    mkpath(char const *dirname, bool fatal=true);

//  OPERATIONS ON FILES

bool    unlink(char const *path, bool fatal=false);
bool    unlink(char const *prefix, char separator, char const *suffix, bool fatal=false);

//  OPERATIONS ON PERMISSIONS

bool    makeReadOnly(char const *prefix, char separator='.', char const *suffix=nullptr);
bool    makeWritable(char const *prefix, char separator='.', char const *suffix=nullptr);

//  OPERATIONS ON METADATA

off_t   sizeOfFile(char const *prefix, char separator='.', char const *suffix=nullptr);
off_t   sizeOfFile(FILE *file);

uint64  timeOfFile(char const *prefix, char separator='.', char const *suffix=nullptr);
uint64  timeOfFile(FILE *file);

bool    setFileTime(uint64 time, char const *prefix, char separator='.', char const *suffix=nullptr);
bool    setFileTime(double time, char const *prefix, char separator='.', char const *suffix=nullptr);
bool    setFileTime(struct timespec time, char const *prefix, char separator='.', char const *suffix=nullptr);

off_t   ftell(FILE *stream);
void    fseek(FILE *stream, off_t offset, int whence);

//  CREATING or USING FILES

FILE   *openInputFile (char const *prefix, char separator='.', char const *suffix=nullptr, bool doOpen=true);
FILE   *openOutputFile(char const *prefix, char separator='.', char const *suffix=nullptr, bool doOpen=true);

void    closeFile(FILE *&F, char const *prefix, char separator, char const *suffix, bool critical=true);
void    closeFile(FILE *&F, char const *pathname=nullptr, bool critical=true);

void    createEmptyFile(char const *prefix, char separator='.', char const *suffix=nullptr);

}  //  namespace merylutil::files::v1

#endif  //  MERYLUTIL_FILES_ACCESSING_V1_H
