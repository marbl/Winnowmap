
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

#include "accessing-v1.H"
#include "system.H"

#include <fcntl.h>
#include <sys/stat.h>

namespace merylutil::inline files::inline v1 {

char const *
constructPathName(char const *prefix, char separator, char const *suffix) {
  uint32  plen = (prefix == nullptr) ? 0 : strlen(prefix);
  uint32  dlen = (suffix == nullptr) ? 0 : 1;
  uint32  slen = (suffix == nullptr) ? 0 : strlen(suffix);
  uint32  tlen = plen + dlen + slen + 1;

  char   *path = new char [tlen];
  uint32  p    = 0;

  while ((prefix) && (*prefix))   //  Copy the prefix to path, if it exists.
    path[p++] = *prefix++;

  if (suffix)                     //  Append the separator, if there is a suffix.
    path[p++] = separator;

  while ((suffix) && (*suffix))   //  Append the suffix, if it exists.
    path[p++] = *suffix++;

  path[p] = 0;                    //  Terminate the string.

  assert(p <= tlen);

  return path;
}

}  //  namespace merylutil::files::v1


////////////////////////////////////////
//
//  Return the basename of a path -- that is, strip off any and all file
//  extensions: anything after the first dot after the last slash is removed.
//    d.1/name.ext.Z -> d.1/d.2/name
//
//  If 'filename' refers to a directory, a copy of 'filename' is returned in
//  'basename' - no extensions are stripped.
//
//  Assumes 'basename' has at least as much storage as 'filename'.
//
namespace merylutil::inline files::inline v1 {

void
findBaseFileName(char *basename, char const *filename) {

  strcpy(basename, filename);               //  Copy in-name to out-name.

  if (directoryExists(basename))            //  Stop if this is an actual
    return;                                 //  directory.

  char  *slash = strrchr(basename, '/');    //  Search backwards for
  char  *dot   = nullptr;                   //  the last slash.

  if (slash)                                //  Search forward from the
    dot = strchr(slash, '.');               //  last slash, or from the
  else                                      //  start of the in-name, for
    dot = strchr(basename, '.');            //  the first dot.

  if (dot)                                  //  If a dot was found, terminate
    *dot = 0;                               //  the out-name there.
}

}  //  namespace merylutil::files::v1


////////////////////////////////////////
//
//  Removes the last dotted suffix from 'filename' and places the
//  result in 'basename'.
//
namespace merylutil::inline files::inline v1 {

void
stripLastSuffix(char *basename, char const *filename) {

  strcpy(basename, filename);               //  Copy in-name to out-name.

  if (directoryExists(basename))            //  Stop if this is an actual
    return;                                 //  directory.

  char  *slash = strrchr(basename, '/');    //  Search backwards for the last
  char  *dot   = strrchr(basename, '.');    //  slash and last dot.

  if (dot == nullptr)                       //  If not dot, we're done.
    return;

  if (slash < dot)                          //  If the dot is after the slash,
    *dot = 0;                               //  remove the suffix.
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Rename 'oldname' to 'newname'.  Silently succeeds if 'oldname'
//  doesn't exist.
//
namespace merylutil::inline files::inline v1 {

bool
rename(char const *oldname, char const *newname, bool fatal) {

  if (pathExists(oldname) == false)
    return false;

  errno = 0;
  if (::rename(oldname, newname) == -1)
    return fatalError(fatal, "renane()--  Failed to rename file '%s' to '%s': %s\n",
                      oldname, newname, strerror(errno));

  return true;
}

bool
rename(char const *oldprefix, char oldseparator, char const *oldsuffix,
       char const *newprefix, char newseparator, char const *newsuffix, bool fatal) {
  char const *oldpath = constructPathName(oldprefix, oldseparator, oldsuffix);
  char const *newpath = constructPathName(newprefix, newseparator, newsuffix);
  bool        success = true;

  if (pathExists(oldpath)) {
    errno = 0;
    if (::rename(oldpath, newpath) == -1)
      success = fatalError(fatal, "renane()--  Failed to rename file '%s' to '%s': %s\n",
                           oldpath, newpath, strerror(errno));
  }

  delete [] oldpath;
  delete [] newpath;

  return success;
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Create a symlink:
//   - failing if the source file doesn't exist
//   - silently succeeding if the destination path exists -- even if it is
//     not a link to the source file
//
namespace merylutil::inline files::inline v1 {

void
symlink(char const *pathToFile, char const *pathToLink) {

  if (pathExists(pathToFile) == false)
    fprintf(stderr, "symlink()-- Original file '%s' doesn't exist, won't make a link to nothing.\n",
            pathToFile), exit(1);

  if (pathExists(pathToLink) == true)
    return;

  errno = 0;
  if (::symlink(pathToFile, pathToLink) == -1)
    fprintf(stderr, "symlink()-- Failed to make link '%s' pointing to file '%s': %s\n",
            pathToLink, pathToFile, strerror(errno)), exit(1);
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Tests on existence of files and directories.
//

namespace merylutil::inline files::inline v1 {

bool
pathExists(char const *prefix, char separator, char const *suffix) {
  char const  *pathname = constructPathName(prefix, separator, suffix);
  bool         isp      = false;
  struct stat  s;

  if (stat(pathname, &s) == 0)
    isp = true;

  delete [] pathname;
  return (isp == true);
}

bool
fileExists(char const *prefix, char separator, char const *suffix,
           bool        writable) {
  char const  *pathname = constructPathName(prefix, separator, suffix);
  bool         isd      = true;
  bool         isw      = false;
  mode_t       w        = (S_IWUSR | S_IWGRP | S_IWOTH);
  struct stat  s;

  if (stat(pathname, &s) == 0) {                   //  If path exists...
    isd = (s.st_mode & S_IFDIR);                   //    Is a directory?
    isw = (s.st_mode & w) || (writable == false);  //    Is writable?  (or don't care)
  }

  delete [] pathname;
  return (isd == false) && (isw == true);
}


bool
directoryExists(char const *prefix, char separator, char const *suffix) {
  char const  *pathname = constructPathName(prefix, separator, suffix);
  bool         isd      = false;
  struct stat  s;

  if (stat(pathname, &s) == 0) {                   //  If path exists...
    isd = (s.st_mode & S_IFDIR);                   //    Is a directory?
  }

  delete [] pathname;
  return (isd == true);
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Create or remove directories.
//    mkdir() does nothing if the directory exists.
//    rmdir() does nothing if the directory doesn't exist.
//
namespace merylutil::inline files::inline v1 {

bool
mkdir(char const *dirname, bool fatal) {

  if ((dirname == nullptr) ||             //  Be generous and not fail
      (*dirname == 0))                    //  if nullptr or empty
    return true;                          //  name supplied.

  if (directoryExists(dirname) == true)   //  If it already exists,
    return true;                          //  we're done!

  errno = 0;
  if (::mkdir(dirname, S_IRWXU | S_IRWXG | S_IRWXO) == -1)
    return fatalError(fatal, "mkdir()--  Failed to create directory '%s': %s\n",
                      dirname, strerror(errno));

  return true;
}

bool
mkpath(char const *dirname, bool fatal) {
  char *dircopy = duplicateString(dirname);

  for (char *dir = strchr(dircopy, '/'); dir; dir = strchr(dir+1, '/')) {
    *dir = 0;
    mkdir(dircopy, false);
    *dir = '/';
  }

  delete [] dircopy;

  return mkdir(dirname, fatal);
}

bool
rmdir(char const *dirname, bool fatal) {

  if (directoryExists(dirname) == false)   //  If it doesn't exist,
    return true;                           //  we're done!

  errno = 0;
  if (::rmdir(dirname) == -1)
    return fatalError(fatal, "rmdir()--  Failed to remove directory '%s': %s\n",
                      dirname, strerror(errno));

  return true;
}

}  //  namespace merylutil::files::v1




////////////////////////////////////////
//
//  Remove a file, or do nothing if the file doesn't exist.
//
namespace merylutil::inline files::inline v1 {

bool
unlink(char const *path, bool fatal) {
  return unlink(path, 0, nullptr, fatal);
}

bool
unlink(char const *prefix, char separator, char const *suffix, bool fatal) {
  char const *pathname = constructPathName(prefix, separator, suffix);
  bool        success  = true;

  if (fileExists(pathname)) {
    errno = 0;
    if (::unlink(pathname) == -1)
      success = fatalError(fatal, "unlink()--  Failed to remove file '%s': %s\n",
                           pathname, strerror(errno));
  }

  delete [] pathname;
  return success;
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Permissions.
//
//  Remove ALL write bits from a given path.
//  Set write bits on a given path, relative to what is allowed in the umask.
//
namespace merylutil::inline files::inline v1 {

bool
makeReadOnly(char const *prefix, char separator, char const *suffix) {
  char const  *pathname = constructPathName(prefix, separator, suffix);
  bool         success  = false;
  struct stat  s;

  if (stat(pathname, &s) == 0) {              //  If file exists...
    mode_t w = S_IWUSR | S_IWGRP | S_IWOTH;   //    Create write-enable mask.
    mode_t m = (s.st_mode) & ~w;              //    Turn off write-enable bits.

    errno = 0;
    if (::chmod(pathname, m) == -1)
      fprintf(stderr, "WARNING: Failed to remove write permission from file '%s': %s\n", pathname, strerror(errno));
    else
      success = true;
  }

  delete [] pathname;
  return success;
}

bool
makeWritable(char const *prefix, char separator, char const *suffix) {
  char const  *pathname = constructPathName(prefix, separator, suffix);
  bool         success  = false;
  struct stat  s;

  if (stat(pathname, &s) == 0) {              //  If file exists...
    mode_t u = umask(0);                      //    Destructively read the current umask.
    mode_t w = S_IWUSR | S_IWGRP | S_IWOTH;   //    Create write-enable mask.
    mode_t m = (s.st_mode) | (w & ~u);        //    Turn on write, respecting umask.

    umask(u);                                 //    Restore default umask.

    errno = 0;                                //    Add allowed write bits to the file.
    if (::chmod(pathname, m) == -1)
      fprintf(stderr, "WARNING: Failed to add write permission to file '%s': %s\n", pathname, strerror(errno));
    else
      success = true;
  }

  delete [] pathname;
  return success;
}

}  //  namespace merylutil::files::v1




////////////////////////////////////////
//
//  Metadata
//
namespace merylutil::inline files::inline v1 {

off_t
sizeOfFile(char const *prefix, char separator, char const *suffix) {
  char const  *pathname = constructPathName(prefix, separator, suffix);
  struct stat  s;

  errno = 0;
  if (::stat(pathname, &s) == -1)
    fprintf(stderr, "Failed to stat() file '%s': %s\n", pathname, strerror(errno)), exit(1);

  delete [] pathname;
  return s.st_size;
}

off_t
sizeOfFile(FILE *file) {
  struct stat  s;

  errno = 0;
  if (::fstat(fileno(file), &s) == -1)
    fprintf(stderr, "Failed to stat() FILE*: %s\n", strerror(errno)), exit(1);

  return s.st_size;
}

uint64
timeOfFile(char const *prefix, char separator, char const *suffix) {
  muFileTime mt;

  mt.getTimeOfFile(prefix, separator, suffix);

  return mt.lastModifyTime().getTime_int64();
}

uint64
timeOfFile(FILE *file) {
  muFileTime mt;

  mt.getTimeOfFile(file);

  return mt.lastModifyTime().getTime_int64();
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  File position.
//
//  For fseek(), if the stream is already at the correct position, just
//  return.
//
//      Unless we're on FreeBSD.  For unknown reasons, FreeBSD fails
//      updating the seqStore with mate links.  It seems to misplace the
//      file pointer, and ends up writing the record to the wrong
//      location.  ftell() is returning the correct current location,
//      and so AS_PER_genericStore doesn't seek() and just writes to the
//      current position.  At the end of the write, we're off by 4096
//      bytes.
//          LINK 498318175,1538 <-> 498318174,1537
//          fseek()--  seek to 159904 (whence=0); already there
//          safeWrite()-- write nobj=1x104 = 104 bytes at position 159904
//          safeWrite()-- wrote nobj=1x104 = 104 bytes position now 164000
//          safeWrite()-- EXPECTED 160008, ended up at 164000
//
namespace merylutil::inline files::inline v1 {

off_t
ftell(FILE *stream) {

  errno = 0;
  off_t pos = ::ftello(stream);

  if ((errno == ESPIPE) || (errno == EBADF))   //  Not a seekable stream.
    return ((off_t)1) << 42;                   //  Return some goofy big number.

  if (errno)
    fprintf(stderr, "ftell()--  Failed with %s.\n", strerror(errno)), exit(1);

  return pos;
}

void
fseek(FILE *stream, off_t offset, int whence) {
  off_t   beginpos = merylutil::ftell(stream);

#if !defined __FreeBSD__ && !defined __osf__ && !defined __APPLE__
  if ((whence == SEEK_SET) && (beginpos == offset))
    return;
#endif  //  __FreeBSD__

  errno = 0;
  if (::fseeko(stream, offset, whence) == -1)
    fprintf(stderr, "fseek()--  Failed with %s.\n", strerror(errno)), exit(1);

  if (whence == SEEK_SET)
    assert(::ftell(stream) == offset);
}

}  //  namespace merylutil::files::v1





////////////////////////////////////////
//
//
//

namespace merylutil::inline files::inline v1 {

FILE *
openInputFile(char const *prefix,
              char        separator,
              char const *suffix,
              bool        doOpen) {
  if (doOpen == false)
    return nullptr;

  char const *pathname = constructPathName(prefix, separator, suffix);

  if (strcmp(pathname, "-") == 0) {
    delete [] pathname;
    return stdin;
  }

  errno = 0;
  FILE *F = fopen(pathname, "r");
  if (F == nullptr)
    fprintf(stderr, "Failed to open '%s' for reading: %s\n", pathname, strerror(errno)), exit(1);

  delete [] pathname;
  return F;
}

}  //  namespace merylutil::files::v1



namespace merylutil::inline files::inline v1 {

FILE *
openOutputFile(char const *prefix,      //  The unlink() below is to prevent race conditions
               char        separator,   //  when two processes open the same file; we've seen
               char const *suffix,      //  instances where the outputs are intermingled.  This
               bool        doOpen) {    //  isn't perfect; they could still race, but the window
  if (doOpen == false)                  //  is much smaller now.
    return nullptr;

  char const *pathname = constructPathName(prefix, separator, suffix);
  FILE       *F        = stdout;

  if (strcmp(pathname, "-") != 0) {
    unlink(pathname);

    errno = 0;
    F = fopen(pathname, "w");
    if (F == nullptr)
      fprintf(stderr, "Failed to open '%s' for writing: %s\n", pathname, strerror(errno)), exit(1);
  }

  delete [] pathname;
  return F;
}

void
closeFile(FILE *&F, char const *prefix, char separator, char const *suffix, bool critical) {

  if ((F == nullptr) || (F == stdout) || (F == stderr))
    return;

  errno = 0;
  fclose(F);   F = nullptr;

  if ((critical == false) || (errno == 0))
    return;

  if ((prefix) && (suffix))
    fprintf(stderr, "Failed to close file '%s%c%s': %s\n", prefix, separator, suffix, strerror(errno));
  else if (prefix)
    fprintf(stderr, "Failed to close file '%s': %s\n", prefix, strerror(errno));
  else
    fprintf(stderr, "Failed to close file: %s\n", strerror(errno));

  exit(1);
}

void
closeFile(FILE *&F, char const *pathname, bool critical) {
  closeFile(F, pathname, '.', nullptr, critical);
}

}  //  namespace merylutil::files::v1



////////////////////////////////////////
//
//  Create an empty file.  The local variable is because closeFile() takes a
//  reference to FILE*, which openFile() can't supply.
//
namespace merylutil::inline files::inline v1 {

void
createEmptyFile(char const *prefix, char separator, char const *suffix) {
  FILE *F = openOutputFile(prefix, separator, suffix);
  closeFile(F, prefix, separator, suffix);
}

}  //  namespace merylutil::files::v1


