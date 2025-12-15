
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

#include "time-v1.H"
#include "files/accessing-v1.H"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>

namespace merylutil::inline system::inline v1 {

muTime &
muTime::setTime(void) {
  struct timeval  tp;

  gettimeofday(&tp, nullptr);

  _s = tp.tv_sec;
  _n = tp.tv_usec * 1000;

  return *this;
}





bool
muFileTime::setTimeOfFile(char const *prefix, char separator, char const *suffix) {
  char const *pathname = constructPathName(prefix, separator, suffix);
  bool        success  = false;

  struct timespec times[2] = {{ .tv_sec = _acc._s, .tv_nsec = _acc._n },    //  Access time
                              { .tv_sec = _mod._s, .tv_nsec = _mod._n }};   //  Modification time

  if (::utimensat(AT_FDCWD, pathname, times, 0) == -1)
    fprintf(stderr, "WARNING: Failed to set time of file '%s': %s\n", pathname, strerror(errno));
  else
    success = true;

  delete [] pathname;

  return success;
}


muTime
muFileTime::getTimeOfFile(char const *prefix, char separator, char const *suffix) {
  struct stat  s;

  char const *pathname = constructPathName(prefix, separator, suffix);

  errno = 0;
  if (stat(pathname, &s) == -1)
    fprintf(stderr, "Failed to stat() file '%s': %s\n", pathname, strerror(errno)), exit(1);

  delete [] pathname;

#ifdef __APPLE__
  _mod._s = s.st_mtimespec.tv_sec;    _acc._s = s.st_atimespec.tv_sec;
  _mod._n = s.st_mtimespec.tv_nsec;   _acc._n = s.st_atimespec.tv_nsec;
#else
  _mod._s = s.st_mtim.tv_sec;         _acc._s = s.st_atim.tv_sec;
  _mod._n = s.st_mtim.tv_nsec;        _acc._n = s.st_atim.tv_nsec;
#endif

  return _mod;
}


muTime
muFileTime::getTimeOfFile(FILE *file) {
  struct stat  s;
  off_t        size = 0;

  errno = 0;
  if (fstat(fileno(file), &s) == -1)
    fprintf(stderr, "Failed to stat() FILE*: %s\n", strerror(errno)), exit(1);

#ifdef __APPLE__
  _mod._s = s.st_mtimespec.tv_sec;    _acc._s = s.st_atimespec.tv_sec;
  _mod._n = s.st_mtimespec.tv_nsec;   _acc._n = s.st_atimespec.tv_nsec;
#else
  _mod._s = s.st_mtim.tv_sec;         _acc._s = s.st_atim.tv_sec;
  _mod._n = s.st_mtim.tv_nsec;        _acc._n = s.st_atim.tv_nsec;
#endif

  return _mod;
}

}  //  namespace merylutil::system::v1
