
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

#include <fcntl.h>
#include <sys/mman.h>

#include "files.H"

namespace merylutil::inline files::inline v1 {

memoryMappedFile::memoryMappedFile(const char *name,
                                   mftType     type) {

  strncpy(_name, name, FILENAME_MAX-1);

  _type = type;

  errno = 0;
  _fd = (_type == mftReadOnly) ? open(_name, O_RDONLY | O_LARGEFILE)
                               : open(_name, O_RDWR   | O_LARGEFILE);
  if (_fd < 0)
    fprintf(stderr, "memoryMappedFile()-- Couldn't open '%s' for mapping: %s\n", _name, strerror(errno)), exit(1);

  struct stat  sb;

  if (fstat(_fd, &sb) < 0)
    fprintf(stderr, "memoryMappedFile()-- Couldn't stat '%s' for mapping: %s\n", _name, strerror(errno)), exit(1);

  _length = sb.st_size;
  _offset = 0;

  if (_length == 0)
    fprintf(stderr, "memoryMappedFile()-- File '%s' is empty, can't map.\n", _name), exit(1);

  //  Map the file to memory, or grab some anonymous space for the file to be copied to.

  if (_type == mftReadOnly)
    _data = mmap(0L, _length, PROT_READ,              MAP_FILE | MAP_PRIVATE, _fd, 0);

  if (_type == mftReadOnlyInCore)
    _data = mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_ANON | MAP_PRIVATE, -1, 0);

  if (_type == mftReadWrite)
    _data = mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED, _fd, 0);

  if (_type == mftReadWriteInCore)
    _data = mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);

  if (_data == MAP_FAILED)
    fprintf(stderr, "memoryMappedFile()-- Failed to map file '%s': %s\n", _name, strerror(errno)), exit(1);

  //  If loading into core, read the file into core.

  if ((_type == mftReadOnlyInCore) ||
      (_type == mftReadWriteInCore)) {
    ssize_t ar = read(_fd, _data, _length);

    if (ar < 0)
      fprintf(stderr, "memoryMappedFile()-- Failed to read %lu bytes from mapped file '%s': %s\n", _length, _name, strerror(errno)), exit(1);
    if (ar != _length)
      fprintf(stderr, "memoryMappedFile()-- Read only %lu bytes out of expected %lu bytes from mapped file '%s': %s\n", ar, _length, _name, strerror(errno)), exit(1);
  }

  //  Close the file if we're done with it.

  if (_type != mftReadWriteInCore) {
    if (close(_fd) < 0)
      fprintf(stderr, "memoryMappedFile()-- Failed to close file '%s' after mapping: %s\n", _name, strerror(errno)), exit(1);

    _fd = -1;
  }
}


memoryMappedFile::~memoryMappedFile() {

  errno = 0;

  if (_type == mftReadWrite)
    if (msync(_data, _length, MS_SYNC) < 0)
      fprintf(stderr, "memoryMappedFile()-- Failed to sync mapped file '%s' of length " F_SIZE_T ": %s\n", _name, _length, strerror(errno)), exit(1);

  if (_type == mftReadWriteInCore) {
    ssize_t aw = write(_fd, _data, _length);

    if (aw < 0)
      fprintf(stderr, "memoryMappedFile()-- Failed to write %lu bytes to file '%s': %s\n", _length, _name, strerror(errno)), exit(1);
    if (aw != _length)
      fprintf(stderr, "memoryMappedFile()-- Wrote only %lu bytes out of expected %lu bytes to file '%s': %s\n", aw, _length, _name, strerror(errno)), exit(1);

    if (close(_fd) < 0)
      fprintf(stderr, "memoryMappedFile()-- Failed to close file mapped file '%s': %s\n", _name, strerror(errno)), exit(1);
  }

  if (munmap(_data, _length) < 0)
    fprintf(stderr, "memoryMappedFile()-- Failed to unmap file '%s': %s\n", _name, strerror(errno)), exit(1);
}

}  //  merylutil::files::v1
