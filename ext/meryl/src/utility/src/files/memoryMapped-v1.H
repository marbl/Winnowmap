
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

#ifndef MERYLUTIL_FILES_MEMORYMAPPED_V1_H
#define MERYLUTIL_FILES_MEMORYMAPPED_V1_H

#include "types.H"

#include <sys/stat.h>


//  The BSD's are able to map to an arbitrary position in the file, but the
//  Linux's can only map to multiples of pagesize.  Thus, this class maps the
//  whole file into address space, then passes out pointers to pieces in it.
//  This is slightly unfortunate, because array out-of-bounds will not be
//  caught.  To be fair, on the BSD's the file is mapped to a length that is
//  a multiple of pagesize, so it would take a big out-of-bounds to fail.
//
//  get(size_t offset, size_t length)
//    returns a pointer to 'length' bytes starting at position 'offset'.  The
//    current position is updated to 'offset + length'.
//
//  get(length)
//    returns a pointer to 'length' bytes starting at the current position.
//    The current position is updated to 'current_position + length'.
//
//  When 'length' is zero, the current position remains unchanged:
//    get() and get(0) return the current positon.
//    get(offset, 0) returns 'offset'.
//

namespace merylutil::inline files::inline v1 {

enum mftType {
  mftReadOnly        = 0x00,
  mftReadOnlyInCore  = 0x01,
  mftReadWrite       = 0x02,
  mftReadWriteInCore = 0x03
};


class memoryMappedFile {
public:
  memoryMappedFile(const char *name,
                   mftType     type = mftReadOnly);
  ~memoryMappedFile();

  void      *get(size_t offset,
                 size_t length);
  void      *get(size_t length=0)  { return(get(_offset, length)); };
  size_t     length(void)          { return(_length);              };
  mftType    type(void)            { return(_type);                };

private:
  char       _name[FILENAME_MAX] = {0};

  mftType    _type = mftReadOnly;

  size_t     _length = 0;   //  Length of the mapped file
  size_t     _offset = 0;   //  File pointer for reading

  int32      _fd   = -1;
  void      *_data = nullptr;
};


inline
void *
memoryMappedFile::get(size_t offset, size_t length) {

  if (offset + length > _length)
    fprintf(stderr, "memoryMappedFile()-- Requested " F_SIZE_T " bytes at position " F_SIZE_T " in file '%s', but only " F_SIZE_T " bytes in file.\n",
            length, offset, _name, _length), exit(1);

  _offset = offset + length;

  return((uint8 *)_data + offset);
}

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_FILES_MEMORYMAPPED_V1_H
