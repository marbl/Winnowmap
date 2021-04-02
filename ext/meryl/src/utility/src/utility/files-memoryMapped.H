
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

#ifndef FILES_MEMORY_MAPPED_H
#define FILES_MEMORY_MAPPED_H

//  Do not include directly.  Use 'files.H' instead.


//  The BSD's are able to map to an arbitrary position in the file, but the Linux's can only map to
//  multiples of pagesize.  Thus, this class maps the whole file into address space, then passes out
//  pointers to pieces in it.  This is slightly unfortunate, because array out-of-bounds will not be
//  caught.  To be fair, on the BSD's the file is mapped to a length that is a multiple of pagesize,
//  so it would take a big out-of-bounds to fail.

enum memoryMappedFileType {
  memoryMappedFile_readOnly        = 0x00,
  memoryMappedFile_readOnlyInCore  = 0x01,
  memoryMappedFile_readWrite       = 0x02,
  memoryMappedFile_readWriteInCore = 0x03
};


#ifndef MAP_POPULATE
#define MAP_POPULATE 0
#endif

class memoryMappedFile {
public:
  memoryMappedFile(const char           *name,
                   memoryMappedFileType  type = memoryMappedFile_readOnly);
  ~memoryMappedFile();

  //  get(size_t offset, size_t length) returns 'length' bytes starting starting at position
  //  'offset'.  The current position is updated to 'offset + length'.
  //
  //  get(length) returns 'length' bytes starting at the current position.  The current position
  //  is updated to the byte after the last one returned.
  //
  //  get(), get(0) and get(offset, 0) all return the current position.

  void  *get(size_t offset, size_t length) {

    if (offset + length > _length)
      fprintf(stderr, "memoryMappedFile()-- Requested " F_SIZE_T " bytes at position " F_SIZE_T " in file '%s', but only " F_SIZE_T " bytes in file.\n",
              length, offset, _name, _length), exit(1);

    _offset = offset + length;

    return((uint8 *)_data + offset);
  };

  void                  *get(size_t length=0)  { return(get(_offset, length)); };
  size_t                 length(void)          { return(_length);              };
  memoryMappedFileType   type(void)            { return(_type);                };


private:
  char                    _name[FILENAME_MAX];

  memoryMappedFileType    _type;

  size_t                  _length;  //  Length of the mapped file
  size_t                  _offset;  //  File pointer for reading

  int32                   _fd;
  void                   *_data;
};

#endif  //  FILES_MEMORYMAPPED_H
