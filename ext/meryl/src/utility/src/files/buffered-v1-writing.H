
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

#ifndef MERYLUTIL_FILES_BUFFERED_V1_WRITING_H
#define MERYLUTIL_FILES_BUFFERED_V1_WRITING_H

#include "types.H"

namespace merylutil::inline files::inline v1 {

class writeBuffer {
public:
  writeBuffer(const char *pfx, char sep, const char *sfx, const char *mode, uint64 bMax = 1024 * 1024) { initialize(pfx, sep, sfx, mode, bMax); }
  writeBuffer(const char *filename,                       const char *mode, uint64 bMax = 1024 * 1024) { initialize(filename, '.', nullptr, mode, bMax); }
  ~writeBuffer();

private:
  void                 initialize(const char *pfx, char sep, const char *sfx, const char *mode, uint64 bMax);

public:
  const char          *filename(void) { return(_filename); }
  uint64               tell(void)     { return(_filePos);  }

  void                 write(void const *data, uint64 length);
  void                 flush(void);

  //  If deata length is zero, this chunk is the start of a recursive record.
  //  The chunk header is written (NAME<length>) but length is set to zero.
  //  The position of the length field is is pushed onto an internal stack.
  //  When the chunk is closed, the length field is updated with
  //  the size of blocks added to it.
  //
  //  To get an empty chunk:
  //    writeIFFchunk("NAME", NULL, 0);
  //    closeIFFchunk();
  //
private:
  void   appendIFFdata(void *data, uint32 dataLength);

public:
  template<typename OBJ>
  void   writeIFFarray(char const *name, OBJ *array, uint32 nObjects) {
    writeIFFchunk(name, array, sizeof(OBJ) * nObjects);
  }

  template<typename OBJ>
  void   writeIFFobject(char const *name, OBJ &object) {
    writeIFFchunk(name, &object, sizeof(OBJ));
  }

  template<typename OBJ>
  void   writeIFFobject(char const *name, OBJ *object) {
    writeIFFchunk(name, object, sizeof(OBJ));
  }

  void   writeIFFchunk(char const *name, void *data=NULL, uint32 dataLength=0);

  void   closeIFFchunk(char const *name=NULL);

private:
  void                 open(void);
  void                 writeToDisk(void const *data, uint64 length);

  char                _filename[FILENAME_MAX+1] = {0};
  char                _filemode[17]             = {0};

  FILE               *_file           = nullptr;
  uint64              _filePos        = 0;

  uint64              _bufferLen      = 0;
  uint64              _bufferMax      = 0;
  uint8               *_buffer         = nullptr;

  uint64              _chunkBufferLen = 0;         //  For building up recursive chunks,
  uint64              _chunkBufferMax = 0;         //  another buffer of data.
  uint8              *_chunkBuffer    = nullptr;

  uint32              _chunkStartsLen = 0;         //  A stack of start positions.
  uint32              _chunkStartsMax = 0;
  uint64             *_chunkStarts    = nullptr;
  uint64             *_chunkSizes     = nullptr;
};

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_FILES_BUFFERED_V1_WRITING_H
