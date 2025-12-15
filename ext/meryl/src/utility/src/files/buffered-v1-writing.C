
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

#include "arrays.H"
#include "files.H"

namespace merylutil::inline files::inline v1 {


void
writeBuffer::initialize(const char *pfx, char sep, const char *sfx, const char *mode, uint64 bMax) {

  //  Figure out what file to write to.  Unlike readBuffer, we do not allow
  //  writing to stdout, so our decision is based soley on the existence of
  //  a suffix.

  if (sfx == nullptr)
    strncpy(_filename, pfx, FILENAME_MAX);
  else
    snprintf(_filename, FILENAME_MAX, "%s%c%s", pfx, sep, sfx);

  //  Remember the mode.

  strncpy(_filemode, mode, 16);

  //  If appending, open the file now so we can set the file position.
  //  Otherwise, if the mode isn't for writing, fail.

  if      (mode[0] == 'a')
    open();
  else if (mode[0] != 'w')
    fprintf(stderr, "writeBuffer()--  Unknown mode '%s'\n", mode), exit(1);

  //  Allocate a buffer.

  _bufferMax      = (bMax == 0) ? 32 * 1024 : bMax;
  _buffer         = new uint8 [_bufferMax];
}



writeBuffer::~writeBuffer() {
  flush();

  delete [] _buffer;

  delete [] _chunkBuffer;

  delete [] _chunkStarts;
  delete [] _chunkSizes;

  closeFile(_file, _filename);
}



void
writeBuffer::write(void const *data, uint64 length) {

  if (_bufferMax < _bufferLen + length)           //  Flush the buffer if this
    flush();                                      //  data is too big for it.

  if (_bufferMax < length) {                      //  And if it is still too big
    assert(_bufferLen == 0);                      //  (ensure the buffer is empty)
    writeToDisk(data, length);                    //  and just dump it to disk.
  }

  else {                                          //  Otherwise, copy it to
    memcpy(_buffer + _bufferLen, data, length);   //  our buffer.
    _bufferLen += length;
  }

  assert(_bufferLen <= _bufferMax);

  _filePos += length;
}



//  If dataLength is zero, only a header is written, and the chunk stack
//  is increased.
//
//  If dataLength is non-zero, and the stack is empty, the block
//  is immediately written to disk.
//
//  If dataLength is non-zero, and the stack has entries, the block
//  is copied to our internal chunk buffer and the size added
//  to ALL container chunks.
//
void
writeBuffer::writeIFFchunk(char const *name, void *data, uint32 dataLength) {
  uint8  header [4 * sizeof(uint8) + sizeof(uint32)] = { 0 };
  uint8  padding[4 * sizeof(uint8)]                  = { 0 };

  //  Figure out how much padding we need to add to the data to make it
  //  align on a 32-bit boundary.

  uint32 headLength = 4 * sizeof(uint8) + sizeof(uint32);
  uint32 padLength  = 4 - (dataLength % 4);

  if (padLength == 4)
    padLength = 0;

  //  Create the chunk header.

  dataLength += padLength;

  memcpy(header + 0,  name,       sizeof(uint8) * 4);
  memcpy(header + 4, &dataLength, sizeof(uint32));

  dataLength -= padLength;

  //  If the chunk size is zero, add a new container chunk (adding space for
  //  8 more if needed).
  if (dataLength == 0) {
    increaseArrayPair(_chunkStarts, _chunkSizes, _chunkStartsLen, _chunkStartsMax, 8);

    _chunkStarts[_chunkStartsLen] = _chunkBufferLen;
    _chunkSizes [_chunkStartsLen] = 0;

    _chunkStartsLen++;

    appendIFFdata(header, headLength);
  }

  //  But if there is no container, we can immediately output the chunk.
  else if (_chunkStartsLen == 0) {
    write(header,  headLength);
    write(data,    dataLength);
    write(padding, padLength);
  }

  //  Otherwise, append the chunk to our buffer and increase
  //  the size of each parent container.
  else {
    appendIFFdata(header,  headLength);
    appendIFFdata(data,    dataLength);
    appendIFFdata(padding, padLength);

    for (uint32 ii=0; ii<_chunkStartsLen; ii++)
      _chunkSizes[ii] += headLength + dataLength + padLength;
  }
}



void
writeBuffer::appendIFFdata(void *data, uint32 dataLength) {

  if (dataLength == 0)
    return;

  //  If this data will exceed our current allocation, double what we've
  //  allocated until it'll fit.

  if (_chunkBufferLen + dataLength > _chunkBufferMax) {
    uint64  newMax = (_chunkBufferMax == 0) ? 16384 : _chunkBufferMax;

    while (newMax < _chunkBufferLen + dataLength)
      newMax *= 2;

    resizeArray(_chunkBuffer, _chunkBufferLen, _chunkBufferMax, newMax);
  }

  //  Copy the data into the buffer and update the length.

  memcpy(_chunkBuffer + _chunkBufferLen, data, dataLength);

  _chunkBufferLen += dataLength;
}



//  The user is done with this chunk.
//
//  If the stack has entries, set the length of the last
//  chunk and pop the stack.  If the stack is now
//  empty, write the block.
//
void
writeBuffer::closeIFFchunk(char const *name) {

  //  If no chunk to close, report an error.

  if (_chunkStartsLen == 0) {
    fprintf(stderr, "writeBuffer::closeIFFchunk()-- no chunk to close.\n");
    exit(1);
    return;
  }

  //  Refer to the last chunkStarts entry.

  _chunkStartsLen--;

  //  If a name supplied, check that it's the same as the chunk we're
  //  closing.

  if (name) {
    uint64  cs = _chunkStarts[_chunkStartsLen];

    if ((name[0] != _chunkBuffer[cs + 0]) ||
        (name[1] != _chunkBuffer[cs + 1]) ||
        (name[2] != _chunkBuffer[cs + 2]) ||
        (name[3] != _chunkBuffer[cs + 3])) {
      fprintf(stderr, "writeBuffer::closeIFFchunk()-- requested to close chunk '%c%c%c%c' but current chunk is '%c%c%c%c'.\n",
              name[0], name[1], name[2], name[3],
              _chunkBuffer[cs + 0],
              _chunkBuffer[cs + 1],
              _chunkBuffer[cs + 2],
              _chunkBuffer[cs + 3]);
      exit(1);
    }
  }

  //  Update the size of the chunk container we're in.

  *(uint32 *)(_chunkBuffer + _chunkStarts[_chunkStartsLen] + 4) = _chunkSizes[_chunkStartsLen];

  //  If there are no more containers, write this buffer to disk, and clear it.

  if (_chunkStartsLen == 0) {
    write(_chunkBuffer, _chunkBufferLen);

    _chunkBufferLen = 0;
  }
}




void
writeBuffer::open(void) {
  if (_file != nullptr)
    return;

  errno = 0;
  _file = fopen(_filename, _filemode);
  if (errno)
    fprintf(stderr, "writeBuffer()--  Failed to open file '%s' with mode '%s': %s\n",
            _filename, _filemode, strerror(errno)), exit(1);

  //  If appending, _filePos is zero, and ftell() is non-zero.
  //  If writing, _filePos is non-zero, and ftell() is zero.
  _filePos += merylutil::ftell(_file);
}



void
writeBuffer::writeToDisk(void const *data, uint64 length) {
  if (length == 0)
    return;

  open();
  writeToFile((char *)data, "writeBuffer::writeToDisk", length, _file);
}



void
writeBuffer::flush(void) {
  writeToDisk(_buffer, _bufferLen);
  _bufferLen = 0;
}



}  //  merylutil::files::v1

