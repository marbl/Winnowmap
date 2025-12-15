
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


//  Open a readBuffer from stdin, '{pfx}' or '{pfx}{sep}{sfx}'.

void
readBuffer::initialize(const char *pfx, char sep, const char *sfx, uint64 bMax) {

  if ((pfx == nullptr) ||                                          //  Read from stdin if no pfx supplied
      ((pfx[0] == '-') && (pfx[1] == 0))) {                        //  or if pfx is '-'.
    strcpy(_fn, "(stdin)");
    _stdin = true;
  }
  else if ((pfx != nullptr) && (sfx == nullptr)) {                 //  Read from file '{pfx}' if there
    strncpy(_fn, pfx, FILENAME_MAX);                               //  is no suffix supplied.
    _owned = true;
  }
  else if ((pfx != nullptr) && (sep != 0) && (sfx != nullptr)) {   //  Read from file {pfx}{sep}{sfx}
    snprintf(_fn, FILENAME_MAX, "%s%c%s", pfx, sep, sfx);          //  if both those are supplied.
    _owned = true;
  }
  else {
    fprintf(stderr, "readBuffer()-- Invalid filename pfx='%s' sep='%d' sfx='%s'\n", pfx, sep, sfx);
    exit(1);
  }

  _bMax   = (bMax == 0) ? 32 * 1024 : bMax;                        //  Allocate a buffer.
  _b      = new uint8 [_bMax + 1];

  errno = 0;                                                       //  Open the file, failing if
  _f = (_stdin) ? fileno(stdin)                                    //  it is the console.
                   : ::open(_fn, O_RDONLY | O_LARGEFILE);
  if (_f == -1)
    fprintf(stderr, "readBuffer()-- couldn't open file '%s': %s\n",
            _fn, strerror(errno)), exit(1);

  if (isatty(_f))
    fprintf(stderr, "readBuffer()-- refuse to use the terminal for input; provide a filename or input from a pipe.\n"), exit(1);
  
  fillBuffer();
}



readBuffer::readBuffer(FILE *file, uint64 bMax) {

  strcpy(_fn, "(hidden file)");

  _bMax   = (bMax == 0) ? 32 * 1024 : bMax;                       //  Allocate a buffer.
  _b      = new uint8 [_bMax + 1];

  _f        = fileno(file);                                       //  Get the file handle.  It never fails.

  errno = 0;                                                      //  Rewind the file, allowing failure
  if ((::lseek(_f, 0, SEEK_SET) == -1) && (errno != ESPIPE))      //  if it's stdin or a pipe.
    fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position 0: %s\n",
            _fn, strerror(errno)), exit(1);

  fillBuffer();
}



void
readBuffer::fillBufferImpl(void) {

  assert(_bBgn + _bPos == _fPos);

  if (_bPos < _bLen) {            //  If still have stuff, save it.
    uint64 bl = _bLen - _bPos;

    memmove(_b, _b + _bPos, bl);

    _bBgn += _bPos;               //    Set begin of this buffer to the end of the last.
    _bPos  = 0;                   //    Position in the buffer is zero.
    _bLen  = bl;                  //    Length of the valid data is also zero.
  }
  else {                          //  Empty or exhausted buffer:
    assert(_bPos == _bLen);       //    Must be exactly empty.
    _bBgn += _bLen;               //    Set begin of this buffer to the end of the last.
    _bPos  = 0;                   //    Position in the buffer is zero.
    _bLen  = 0;                   //    Length of the valid data is also zero.
  }

  assert(_bBgn + _bPos == _fPos);

 again:  //  See read() below
  ssize_t  r = ::read(_f, _b + _bLen, _bMax - _bLen);

  if (r < 0) {                    //  Fail if an error, unless that error
    if (errno == EAGAIN)          //  is because no data was ready to be
      goto again;                 //  returned.
    else
      fprintf(stderr, "readBuffer::fillBuffer()-- couldn't read " F_U64 " bytes from '%s': %s\n",
              _bMax, _fn, strerror(errno)), exit(1);
  }
  else if (r == 0) {              //  EOF if no data returned.
    _eof = (_bPos == _bLen);      //  If we were continuing a read from before
                                  //  then we have unprocessed data in the buffer so we're not at the end
                                  //  the next call to this function to re-fill the buffer will mark it as EOF
    assert((!_eof && _bPos < _bLen) || (_eof && _bPos == _bLen));
  }
  else {                          //  Otherwise, we got data.
    _eof = false;                 //  (and therefore, not eof)
    _bLen += (uint64)r;
    assert(_bPos < _bLen);
    assert(_bLen > 0);
    if (_bLen < _bMax)            // when we don't get what we requested for the buffer, try again until we do
       goto again;
  }

  assert(_bBgn + _bPos == _fPos);
}


//  TEST:  What happens if seek to last position in buffer,
//  when that is the end of file?
void
readBuffer::seek(uint64 pos, uint64 extra) {

  if ((_bBgn <= pos) &&                   //  New position is in the
      (pos + extra <= _bBgn + _bLen)) {   //  existing buffer; no need
    _fPos = pos;                          //  to read from disk.
    _bPos = pos - _bBgn;
  }
  else {                                  //  Need more data!
    if (::lseek(_f, pos, SEEK_SET) == -1)
      fprintf(stderr, "readBuffer()-- '%s' couldn't seek to position " F_U64 ": %s\n",
              _fn, pos, strerror(errno)), exit(1);

    _fPos = pos;                          //  Set internal positions,
    _bBgn = pos;                          //  to 'pos', empty our buffer
    _bPos = 0;                            //  and reload it from the
    _bLen = 0;                            //  current file position.

    fillBuffer();
  }
}



uint64
readBuffer::read(void *buf, uint64 len) {
  char  *bufchar = (char *)buf;

  //  Easy case; the next len bytes are already in the buffer; just
  //  copy and move the position.

  if (_bPos + len <= _bLen) {
    memcpy(bufchar, _b + _bPos, len);

    _fPos += len;
    _bPos += len;

    fillBuffer();

    return(len);
  }

  //  Existing buffer not big enough.  Copy what's there, then finish
  //  with a direct read from disk and then refill the buffer.

  uint64   bCopied = 0;     //  Number of bytes copied into the buffer
  uint64   bAct    = 0;     //  Number of bytes actually read from disk

  bCopied     = _bLen - _bPos;

  memcpy(bufchar, _b + _bPos, bCopied);

  while (bCopied < len) {  //  See fillBuffer() above.
    ssize_t  r = ::read(_f, bufchar + bCopied, len - bCopied);

    if (r < 0) {
      if (errno == EAGAIN)
        continue;
      else
        fprintf(stderr, "readBuffer::fillBuffer()-- couldn't read " F_U64 " bytes from '%s': %s\n",
                len, _fn, strerror(errno)), exit(1);
    }
    else if (r == 0)        //  EOF if no data left.  Return whatever we
      len = 0;              //  read by declaring the desired length to be zero.
    else                    //  Otherwise, we got data.
      bCopied += (uint64)r;
  }

  _fPos += bCopied;         //  Advance the actual file position to however much we just read.
  _bBgn  = _fPos;           //  And set the buffer begin to that too.
  _bPos  = 0;
  _bLen  = 0;               //  Set the buffer as empty, so we fill it again.

  fillBuffer();

  return(bCopied);
}






bool
readBuffer::peekIFFchunk(char name[4], uint32 &dataLen) {

  //  Seek to the current position, making sure there are at least
  //  8 bytes still in the buffer.

  seek(_fPos, 8);

  //  If there's space for a valid IFF header, return the name and length.

  if (_bPos + 8 < _bLen) {
    memcpy( name,    _b + _bPos,     sizeof(char) * 4);
    memcpy(&dataLen, _b + _bPos + 4, sizeof(uint32));

    return(true);
  }

  //  If not, return an empty name and length of zero.

  name[0] = 0;
  name[1] = 0;
  name[2] = 0;
  name[3] = 0;

  dataLen = 0;

  return(false);
}



//  Read a specific chunk from the buffer.
//  Return false if the chunk isn't named 'name' and of length 'dataLen'.
//
bool
readBuffer::readIFFchunk(char const *name,
                         void       *data,
                         uint32      dataLen) {
  char    dtag[4] = {0};
  uint32  dlen    =  0;

  if (peekIFFchunk(dtag, dlen) == false)
    return(false);

  if ((dtag[0] != name[0]) ||
      (dtag[1] != name[1]) ||
      (dtag[2] != name[2]) ||
      (dtag[3] != name[3]) ||
      (dlen    != dataLen))
    return(false);

  //  It's the one we want, so read the data for real.

  uint32   rl = 0;

  rl += read( dtag, 4);
  rl += read(&dlen, sizeof(uint32));
  rl += read( data, dataLen);

  return(rl == 4 + sizeof(uint32) + dataLen);
}



void
readBuffer::readIFFchunk(char*name, uint8 *&data, uint32 &dataLen, uint32 &dataMax) {

  //  Read the name and data length.

  read( name,    4);
  read(&dataLen, sizeof(uint32));

  //  Allocate space for the data.

  resizeArray(data, 0, dataMax, dataLen);

  //  Copy the data to 'data'.

  read(data, dataLen);
}




}  //  merylutil::files::v1

