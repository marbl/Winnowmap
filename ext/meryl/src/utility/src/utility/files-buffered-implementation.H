
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

//  Returns the next letter in the buffer, but DOES NOT advance past
//  it.  Might have some wierd interaction with EOF -- if you peek()
//  and the next thing is eof , the _eof flag might get set.
//
inline
char
readBuffer::peek(void) {

  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  return(_buffer[_bufferPos]);
}


//  Returns the next letter in the buffer and advances file position by one.
//  Returns 0 on EOF.
//
inline
char
readBuffer::read(void) {

 readAgain:
  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  _bufferPos++;
  _filePos++;

  if ((_ignoreCR) && (_buffer[_bufferPos-1] == '\r'))
    goto readAgain;

  return(_buffer[_bufferPos-1]);
}



//  Like read(), but does not advance file position if
//  the 'stop' character is encountered.
//  Returns 'stop' on EOF.
//
inline
char
readBuffer::readuntil(char stop) {
  char  ch;

 readuntilAgain:
  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  ch = _buffer[_bufferPos];

  if (ch == stop)
    return(stop);

  _bufferPos++;
  _filePos++;

  if ((_ignoreCR) && (ch == '\r'))
    goto readuntilAgain;

  return(ch);
}


//  Advances the file position to the next 'stop' character and returns
//  'stop', or 0 if eof.
//
//  The file is left on the 'stop' character if after == false,
//  otherwise, the file is positioned after the 'stop' character.
//
inline
void
readBuffer::skipAhead(char stop, bool after) {

  while (_eof == false) {

    //  Insert a fake stop at the end of the buffer.
    _buffer[_bufferLen] = stop;

    //  Scan the buffer for the next real or fake stop.
    while (_buffer[_bufferPos] != stop) {
      _bufferPos++;
      _filePos++;
    }

    //  If we hit the fake stop, fill the buffer again and continue.
    if (_bufferPos == _bufferLen) {
      fillBuffer();
      continue;
    }

    //  Otherwise, we've found a real stop.  Skip it if desired.
    if (after) {
      _bufferPos++;
      _filePos++;
    }

    //  Fill the buffer again, if needed.
    if (_bufferPos == _bufferLen)
      fillBuffer();

    return;
  }
}


//  Advances the file position to the next 'stop' character, copying
//  the characters into a user buffer.
//
//  The file is left on the 'stop' character.
//
//  Returns the number of letters copied.
//
inline
uint64
readBuffer::copyUntil(char stop, char *dest, uint64 destLen) {
  uint64  copied = 0;

  while (_eof == false) {
    while ((_buffer[_bufferPos] != stop) &&
           (_bufferPos < _bufferLen) &&
           (copied < destLen)) {
      dest[copied] = _buffer[_bufferPos];

      copied++;
      _bufferPos++;
      _filePos++;
    }

    if (_bufferPos < _bufferLen)
      return(copied);

    fillBuffer();
  }

  return(copied);
}

