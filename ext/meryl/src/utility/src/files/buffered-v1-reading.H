
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

#ifndef MERYLUTIL_FILES_BUFFERED_V1_READING_H
#define MERYLUTIL_FILES_BUFFERED_V1_READING_H

#include "types.H"
#include "strings.H"  //  debug

//
//  Reads a file using a small-ish buffer.
//  Allows 'free' seeks within the buffer.
//  Provides useful function for copying data out and skipping over data.
//
//
//    [ descriptions largely correct but function names incorrect ]
//
//
//  skip(charclass, sense, next)
//   - advances buffer position while current letter is in the
//     supplied character class (sense=true) or until the letter is in
//     the character class (sense=false; alternative explanation is to advance
//     while the letter is not in the class).
//   - if next is true, will then continue skipping while the sense is
//     inverted.
//
//  skipUntil(stop, next=false)
//   - advances buffer position while the current letter is not a 'stop'
//     letter.  does nothing if current letter is a 'stop' letter.
//   - returns letter of resulting buffer position, or NUL if EOF.
//   - if 'next', will further advance buffer position to the next
//     non-'stop' letter.
//
//     skipUntil(isWhitSspace)
//      - moves to the next whitespace
//     skipUntil(isWhiteSpace, true)
//      - moves to the next whitespace, then to the next non-whitespace
//      - is effectively (assuming no non-printable bytes):
//          skipUntil(isWhiteSpace)
//          skipUntil(isVisible)
//
//  skipLine()
//   - advances to the end of the line and over any line termination.
//   - line termination is one instance of
//       \n   - on POSIX systems
//       \r\n - on DOS systems (also most other early OSs)
//       \r   - on 1980s systems (also Lisp Machine)
//       \n\r - on nothing (makes no sense; the original use of \r\n was so
//                          that during the very slow return of the print
//                          head to the origin (\r) the paper could be
//                          advanced to the next line (\n).  \n\r is thus
//                          slower by one character transmission time).
//   - is NOT equivalent to skipUntil(isWhiteSpace, true) (above) as
//     that will skip any amount of whitespace at the end.
//
//  skipVisible()    - alias for skip(isVisible)
//  skipWhitespace() - alias for skip(isWhiteSpace)
//
//  The read variants are the same, except they copy the skipped letters
//  into the output buffer, however, readUntil(A, true) does not copy the
//  letters during the 'next' phase.
//


//  State of buffer is somewhat inconsistent;
//  should functions leave _bPos == _bLen (and so they must call fillBuffer first)
//  or should they ensure the buffer is valid (call at the end)?


namespace merylutil::inline files::inline v1 {

class readBuffer {
private:
  void          initialize(const char *pfx, char sep, const char *sfx, uint64 bMax);

public:
  readBuffer(const char *pfx, char sep, const char *sfx, uint64 bMax = 32 * 1024) { initialize(pfx, sep, sfx, bMax); }
  readBuffer(const char *filename,                       uint64 bMax = 32 * 1024) { initialize(filename, '.', nullptr, bMax); }
  readBuffer(FILE *F,                                    uint64 bMax = 32 * 1024);
  ~readBuffer() {
    delete [] _b;
    if (_owned == true)
      ::close(_f);
  };

public:
  const char   *filename(void) { return(_fn); }

public:
  void          seek(uint64 pos, uint64 extra=0);
  uint64        tell(void)   { return(_fPos); }

  bool          eof(void)    { return(_eof); }             //  True if next read will hit EOF.

  //
  //  Character support.
  //

public:
  char          peek(void) {                               //  Return the current letter.
    return (_eof == false) ? _b[_bPos] : 0;
  }
  char          next(void) {                               //  Return current then maybe move to next.
    char ch = (_eof == false) ? _b[_bPos] : 0;
    advance();
    return ch;
  }

  uint64        read(void *buf, uint64 len);               //  Copy 'len' letters into 'buf'.

  //
  //  Word support.
  //

  template<typename FN>
  char           skip(FN charclass, bool sense=true, bool next=false);

  char           skipLine(void);
  char           skipVisible(void)     { return skip(isVisible);    };
  char           skipWhitespace(void)  { return skip(isWhiteSpace); };

  template<typename N, typename FN>
  char           copy(char *s, N &sLen, N &sMax, FN charclass, bool sense=true, bool next=false);

#if 0
  template<typename N, typename FN>
  void           copyUntil(char *s, N &sLen, N &sMax, FN charclass, bool next=false);
  template<typename N>
  void           copyLine(char *s, N &sLen, N &sMax);
#endif

  template<typename N>
  void           copyVisible(char *s, N &sLen, N &sMax)     { copy(s, sLen, sMax, isVisible);   };

  //
  //  IFF support.
  //

public:
  bool          peekIFFchunk(char *name, uint32 &dataLen);
  void          readIFFchunk(char   *name,                 //  Read a generic chunk into re-allocatable memory.
                             uint8 *&data,
                             uint32 &dataLen,
                             uint32 &dataMax);
  bool          readIFFchunk(char const *name,             //  Read a specific chunk into a pre-allocated buffer.
                             void       *data,
                             uint32      dataLen);

  template<typename OBJ>
  bool          readIFFarray(char const *name, OBJ *array, uint32 nObjects) {  //  Read specific chunk
    return(readIFFchunk(name, array, sizeof(OBJ) * nObjects));                 //  into an array.
  }

  template<typename OBJ>
  bool          readIFFobject(char const *name, OBJ &object) {                 //  Read specific chunk
    return(readIFFchunk(name, &object, sizeof(OBJ)));                          //  into object via reference.
  }

  template<typename OBJ>
  bool          readIFFobject(char const *name, OBJ *object) {                 //  Read specific chunk
    return(readIFFchunk(name, object, sizeof(OBJ)));                           //  into object via pointer.
  }


private:
  void          fillBufferImpl(void);            //  Fill the buffer if it is completely empty.

  void          fillBuffer(void) {               //  Fill the buffer if it is completely empty.
    if ((_eof == false) &&                       //  The usual case is 'not empty' and we inline
        (_bPos == _bLen))                        //  that test only calling a function if
      fillBufferImpl();                          //  absolutely necessary.
  }

  void          advance(void) {                  //  Advance position.
    if (_eof == false) {
      _bPos++;
      _fPos++;
      fillBuffer();
    }
  }

  template<typename FN>
  void          advanceif(FN charclass) {        //  Advance position only if current letter is in charclass.
    if (_eof == false) {
      if (charclass(_b[_bPos]) == false)
        return;
      advance();
    }
  }


private:
  char         _fn[FILENAME_MAX+1] = {0};        //  Actual or symbolic name of opened file.

  int          _f      = -1;                     //  Handle for file.
  uint64       _fPos   = 0;                      //  Position we're at in the file.

  bool         _owned  = false;                  //  True if we opened the file.
  bool         _stdin  = false;                  //  True if reading from stdin.
  bool         _eof    = false;                  //  True if we hit EOF.

  uint64       _bBgn  = 0;                       //  File position where this buffer is from.

  uint64       _bPos  = 0;                       //  Position in the buffer we're at.
  uint64       _bLen  = 0;                       //  Length of the valid data in the buffer.
  uint64       _bMax  = 0;                       //  Size of _b allocation.
  uint8       *_b     = nullptr;                 //  Data!
};





template<typename FN>
char
readBuffer::skip(FN cc, bool sense, bool next) {

  while (_eof == false) {
    while ((_bPos < _bLen) &&          //  While stuff in the input buffer, and
           (cc(_b[_bPos]) == sense))   //  letter is in the class (sense==true)
      _bPos++, _fPos++;                //  skip it.

    if (_bPos == _bLen)                //  Get more data if we hit the end
      fillBuffer();                    //  of the buffer,
    else                               //  or stop if we found a stop letter.
      break;
  }

  if (_eof == true)    return 0;
  if (next == false)   return _b[_bPos];

  while (_eof == false) {              //  Skip again, but this opposite to
    while ((_bPos < _bLen) &&          //  what we did above.
           (cc(_b[_bPos]) != sense))
      _bPos++, _fPos++;

    if (_bPos == _bLen)
      fillBuffer();
    else
      break;
  }

  return (_eof == true) ? 0 : _b[_bPos];
}



template<typename N, typename FN>
char
readBuffer::copy(char *s, N &sLen, N &sMax, FN cc, bool sense, bool next) {

  s[sLen] = 0;

  while (_eof == false) {
    while ((sLen < sMax) &&            //  While space in the output, and
           (_bPos < _bLen) &&          //  stuff in the input buffer, and
           (cc(_b[_bPos]) == sense)) { //  letter is in the class (sense==true)
      s[sLen] = _b[_bPos];             //  copy to output.
      sLen++, _bPos++, _fPos++;        //
    }

    if (_bPos == _bLen)                //  Get more data if we hit the end
      fillBuffer();                    //  of the buffer,
    else                               //  or stop if we found a stop letter.
      break;
  }

  if (_eof == true)    return 0;
  if (next == false)   return _b[_bPos];

  while (_eof == false) {              //  Skip again, but this opposite to
    while ((_bPos < _bLen) &&          //  what we did above.
           (cc(_b[_bPos]) != sense))
      _bPos++, _fPos++;

    if (_bPos == _bLen)
      fillBuffer();
    else
      break;
  }

  return (_eof == true) ? 0 : _b[_bPos];
}



inline
char
readBuffer::skipLine(void) {
  skip(isEndOfLine, false, false);       //  Skip non-EOL letter.

  if      (_eof == true)                 //  If at EOF now,
    ;                                    //    do nothing.
  else if (isCR(_b[_bPos]))              //  If a CR,
    advance(), advanceif(isLF);          //    skip it and any single next LF.
  else if (isLF(_b[_bPos]))              //  If a LF,
    advance();                           //    skip it.
  else                                   //  Otherwise,
    assert(0);                           //    busted, expecting CR or LF!

  return (_eof == true) ? 0 : _b[_bPos];
}

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_FILES_BUFFERED_V1_READING_H
