
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

#include "dnaSeqFile-v1.H"

#include "arrays.H"
#include "strings.H"

//
//  Utility function for finding index files.
//
const uint64 dnaSeqVersion01 = 0x3130716553616e64;   //  dnaSeq01
const uint64 dnaSeqVersion02 = 0x3230716553616e64;   //  dnaSeq02 - not used yet

static
char const *
makeIndexName(char const *prefix) {
  char const *suffix = ".dnaSeqIndex";
  uint32      plen   = strlen(prefix);
  uint32      slen   = strlen(suffix);
  char       *iname  = new char [plen + slen + 1];

  memcpy(iname,        prefix, plen + 1);   //  +1 for the NUL byte.
  memcpy(iname + plen, suffix, slen + 1);

  return iname;
}



namespace merylutil::inline sequence::inline v1 {

//
//  File management.
//

bool
bufSeqFile::open(char const *fn, bool indexed) {
  dnaSeqFile::filename(fn);

  _file       = new compressedFileReader(filename());
  _buffer     = new readBuffer(_file->file(), 128 * 1024);
  _indexable  = false;
  _reopenable = false;
  _compressed = false;

  _buffer->skipWhitespace();

  if ((_buffer->eof() == false) &&
      (_buffer->peek() != '>') &&          //  Close file if
      (_buffer->peek() != '@'))            //  not FASTA/Q.
    return close(); 

  //  Decide FASTQ or SAM, heuristically.
  //   - Detail on the seek(0) below: if the input is from a pipe and
  //     the buffer needs to refresh itself, we cannot seek back to
  //     the start.  To prevent this, make sure the buffer size above
  //     is more than the 'first N letters' limnit below.

  if (_buffer->peek() == '@') {
    _buffer->skipLine();                   //  Skip rest of first line,
    _buffer->skipWhitespace();             //  and any stray whitespace.

    if (_buffer->peek() == '@')            //  If second line is another '@',
      return close();                      //  we're definitely in a SAM header.

    while ((_buffer->peek() != '\n') &&    //  Check for tabs in second line; if
           (_buffer->tell() < 1000))       //  found, we're _likely_ a SAM.  But if
      if (_buffer->next() == '\t')         //  not found in first 1000 letters, either
        return close();                    //  fastq or a very long read name.

    _buffer->skipWhitespace();             //  Skip any stray whitespace.

    if ((_buffer->tell() == 1000) ||       //  If no tabs in first 1000 letters or
        (_buffer->peek() == '+'))          //  third line is a '+', we're probably
      _buffer->seek(0);                    //  FASTQ; go back to the start of the file.
    else                                   //
      return close();                      //  Otherwise, give up and declare it SAM.
  }

  _indexable    = _file->isSeekable();
  _reopenable   = _file->isReopenable();
  _compressed   = _file->isCompressed();

  loadIndex(indexed);

  return true;
}


bool
bufSeqFile::close(void) {
  if (_file) {
    delete _file;    _file   = nullptr;
    delete _buffer;  _buffer = nullptr;
    unloadIndex();
  }
  return false;
}


bool
bufSeqFile::reopen(void) {
  close();
  if (_reopenable == false)
    return false;
  return open(filename());
}

//
//  Index loading/saving/creation
//

bool
bufSeqFile::loadIndex(bool create) {
  char const  *indexName = makeIndexName(filename());
  FILE        *indexFile = nullptr;

  if (fileExists(indexName) == true) {
    FILE   *indexFile = merylutil::openInputFile(indexName);
    uint64  magic;
    uint64  size;
    uint64  date;

    loadFromFile(magic,     "bufSeqFile::magic",    indexFile);
    loadFromFile(size,      "bufSeqFile::size",     indexFile);
    loadFromFile(date,      "bufSeqFile::date",     indexFile);
    loadFromFile(_indexLen, "bufSeqFile::indexLen", indexFile);

    if (magic != dnaSeqVersion01) {
      fprintf(stderr, "ERROR: file '%s' isn't a dnaSeqIndex; manually remove this file.\n", indexName);
      exit(1);
    }

    if ((size == merylutil::sizeOfFile(filename())) &&
        (date == merylutil::timeOfFile(filename()))) {
      _indexMax = _indexLen;
      _index    = new dsfIdxEnt [_indexMax];

      loadFromFile(_index, "bufSeqFile::index", _indexLen, indexFile);
    }
    else {
      fprintf(stderr, "WARNING: file '%s' disagrees with index; recreating index.\n", filename());

      _index    = nullptr;
      _indexLen = 0;
      _indexMax = 0;
    }

    merylutil::closeFile(indexFile, indexName);
  }

  delete [] indexName;

  if ((_index == nullptr) && (create == true)) {
    createIndex();
    loadIndex(false);
  }

  _numSeqs = _indexLen;

  return _numSeqs > 0;    //  Return true if we have an index.
}


void
bufSeqFile::unloadIndex(void) {
  delete [] _index;
  _indexLen = 0;
  _indexMax = 0;
  _index    = nullptr;
}


void
bufSeqFile::createIndex(void) {
  dnaSeq     seq;

  //  Fail if an index is requested for a compressed file.

  if (_file->isCompressed() == true)
    fprintf(stderr, "ERROR: cannot index compressed input '%s'.\n", filename()), exit(1);

  if (_file->isNormal() == false)
    fprintf(stderr, "ERROR: cannot index pipe input.\n"), exit(1);

  //  Rewind the buffer to make sure we're at the start of the file.

  _buffer->seek(0);

  //  Allocate space for the index, set the first entry to the current
  //  position of the file.

  _indexLen = 0;
  _indexMax = 1048576;
  _index    = new dsfIdxEnt [_indexMax];

  _index[0]._fileOffset     = _buffer->tell();
  _index[0]._sequenceLength = 0;

  //  While we read sequences:
  //    update the length of the sequence (we've already saved the position)
  //    make space for more sequences
  //    save the position of the next sequence

  while (dnaSeqFile::loadSequence(seq) == true) {
    if (seq.wasError()) {
      fprintf(stderr, "WARNING: error reading sequence at/before '%s'\n", seq.ident());
    }

    if (seq.wasReSync()) {
      fprintf(stderr, "WARNING: lost sync reading before sequence '%s'\n", seq.ident());
    }

    _index[_indexLen]._sequenceLength = seq.length();

    increaseArray(_index, _indexLen, _indexMax, 1048576);

    _indexLen++;

    _index[_indexLen]._fileOffset     = _buffer->tell();
    _index[_indexLen]._sequenceLength = 0;
  }

  //  Save whatever index we made.

  char const *indexName = makeIndexName(filename());
  FILE       *indexFile = merylutil::openOutputFile(indexName);

  uint64  magic = dnaSeqVersion01;
  uint64  size  = merylutil::sizeOfFile(filename());
  uint64  date  = merylutil::timeOfFile(filename());

  writeToFile(magic,     "bufSeqFile::magic",    indexFile);
  writeToFile(size,      "bufSeqFile::size",     indexFile);
  writeToFile(date,      "bufSeqFile::date",     indexFile);
  writeToFile(_indexLen, "bufSeqFile::indexLen", indexFile);
  writeToFile(_index,    "bufSeqFile::index",    _indexLen, indexFile);

  merylutil::closeFile(indexFile, indexName);

  delete [] indexName;
}



void
bufSeqFile::destroyIndex(void) {
  merylutil::unlink(makeIndexName(filename()));
}

//
//  Sequence accessors
//

bool                                    //  Position file at start of sequence 'i' and
bufSeqFile::findSequence(uint64 i) {    //  remember the index next to be loaded.
  if (_indexLen == 0)
    return false;                       //  No index or index i not in this file.
  _buffer->seek(_index[_seqIdx = i]._fileOffset);
  return true;
}

uint64
bufSeqFile::sequenceLength(uint64 i) {
  if (_indexLen <= i)
    return uint64max;                   //  No index or index i not in this file.
  return _index[i]._sequenceLength;
}

//
//
//

bool
bufSeqFile::loadFASTA(char  *&name, uint32 &nameMax,
                      char  *&seq,
                      uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint64 &qltLen) {
  uint64  nameLen = 0;
  char    ch      = _buffer->next();

  //  Skip any whitespace.

  while (isWhiteSpace(ch))
    ch = _buffer->next();

  //  Fail rather ungracefully if we aren't at a sequence start.

  if (ch != '>')
    return false;

  //  Read the header line into the name string.  We cannot skip whitespace
  //  here, but we do allow DOS to insert a \r before any \n.

  for (ch=_buffer->next(); (ch != '\n') && (ch != 0); ch=_buffer->next()) {
    if (ch == '\r')
      continue;
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }
  name[nameLen] = 0;

  //  Read sequence, skipping whitespace, until we hit a new sequence or eof.

  seqLen = 0;
  qltLen = 0;

  for (ch = _buffer->peek(); ((ch != '>') &&
                              (ch != '@') &&
                              (ch !=  0)); ch = _buffer->peek()) {
    assert(_buffer->eof() == false);

    ch = _buffer->next();

    if (isWhiteSpace(ch))
      continue;

    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);

    seq[seqLen++] = ch;
    qlt[qltLen++] = 0;
  }

  seq[seqLen] = 0;
  qlt[qltLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);
  assert(qltLen  < seqMax);

  _seqIdx++;

  return true;
}



bool
bufSeqFile::loadFASTQ(char  *&name, uint32 &nameMax,
                      char  *&seq,
                      uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint64 &qltLen) {
  uint32  nameLen = 0;
  char    ch      = _buffer->next();

  //  Skip any whitespace.

  while (isWhiteSpace(ch))
    ch = _buffer->next();

  //  Fail rather ungracefully if we aren't at a sequence start.

  if (ch != '@')
    return false;

  //  Read the header line into the name string.  We cannot skip whitespace
  //  here, but we do allow DOS to insert a \r before any \n.

  for (ch=_buffer->next(); (ch != '\n') && (ch != 0); ch=_buffer->next()) {
    if (ch == '\r')
      continue;
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Trim back the header line to remove white space at the end.

  while ((nameLen > 0) && (isWhiteSpace(name[nameLen-1])))
    nameLen--;

  name[nameLen] = 0;

  //  Skip any whitespace, again.  Once we hit non-whitespace we'll suck in
  //  the whole line.

  while (isWhiteSpace(ch))
    ch = _buffer->next();

  //  Read sequence.  Pesky DOS files end with \r\n, and it suffices
  //  to stop on the \n and ignore all the rest.

  seqLen = 0;
  qltLen = 0;

  for (; (ch != '\n') && (ch != 0); ch=_buffer->next()) {
    if (isWhiteSpace(ch))
      continue;
    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);
    seq[seqLen++] = ch;
  }

  //  Skip any more whitespace, fail if we're not at a quality start, then
  //  suck in the quality line.  And then skip more whitespace.

  while (isWhiteSpace(ch))
    ch = _buffer->next();

  if (ch != '+')
    return false;

  for (ch=_buffer->next(); (ch != '\n') && (ch != 0); ch=_buffer->next()) {
    ;
  }

  while (isWhiteSpace(ch))
    ch = _buffer->next();

  //  Read qualities and convert to integers.

  for (; (ch != '\n') && (ch != 0); ch=_buffer->next()) {
    if (isWhiteSpace(ch))
      continue;
    if (qltLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, qltLen, seqMax, 3 * seqMax / 2);
    qlt[qltLen++] = ch - '!';
  }

  //  Skip whitespace after the sequence.  This one is a little weird.  It
  //  tests if the _next_ letter is whitespace, and if so, gets it from the
  //  buffer.  After this loop, the _next_ letter in the buffer should be
  //  either a '>' or a '@'.

  while (isWhiteSpace(_buffer->peek()))
    _buffer->next();

  seq[seqLen] = 0;
  qlt[qltLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);
  assert(qltLen  < seqMax);

  _seqIdx++;

  return true;
}



bool
bufSeqFile::loadSequence(char  *&name, uint32 &nameMax,
                         char  *&seq,
                         uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint32 &error) {
  uint64 qltLen = 0;

  //  Allocate space for the arrays, if they're currently unallocated.

  if (nameMax == 0)
    resizeArray(name, 0, nameMax, (uint32)1024);

  if (seqMax == 0)
    resizeArrayPair(seq, qlt, 0, seqMax, (uint64)65536);

  //  Clear our return values.

  bool   loadSuccess = false;

  _isFASTA = false;
  _isFASTQ = false;

  name[0] = 0;
  seq[0]  = 0;
  qlt[0]  = 0;
  seqLen  = 0;

  error   = 0;

  //  Skip any whitespace at the start of the file, or before the next FASTQ
  //  sequence (the FASTA reader will automagically skip whitespace at the
  //  end of the sequence).  (minor duplication in the constructor)

  while (isWhiteSpace(_buffer->peek()))
    _buffer->next();

  //  If we're not at a sequence start, scan ahead to find the next one.
  //  Not bulletproof; FASTQ qv's can match this.

  if ((_buffer->peek() != '>') &&
      (_buffer->peek() != '@') &&
      (_buffer->peek() !=  0)) {
    //fprintf(stderr, "bufSeqFile::loadSequence()-- sequence sync lost at position %lu, attempting to find the next sequence.\n", _buffer->tell());
    error |= 0x02;
  }

  bool  lastWhite = isWhiteSpace(_buffer->peek());

  while ((_buffer->peek() != '>') &&
         (_buffer->peek() != '@') &&
         (_buffer->peek() !=  0)) {
    _buffer->next();
  }

  //  Peek at the file to decide what type of sequence we need to read.

  if      (_buffer->peek() == '>') {
    _isFASTA    = true;
    loadSuccess = loadFASTA(name, nameMax, seq, qlt, seqMax, seqLen, qltLen);
  }

  else if (_buffer->peek() == '@') {
    _isFASTQ    = true;
    loadSuccess = loadFASTQ(name, nameMax, seq, qlt, seqMax, seqLen, qltLen);
  }

  else {
    _isFASTA = false;
    _isFASTQ = false;

    return false;
  }

  //  If we failed to load a sequence, report an error message and zero out
  //  the sequence.  Leave the name as-is so we can at least return a length
  //  zero sequence.  If we failed to load a name, it'll still be set to NUL.

  if (loadSuccess == false) {
    //if (name[0] == 0)
    //  fprintf(stderr, "bufSeqFile::loadSequence()-- failed to read sequence correctly at position %lu.\n", _buffer->tell());
    //else
    //  fprintf(stderr, "bufSeqFile::loadSequence()-- failed to read sequence '%s' correctly at position %lu.\n", name, _buffer->tell());

    error |= 0x01;

    seq[0]  = 0;
    qlt[0]  = 0;
    seqLen  = 0;
  }

  return true;
}



bool
bufSeqFile::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength,
                      bool    &endOfSequence) {

  endOfSequence = true;                  //  Assume we'll read the whole sequence.

  _buffer->skipWhitespace();             //  Advance to the first real character.

  if (_buffer->eof() == true)            //  If at EOF, not much to do.
    return false;

  if (_buffer->peek() == '+') {          //  If unlucky enough to be at a QV line:
    _buffer->skipLine();                 //    Skip the + line.
    _buffer->skipLine();                 //    Skip the QV line.
  }

  if ((_buffer->peek() == '>') ||        //  Usually, we're at a FASTA or
      (_buffer->peek() == '@')) {        //  FASTQ name line:
    _buffer->skipLine();                 //    Skip the name line.
  }

  while (_buffer->eof() == false) {
    _buffer->skipWhitespace();
    _buffer->copyVisible(seq, seqLength, maxLength);
    _buffer->skipWhitespace();

    assert(seqLength <= maxLength);

    if ((_buffer->peek() == '>') ||      //  If at the next sequence,
        (_buffer->peek() == '@'))        //  stop reading and return success.
      return true;

    if (seqLength == maxLength) {        //  But if out of space, note that
      endOfSequence = false;             //  we're in the middle of a sequence,
      return true;                       //  and return success.
    }

    _buffer->next();
  }

  return seqLength > 0;                  //  EOF.  Return true if bases exist.
}



}  //  namespace merylutil::sequence::v1
