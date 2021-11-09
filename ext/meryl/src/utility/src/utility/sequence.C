
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

#include "sequence.H"
#include "arrays.H"



static
const
char
inv[256] = {
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x00 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x08 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x10 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x18 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x20 -  !"#$%&'
   0,  0,  0,  0,  0, '-', 0,  0,  //  0x28 - ()*+,-./
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x30 - 01234567
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x38 - 89:;<=>?
   0, 'T', 0, 'G', 0,  0,  0, 'C', //  0x40 - @ABCDEFG
   0,  0,  0,  0,  0,  0, 'N', 0,  //  0x48 - HIJKLMNO
   0,  0,  0,  0, 'A', 0,  0,  0,  //  0x50 - PQRSTUVW
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x58 - XYZ[\]^_
   0, 't', 0, 'g', 0,  0,  0, 'c', //  0x60 - `abcdefg
   0,  0,  0,  0,  0,  0, 'n', 0,  //  0x68 - hijklmno
   0,  0,  0,  0, 'a', 0,  0,  0,  //  0x70 - pqrstuvw
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x78 - xyz{|}~
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x80 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x88 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x90 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0x98 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xa0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xa8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xb0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xb8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xc0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xc8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xd0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xd8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xe0 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xe8 -
   0,  0,  0,  0,  0,  0,  0,  0,  //  0xf0 -
   0,  0,  0,  0,  0,  0,  0,  0   //  0xf8 -
};

static
const
char
Dacgtn[5] = { 'A',
              'C',
              'G',
              'T',
              'N' };

static
const
uint8
Eacgtn[256] = {
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x00 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x08 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x10 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x18 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x20 -  !"#$%&'
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x28 - ()*+,-./
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x30 - 01234567
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x38 - 89:;<=>?
   0, 0x00,    0, 0x01,    0,    0,    0, 0x02,    //  0x40 - @ABCDEFG
   0,    0,    0,    0,    0,    0, 0x04,    0,    //  0x48 - HIJKLMNO
   0,    0,    0,    0, 0x03,    0,    0,    0,    //  0x50 - PQRSTUVW
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x58 - XYZ[\]^_
   0, 0x00,    0, 0x01,    0,    0,    0, 0x02,    //  0x60 - `abcdefg
   0,    0,    0,    0,    0,    0, 0x04,    0,    //  0x68 - hijklmno
   0,    0,    0,    0, 0x03,    0,    0,    0,    //  0x70 - pqrstuvw
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x78 - xyz{|}~
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x80 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x88 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x90 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0x98 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xa0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xa8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xb0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xb8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xc0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xc8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xd0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xd8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xe0 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xe8 -
   0,    0,    0,    0,    0,    0,    0,    0,    //  0xf0 -
   0,    0,    0,    0,    0,    0,    0,    0     //  0xf8 -
};




void
reverseComplementSequence(char *seq, int len) {
  char   c=0;
  char  *s=seq,  *S=seq+len-1;

  if (len == 0) {
    len = strlen(seq);
    S = seq + len - 1;
  }

  while (s < S) {
    c    = *s;
    *s++ =  inv[*S];
    *S-- =  inv[c];
  }

  if (s == S)
    *s = inv[*s];
}



char *
reverseComplementCopy(char *seq, int len) {
  char  *rev = new char [len+1];

  assert(len > 0);

  for (int32 p=len, q=0; p>0; )
    rev[q++] = inv[seq[--p]];

  rev[len] = 0;

  return(rev);
}



template<typename qvType>
void
reverseComplement(char *seq, qvType *qlt, int len) {
  char    c=0;
  char   *s=seq,  *S=seq+len-1;
  qvType *q=qlt,  *Q=qlt+len-1;

  if (qlt == NULL) {
    reverseComplementSequence(seq, len);
    return;
  }

  assert(len > 0);

  if (len == 0) {
    len = strlen(seq);
    S = seq + len - 1;
    Q = qlt + len - 1;
  }

  while (s < S) {
    c    = *s;
    *s++ =  inv[*S];
    *S-- =  inv[c];

    c    = *q;
    *q++ = *Q;
    *Q-- =  c;
  }

  if (s == S)
    *s = inv[*s];
}

template void reverseComplement<char> (char *seq, char  *qlt, int len);   //  Give the linker
template void reverseComplement<uint8>(char *seq, uint8 *qlt, int len);   //  something to link



//  Compress homopolymer runs to a single letter.  Returns the length of the
//  compressed sequence.
//
//  'bases' does not need to be NUL terminated.
//
//  If 'compr' is supplied, the compressed sequence is returned.  'compr' can be
//  the same array as 'bases'.  The output string is always NUL terminated.
//
//  If 'ntoc' is supplied, a mapping from 'normal' position to 'compressed'
//  position is returned.  This output is INCORRECT if 'skip' is set (the
//  values of any skipped bases are undefined).
//
//  If 'skip' is set, ignore these bases at the start of the sequence that
//  are the same.  This is to allow one to homopoly compress a sequence in
//  chunks; compress the first 1000 bases, the the next 1000, and so on.
//  Each pass sets skip to the last base of the previous chunk.
//
uint32
homopolyCompress(char *bases, uint32 basesLen, char *compr, uint32 *ntoc, char skip) {
  uint32  cc = 0;  //  position of the start of the run
  uint32  rr = 1;  //  position of the scan head
  uint32  sl = 0;  //  length of the compressed sequence

  while ((bases[cc] == skip) &&   //  If 'skip' is set, ignore these bases
         (cc < basesLen)) {        //  at the start of 'bases'.
    cc++;
    rr++;
  }

  if (compr)                      //  Either the first base, or
    compr[sl] = bases[cc];        //  the terminating NUL.

  if (ntoc)                       //  Save the mapping from the first
    ntoc[cc] = sl;                //  normal to first compressed base.

  if (basesLen == 0)
    return(0);

  sl++;

  while (rr < basesLen) {

    //  In a run, move the scan head one position, and set the
    //  mapping to the last compressed base.
    if ((bases[cc] | 0x20) == (bases[rr] | 0x20)) {    //  Convert to lowercase before comparing.
      if (ntoc)
        ntoc[rr] = sl - 1;
      rr++;
    }

    //  Just ended a run.  Set the start of the (next) run
    //  to the current position, move the current position
    //  to the next base, and increase the length of the
    //  compressed sequence.
    else {
      if (compr)
        compr[sl] = bases[rr];
      if (ntoc)
        ntoc[rr] = sl;
      cc = rr;
      rr++;
      sl++;
    }
  }

  //  Terminate the compressed string.
  if (compr)
    compr[sl] = 0;

  //  The 'space' after the end of the bases maps to the 'space'
  //  after the compressed bases.
  if (ntoc)
    ntoc[rr]  = sl;

  return(sl);
}



void
decode2bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {
  uint32       chunkPos = 0;

  assert(seq != NULL);

  for (uint32 ii=0; ii<seqLen; ) {
    if (chunkPos == chunkLen) {
      fprintf(stderr, "decode2bit()-- ran out of chunk (length %u) before end of sequence (at %u out of %u)\n",
              chunkLen, ii, seqLen);
    }
    assert(chunkPos < chunkLen);

    uint8  byte = chunk[chunkPos++];

    if (ii + 4 < seqLen) {
      seq[ii++] = Dacgtn[((byte >> 6) & 0x03)];
      seq[ii++] = Dacgtn[((byte >> 4) & 0x03)];
      seq[ii++] = Dacgtn[((byte >> 2) & 0x03)];
      seq[ii++] = Dacgtn[((byte >> 0) & 0x03)];
    }

    else {
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 6) & 0x03)];  //  This if is redundant,
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 4) & 0x03)];  //    but pretty.
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 2) & 0x03)];
      if (ii < seqLen)  seq[ii++] = Dacgtn[((byte >> 0) & 0x03)];
    }
  }

  seq[seqLen] = 0;
}



uint32
encode2bitSequence(uint8 *&chunk, char *seq, uint32 seqLen) {

  for (uint32 ii=0; ii<seqLen; ii++) {       //  If non-ACGT present, return
    char  base = seq[ii];                    //  0 to indicate we can't encode.

    if ((base != 'a') && (base != 'A') &&
        (base != 'c') && (base != 'C') &&
        (base != 'g') && (base != 'G') &&
        (base != 't') && (base != 'T')) {
      fprintf(stderr, "Invalid base %c detected at position %u\n", base, ii);
      return(0);
    }
  }

  uint32 chunkLen = 0;

  if (chunk == NULL)
    chunk = new uint8 [ seqLen / 4 + 1];

  for (uint32 ii=0; ii<seqLen; ) {
    uint8  byte = 0;

    if (ii + 4 < seqLen) {
      byte  = Eacgtn[seq[ii++]];  byte <<= 2;
      byte |= Eacgtn[seq[ii++]];  byte <<= 2;
      byte |= Eacgtn[seq[ii++]];  byte <<= 2;
      byte |= Eacgtn[seq[ii++]];
    }

    else {
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }   byte <<= 2;  //  This if is redundant, but pretty.
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }   byte <<= 2;  //  Yes, all three always shift,
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }   byte <<= 2;  //  not conditionally.
      if (ii < seqLen)  { byte |= Eacgtn[seq[ii++]]; }
    }

    chunk[chunkLen++] = byte;
  }

  return(chunkLen);
}



void
decode3bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {
  uint32       chunkPos = 0;

  assert(seq != NULL);

  for (uint32 ii=0; ii<seqLen; ) {
    if (chunkPos == chunkLen) {
      fprintf(stderr, "decode3bit()-- ran out of chunk (length %u) before end of sequence (at %u out of %u)\n",
              chunkLen, ii, seqLen);
    }
    assert(chunkPos < chunkLen);

    uint8  byte = chunk[chunkPos++];
    uint8  c1   = 0;
    uint8  c2   = 0;
    uint8  c3   = 0;

    if (ii + 3 < seqLen) {
      uint8 c1 = byte / 5 / 5;    byte -= c1 * 5 * 5;
      uint8 c2 = byte / 5;        byte -= c2 * 5;
      uint8 c3 = byte;

      seq[ii++] = Dacgtn[c1];
      seq[ii++] = Dacgtn[c2];
      seq[ii++] = Dacgtn[c3];
    }

    else {
      uint8 c1 = byte / 5 / 5;    byte -= c1 * 5 * 5;
      uint8 c2 = byte / 5;        byte -= c2 * 5;
      uint8 c3 = byte;

      if (ii < seqLen)  seq[ii++] = Dacgtn[c1];  //  This if is redundant,
      if (ii < seqLen)  seq[ii++] = Dacgtn[c2];  //    but pretty.
      if (ii < seqLen)  seq[ii++] = Dacgtn[c3];
    }
  }

  seq[seqLen] = 0;
}



uint32
encode3bitSequence(uint8 *&chunk, char *seq, uint32 seqLen) {

  for (uint32 ii=0; ii<seqLen; ii++) {       //  If non-ACGTN present, return
    char  base = seq[ii];                    //  0 to indicate we can't encode.

    if ((base != 'a') && (base != 'A') &&
        (base != 'c') && (base != 'C') &&
        (base != 'g') && (base != 'G') &&
        (base != 't') && (base != 'T') &&
        (base != 'n') && (base != 'N')) {
      fprintf(stderr, "Invalid base %c detected at position %u\n", base, ii);
      return(0);
    }
  }

  uint32 chunkLen = 0;

  if (chunk == NULL)
    chunk = new uint8 [ seqLen / 3 + 1];

  for (uint32 ii=0; ii<seqLen; ) {
    uint8  byte = 0;

    if (ii + 3 < seqLen) {
      byte += Eacgtn[seq[ii++]] * 5 * 5;
      byte += Eacgtn[seq[ii++]] * 5;
      byte += Eacgtn[seq[ii++]];
    }

    else {
      if (ii < seqLen)   byte += Eacgtn[seq[ii++]] * 5 * 5;
      if (ii < seqLen)   byte += Eacgtn[seq[ii++]] * 5;
      if (ii < seqLen)   byte += Eacgtn[seq[ii++]];
    }

    chunk[chunkLen++] = byte;
  }

  return(chunkLen);
}



void
decode8bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {

  assert(seq != NULL);

  for (uint32 ii=0; ii<seqLen; ii++)
    seq[ii] = chunk[ii];

  seq[seqLen] = 0;
}



uint32
encode8bitSequence(uint8 *&chunk, char *seq, uint32 seqLen) {

  if (chunk == NULL)
    chunk = new uint8 [ seqLen ];

  for (uint32 ii=0; ii<seqLen; ii++)
    chunk[ii] = seq[ii];

  return(seqLen);
}



////////////////////////////////////////
//  dnaSeq functions
//

dnaSeq::dnaSeq() {
};


dnaSeq::~dnaSeq() {
  delete [] _name;
  delete [] _seq;
  delete [] _qlt;
};


void
dnaSeq::releaseAll(void) {
  delete [] _name;    _name = _ident = _flags = nullptr;
  delete [] _seq;     _seq                    = nullptr;
  delete [] _qlt;     _qlt                    = nullptr;

  _nameMax = 0;
  _seqMax  = 0;
  _seqLen  = 0;
}


void
dnaSeq::releaseBases(void) {
  delete [] _seq;     _seq                    = nullptr;
  delete [] _qlt;     _qlt                    = nullptr;

  _seqMax  = 0;
  _seqLen  = 0;
}


bool
dnaSeq::copy(char  *bout,
             uint32 bgn, uint32 end, bool terminate) {

  if ((end < bgn) || (_seqLen < end))
    return(false);

  for (uint32 ii=bgn; ii<end; ii++)
    bout[ii-bgn] = _seq[ii];

  if (terminate)
    bout[end-bgn] = 0;

  return(true);
}


bool
dnaSeq::copy(char  *bout,
             uint8 *qout,
             uint32 bgn, uint32 end, bool terminate) {

  if ((end < bgn) || (_seqLen < end))
    return(false);

  for (uint32 ii=bgn; ii<end; ii++) {
    bout[ii-bgn] = _seq[ii];
    qout[ii-bgn] = _qlt[ii];
  }

  if (terminate) {
    bout[end-bgn] = 0;
    qout[end-bgn] = 0;
  }

  return(true);
}


void
dnaSeq::findNameAndFlags(void) {
  uint32 ii=0;

  while (isWhiteSpace(_name[ii]) == true)   //  Skip white space before the name.
    ii++;                                   //  Why do you torture us?

  _ident = _name + ii;                      //  At the start of the name.

  while (isVisible(_name[ii]) == true)      //  Skip over the name.
    ii++;

  if (isNUL(_name[ii]) == true) {           //  If at the end of the string,
    _flags = _name + ii;                    //  there are no flags,
    return;                                 //  so just return.
  }

  _name[ii++] = 0;                          //  Terminate the name, move ahead.

  while (isWhiteSpace(_name[ii]) == true)   //  Otherwise, skip whitespace
    ii++;                                   //  to get to the flags.

  _flags = _name + ii;                      //  Flags are here or NUL.
}



////////////////////////////////////////
//  dnaSeqFile functions
//

dnaSeqFile::dnaSeqFile(char const *filename, bool indexed) {
  _filename = duplicateString(filename);

  reopen(indexed);
}



dnaSeqFile::~dnaSeqFile() {
  delete [] _filename;
  delete    _file;
  delete    _buffer;
  delete [] _index;
}



//  Open, or reopen, an input file.
//
void
dnaSeqFile::reopen(bool indexed) {

  //  If a _file exists already, reopen it, otherwise, make a new one.
  if (_file)
    _file->reopen();
  else
    _file = new compressedFileReader(_filename);

  //  Since the file object is always new, we need to make a new read buffer.
  //  gzip inputs seem to be (on FreeBSD) returning only 64k blocks
  //  regardless of the size of our buffer; but uncompressed inputs will
  //  benefit slightly from a bit larger buffer.
  delete _buffer;

  _buffer = new readBuffer(_file->file(), 128 * 1024);

  //  If we have an index already or one is requested, (re)generate it.

  if ((_index != nullptr) || (indexed == true))
    generateIndex();
}



bool
dnaSeqFile::findSequence(uint64 i) {

  if (_indexLen == 0)   return(false);
  if (_indexLen <= i)   return(false);

  _buffer->seek(_index[i]._fileOffset);

  _seqIdx = i;

  return(true);
}



uint64
dnaSeqFile::sequenceLength(uint64 i) {

  if (_indexLen == 0)   return(UINT64_MAX);
  if (_indexLen <= i)   return(UINT64_MAX);

  return(_index[i]._sequenceLength);
}




////////////////////////////////////////
//  dnaSeqFile indexing
//

const uint64 dnaSeqVersion01 = 0x3130716553616e64;   //  dnaSeq01
const uint64 dnaSeqVersion02 = 0x3230716553616e64;   //  dnaSeq02 - not used yet


char const *
makeIndexName(char const *prefix) {
  char const *suffix = ".dnaSeqIndex";
  uint32      plen   = strlen(prefix);
  uint32      slen   = strlen(suffix);
  char       *iname  = new char [plen + slen + 1];

  memcpy(iname,        prefix, plen + 1);   //  +1 for the NUL byte.
  memcpy(iname + plen, suffix, slen + 1);

  return(iname);
}


//  Load an index.  Returns true if one was loaded.
bool
dnaSeqFile::loadIndex(void) {
  char const  *indexName = makeIndexName(_filename);
  FILE        *indexFile = nullptr;

  if (fileExists(indexName) == true) {
    FILE   *indexFile = AS_UTL_openInputFile(indexName);
    uint64  magic;
    uint64  size;
    uint64  date;

    loadFromFile(magic,     "dnaSeqFile::magic",    indexFile);
    loadFromFile(size,      "dnaSeqFile::size",     indexFile);
    loadFromFile(date,      "dnaSeqFile::date",     indexFile);
    loadFromFile(_indexLen, "dnaSeqFile::indexLen", indexFile);

    if (magic != dnaSeqVersion01) {
      fprintf(stderr, "ERROR: file '%s' isn't a dnaSeqIndex; manually remove this file.\n", indexName);
      exit(1);
    }

    if ((size == AS_UTL_sizeOfFile(_filename)) &&
        (date == AS_UTL_timeOfFile(_filename))) {
      _index = new dnaSeqIndexEntry [_indexLen];

      loadFromFile(_index, "dnaSeqFile::index", _indexLen, indexFile);

    } else {
      fprintf(stderr, "WARNING: file '%s' disagrees with index; recreating index.\n", _filename);

      _index    = nullptr;
      _indexLen = 0;
      _indexMax = 0;
    }

    AS_UTL_closeFile(indexFile, indexName);
  }

  delete [] indexName;

  return(_index != nullptr);   //  Return true if we have an index.
}



void
dnaSeqFile::saveIndex(void) {
  char const *indexName = makeIndexName(_filename);
  FILE       *indexFile = AS_UTL_openOutputFile(indexName);

  uint64  magic = dnaSeqVersion01;
  uint64  size  = AS_UTL_sizeOfFile(_filename);
  uint64  date  = AS_UTL_timeOfFile(_filename);

  writeToFile(magic,     "dnaSeqFile::magic",    indexFile);
  writeToFile(size,      "dnaSeqFile::size",     indexFile);
  writeToFile(date,      "dnaSeqFile::date",     indexFile);
  writeToFile(_indexLen, "dnaSeqFile::indexLen", indexFile);
  writeToFile(_index,    "dnaSeqFile::index",    _indexLen, indexFile);

  AS_UTL_closeFile(indexFile, indexName);

  delete [] indexName;
}



void
dnaSeqFile::generateIndex(void) {
  dnaSeq     seq;

  //  Fail if an index is requested for a compressed file.

  if (_file->isCompressed() == true)
    fprintf(stderr, "ERROR: cannot index compressed input '%s'.\n", _filename), exit(1);

  if (_file->isNormal() == false)
    fprintf(stderr, "ERROR: cannot index pipe input.\n"), exit(1);

  //  If we can load an index, do it and return.

  if (loadIndex() == true)
    return;

  //  Rewind the buffer to make sure we're at the start of the file.

  _buffer->seek(0);

  //  Allocate space for the index, set the first entry to the current
  //  position of the file.

  _indexLen = 0;
  _indexMax = 1048576;
  _index    = new dnaSeqIndexEntry [_indexMax];

  _index[0]._fileOffset     = _buffer->tell();
  _index[0]._sequenceLength = 0;

  //  While we read sequences:
  //    update the length of the sequence (we've already saved the position)
  //    make space for more sequences
  //    save the position of the next sequence

  while (loadSequence(seq) == true) {
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

  saveIndex();
}



void
dnaSeqFile::removeIndex(void) {

  delete [] _index;

  _indexLen = 0;
  _indexMax = 0;
  _index    = nullptr;
}



bool
dnaSeqFile::loadFASTA(char  *&name, uint32 &nameMax,
                      char  *&seq,
                      uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint64 &qltLen) {
  uint64  nameLen = 0;
  char    ch      = _buffer->read();

  //  Skip any whitespace.

  while (isWhiteSpace(ch))
    ch = _buffer->read();

  //  Fail rather ungracefully if we aren't at a sequence start.

  if (ch != '>')
    return(false);

  //  Read the header line into the name string.  We cannot skip whitespace
  //  here, but we do allow DOS to insert a \r before any \n.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (ch == '\r')
      continue;
    if (nameLen+1 >= nameMax)
      resizeArray(name, nameLen, nameMax, 3 * nameMax / 2);
    name[nameLen++] = ch;
  }

  //  Trim back the header line to remove white space at the end.  The
  //  terminating NUL is tacked on at the end.

  while ((nameLen > 0) && (isWhiteSpace(name[nameLen-1])))
    nameLen--;

  name[nameLen] = 0;

  //  Read sequence, skipping whitespace, until we hit a new sequence or eof.

  seqLen = 0;
  qltLen = 0;

  for (ch = _buffer->peek(); ((ch != '>') &&
                              (ch != '@') &&
                              (ch !=  0)); ch = _buffer->peek()) {
    assert(_buffer->eof() == false);

    ch = _buffer->read();

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

  return(true);
}



bool
dnaSeqFile::loadFASTQ(char  *&name, uint32 &nameMax,
                      char  *&seq,
                      uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint64 &qltLen) {
  uint32  nameLen = 0;
  char    ch      = _buffer->read();

  //  Skip any whitespace.

  while (isWhiteSpace(ch))
    ch = _buffer->read();

  //  Fail rather ungracefully if we aren't at a sequence start.

  if (ch != '@')
    return(false);

  //  Read the header line into the name string.  We cannot skip whitespace
  //  here, but we do allow DOS to insert a \r before any \n.

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
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
    ch = _buffer->read();

  //  Read sequence.  Pesky DOS files end with \r\n, and it suffices
  //  to stop on the \n and ignore all the rest.

  seqLen = 0;
  qltLen = 0;

  for (; (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    if (isWhiteSpace(ch))
      continue;
    if (seqLen+1 >= seqMax)
      resizeArrayPair(seq, qlt, seqLen, seqMax, 3 * seqMax / 2);
    seq[seqLen++] = ch;
  }

  //  Skip any more whitespace, fail if we're not at a quality start, then
  //  suck in the quality line.  And then skip more whitespace.

  while (isWhiteSpace(ch))
    ch = _buffer->read();

  if (ch != '+')
    return(false);

  for (ch=_buffer->read(); (ch != '\n') && (ch != 0); ch=_buffer->read()) {
    ;
  }

  while (isWhiteSpace(ch))
    ch = _buffer->read();

  //  Read qualities and convert to integers.

  for (; (ch != '\n') && (ch != 0); ch=_buffer->read()) {
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
    _buffer->read();

  seq[seqLen] = 0;
  qlt[qltLen] = 0;

  assert(nameLen < nameMax);
  assert(seqLen  < seqMax);
  assert(qltLen  < seqMax);

  _seqIdx++;

  return(true);
}



bool
dnaSeqFile::loadSequence(char  *&name, uint32 &nameMax,
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
  //  end of the sequence).

  while (isWhiteSpace(_buffer->peek()))
    _buffer->read();

  //  If we're not at a sequence start, scan ahead to find the next one.
  //  Not bulletproof; FASTQ qv's can match this.

  if ((_buffer->peek() != '>') &&
      (_buffer->peek() != '@') &&
      (_buffer->peek() !=  0)) {
    //fprintf(stderr, "dnaSeqFile::loadSequence()-- sequence sync lost at position %lu, attempting to find the next sequence.\n", _buffer->tell());
    error |= 0x02;
  }

  bool  lastWhite = isWhiteSpace(_buffer->peek());

  while ((_buffer->peek() != '>') &&
         (_buffer->peek() != '@') &&
         (_buffer->peek() !=  0)) {
    _buffer->read();
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

    return(false);
  }

  //  If we failed to load a sequence, report an error message and zero out
  //  the sequence.  Leave the name as-is so we can at least return a length
  //  zero sequence.  If we failed to load a name, it'll still be set to NUL.

  if (loadSuccess == false) {
    //if (name[0] == 0)
    //  fprintf(stderr, "dnaSeqFile::loadSequence()-- failed to read sequence correctly at position %lu.\n", _buffer->tell());
    //else
    //  fprintf(stderr, "dnaSeqFile::loadSequence()-- failed to read sequence '%s' correctly at position %lu.\n", name, _buffer->tell());

    error |= 0x01;

    seq[0]  = 0;
    qlt[0]  = 0;
    seqLen  = 0;
  }

  return(true);
}



bool
dnaSeqFile::loadSequence(dnaSeq &seq) {
  bool result = loadSequence(seq._name, seq._nameMax,
                             seq._seq,
                             seq._qlt,  seq._seqMax, seq._seqLen, seq._error);

  if (result)
    seq.findNameAndFlags();

  return(result);
}



bool
dnaSeqFile::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength,
                      bool    &endOfSequence) {

  seqLength     = 0;
  endOfSequence = false;

  if (_buffer->eof() == true)
    return(false);

  //  If this is a new file, skip the first name line.

  if (_buffer->tell() == 0) {
    while (_buffer->peek() == '\n')    //  Skip whitespace before the first name line.
      _buffer->read();

    _buffer->skipAhead('\n', true);
  }

  //  Skip whitespace.

  while (_buffer->peek() == '\n')
    _buffer->read();

  //  If now at EOF, that's it.

  if  (_buffer->eof() == true)
    return(false);

  //  Otherwise, we must be in the middle of sequence, so load
  //  until we're not in sequence or out of space.

  while (_buffer->eof() == false) {

    //  If we're at the start of a new sequence, skip over any QV's and
    //  the next name line, set endOfSequence and return.

    if (_buffer->peek() == '>') {
      _buffer->skipAhead('\n', true);      //  Skip the name line.
      endOfSequence = true;
      return(true);
    }

    if (_buffer->peek() == '+') {
      _buffer->skipAhead('\n', true);      //  Skip the + line.
      _buffer->skipAhead('\n', true);      //  Skip the QV line.
      _buffer->skipAhead('\n', true);      //  Skip the @ line for the next sequence.
      endOfSequence = true;
      return(true);
    }

    //  Read some bases.

    seqLength += _buffer->copyUntil('\n', seq + seqLength, maxLength - seqLength);

    if (seqLength == maxLength)
      return(true);

    //  We're at a newline (or end of file), either way, suck in the next letter
    //  (or nothing) and keep going.

    _buffer->read();
  }

  //  We hit EOF.  If there are bases loaded, then we're at the end of
  //  a sequence, and should return that we loaded bases.

  endOfSequence = (seqLength > 0);

  return(endOfSequence);
}
