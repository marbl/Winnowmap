
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

#include "sequence-v1.H"

namespace merylutil::inline sequence::inline v1 {

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


//  Decode an encoded sequence (in chunk) of length chunkLen.
//  seq must be allocated to have seqLen+1 bytes.
//  seqLen must be the length of the sequence to decode.
//
void
decode2bitSequence(uint8 const *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {
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
encode2bitSequence(uint8 *&chunk, char const *seq, uint32 seqLen) {

  assert(seq[seqLen-1] != 0);
  assert(seq[seqLen]   == 0);

  for (uint32 ii=0; ii<seqLen; ii++) {       //  If non-ACGT present, return
    char  base = seq[ii];                    //  0 to indicate we can't encode.

    if ((base < 0x20) ||
        (base > 0x7e)) {
      fprintf(stderr, "Invalid base 0x%02x detected at position %u\n", base, ii);
      return(0);
    }

    if ((base != 'a') && (base != 'A') &&
        (base != 'c') && (base != 'C') &&
        (base != 'g') && (base != 'G') &&
        (base != 't') && (base != 'T')) {
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
decode3bitSequence(uint8 const *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {
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


//  Encode a sequence into chunk.  The length of the chunk in bytes is returned.
//  If the chunk is NULL, it is allocated.  Otherwise, it must be
//  at least seqLen bytes in length.
//
uint32
encode3bitSequence(uint8 *&chunk, char const *seq, uint32 seqLen) {

  assert(seq[seqLen-1] != 0);
  assert(seq[seqLen]   == 0);

  for (uint32 ii=0; ii<seqLen; ii++) {       //  If non-ACGTN present, return
    char  base = seq[ii];                    //  0 to indicate we can't encode.

    if ((base < 0x20) ||
        (base > 0x7e)) {
      fprintf(stderr, "Invalid base 0x%02x detected at position %u\n", base, ii);
      return(0);
    }

    if ((base != 'a') && (base != 'A') &&
        (base != 'c') && (base != 'C') &&
        (base != 'g') && (base != 'G') &&
        (base != 't') && (base != 'T') &&
        (base != 'n') && (base != 'N')) {
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
decode8bitSequence(uint8 const *chunk, uint32 chunkLen, char *seq, uint32 seqLen) {

  assert(seq != NULL);

  for (uint32 ii=0; ii<seqLen; ii++)
    seq[ii] = chunk[ii];

  seq[seqLen] = 0;
}


uint32
encode8bitSequence(uint8 *&chunk, char const *seq, uint32 seqLen) {

  if (chunk == NULL)
    chunk = new uint8 [ seqLen ];

  for (uint32 ii=0; ii<seqLen; ii++)
    chunk[ii] = seq[ii];

  return(seqLen);
}

}  //  merylutil::sequence::v1
