
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

#include "bits.H"


//
//  At the default segmentSize of 64 KB = 524288 bits, we'll allocate 4096
//  128-bit words per segment.  With _wordsPerLock = 64, we'll then have
//  4096 / 64 = 64+1 locks per segment.
//
//  Note that 'values' refers to the user-supplied data of some small size,
//  while 'words' are the 128-bit machine words used to store the data.
//

wordArray::wordArray(uint32 valueWidth, uint64 segmentSizeInBits, bool useLocks) {

  if (valueWidth == 0)
    fprintf(stderr, "wordArray::wordArray()-- valueWidth=%u too large; must greater than zero.\n", valueWidth), exit(1);
  assert(valueWidth > 0);

  if (valueWidth > 128)
    fprintf(stderr, "wordArray::wordArray()-- valueWidth=%u too large; must be at most 128.\n", valueWidth), exit(1);
  assert(valueWidth <= 128);

  _valueWidth       = valueWidth;          //  In bits.
  _valueMask        = buildLowBitMask<uint128>(_valueWidth);
  _segmentSize      = segmentSizeInBits;   //  In bits.

  _valuesPerSegment = _segmentSize / _valueWidth;

  _wordsPerSegment  = _segmentSize / 128;
  _wordsPerLock     = (useLocks == false) ? (0) : (64);
  _locksPerSegment  = (useLocks == false) ? (0) : (_segmentSize / 128 / _wordsPerLock + 1);

  _numValuesAlloc   = 0;
  _validData        = 0;

  _lock.clear();

  _segmentsLen      = 0;
  _segmentsMax      = 16;
  _segments         = new uint128 *          [_segmentsMax];
  _segLocks         = new std::atomic_flag * [_segmentsMax];

  for (uint32 ss=0; ss<_segmentsMax; ss++) {
    _segments[ss] = nullptr;
    _segLocks[ss] = nullptr;
  }
}



wordArray::~wordArray() {
  for (uint32 i=0; i<_segmentsLen; i++) {
    delete [] _segments[i];
    delete [] _segLocks[i];
  }

  delete [] _segments;
  delete [] _segLocks;
}



//  Erase ALL allocated space to the given constant,
//  then sets the array size to be maxElt.
//
void
wordArray::erase(uint8 c, uint64 maxElt) {

  allocate(maxElt);

  for (uint32 seg=0; seg<_segmentsLen; seg++)
    memset(_segments[seg], c, sizeof(uint128) * _wordsPerSegment);

  _validData = maxElt;
}



//  Allocate space for at least nElements.  Space will be allocated for an
//  integer number of blocks, each block having _valuesPerSegment elements.
//
//  Does not ever shrink space or change the number of elements in the array.
//
void
wordArray::allocate(uint64 nElements) {
  uint64 segmentsNeeded = nElements / _valuesPerSegment + 1;

  //  Allocate more space for segment pointers.  Does nothing
  //  if segmentsNeeded <= _segmentsMax.

  resizeArrayPair(_segments,
                  _segLocks,
                  _segmentsLen, _segmentsMax, segmentsNeeded,
                  _raAct::copyData | _raAct::clearNew);

  //  Allocate segments and locks, then open the locks.

  for (uint32 seg=_segmentsLen; seg<segmentsNeeded; seg++) {
    if (_segments[seg] != nullptr)
      continue;

    _segments[seg] = new uint128 [ _wordsPerSegment ];

    if (_locksPerSegment > 0) {
      _segLocks[seg] = new std::atomic_flag [ _locksPerSegment ];

      for (uint32 ll=0; ll<_locksPerSegment; ll++)
        _segLocks[seg][ll].clear();
    }
  }

  //  Update the number of values/segments we have allocated.

  if (segmentsNeeded > _segmentsLen) {
    _numValuesAlloc = segmentsNeeded * _valuesPerSegment;
    _segmentsLen    = segmentsNeeded;
  }
}



void
wordArray::show(void) {
  uint64  lastBit = _validData * _valueWidth;

  fprintf(stderr, "wordArray:\n");
  fprintf(stderr, "  validData        %10lu values\n", _validData);
  fprintf(stderr, "  valueWidth       %10lu bits\n",   _valueWidth);
  fprintf(stderr, "  segmentSize      %10lu bits\n",   _segmentSize);
  fprintf(stderr, "  valuesPerSegment %10lu values\n", _valuesPerSegment);
  fprintf(stderr, "\n");

  //  For each segment, dump full words, until we hit the end of data.

  for (uint64 ss=0; ss<_segmentsLen; ss++) {
    fprintf(stderr, "Segment %lu:\n", ss);

    uint64 bitPos = ss * _valuesPerSegment * _valueWidth;

    for (uint64 ww=0; (ww < _wordsPerSegment) && (bitPos < lastBit); ww += 4) {
      fprintf(stderr, "%5lu: %s %s %s %s\n",
              ww,
              (bitPos + 128 * 0 < lastBit) ? toHex(_segments[ss][ww+0]) : "",
              (bitPos + 128 * 1 < lastBit) ? toHex(_segments[ss][ww+1]) : "",
              (bitPos + 128 * 2 < lastBit) ? toHex(_segments[ss][ww+2]) : "",
              (bitPos + 128 * 3 < lastBit) ? toHex(_segments[ss][ww+3]) : "");

      bitPos += 128 * 4;
    }
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
}
