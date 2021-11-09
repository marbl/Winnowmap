
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

//  To be included only by bits.H
#ifndef BITS_IMPLEMENTATIONS
#error Include bits.H instead of bits-wordArray.H
#endif


inline
uint128
wordArray::get(uint64 eIdx) {
  uint64  seg =                eIdx / _valuesPerSegment;    //  Which segment are we in?
  uint64  pos = _valueWidth * (eIdx % _valuesPerSegment);   //  Bit position of the start of the value.

  uint64  wrd = pos / 128;   //  The word we start in.
  uint64  bit = pos % 128;   //  Starting at this bit.

  uint128 val = 0;

  if (eIdx >= _validData)
    fprintf(stderr, "wordArray::get()-- eIdx %lu >= _validData %lu\n", eIdx, _validData);
  assert(eIdx < _validData);

  //  If the value is all in one word, just shift that word to the right to
  //  put the proper bits in the proper position.
  //
  //  Otherwise, the value spans two words.
  //   - Shift the first word left to place the right-most bits at the left end of the return value.
  //   - Shift the second word right so the left-most bits are at the right end of the return value.
  //
  //                       ssssssssssss <- second shift amount
  //  [--first-word--][--second-word--]
  //             [--value--]
  //                   fffff <- first shift amount

  if (bit + _valueWidth <= 128) {
    val = _segments[seg][wrd] >> (128 - _valueWidth - bit);
  }
  else {
    uint32  fShift = _valueWidth - (128 - bit);
    uint32  sShift = 128 - fShift;

    val  = _segments[seg][wrd+0] << fShift;
    val |= _segments[seg][wrd+1] >> sShift;
  }

  //  Finally, mask off the stuff we don't care about.

  val &= _valueMask;

  return(val);
}



inline
void
wordArray::setLock(void) {
  while (_lock.test_and_set(std::memory_order_relaxed) == true)
    ;
}


inline
void
wordArray::relLock(void) {
  _lock.clear();
}


inline
void
wordArray::setLock(uint64 seg, uint64 lockW1, uint64 lockW2) {

  if (lockW1 == lockW2) {
    while (_segLocks[seg][lockW1].test_and_set(std::memory_order_relaxed) == true)
      ;
  }

  else {
    while (_segLocks[seg][lockW1].test_and_set(std::memory_order_relaxed) == true)
      ;
    while (_segLocks[seg][lockW2].test_and_set(std::memory_order_relaxed) == true)
      ;
  }
}


inline
void
wordArray::relLock(uint64 seg, uint64 lockW1, uint64 lockW2) {

  if (lockW1 == lockW2) {
    _segLocks[seg][lockW1].clear();
  }
  else {
    _segLocks[seg][lockW2].clear();
    _segLocks[seg][lockW1].clear();
  }
}



inline
void
wordArray::set(uint64 eIdx, uint128 value) {
  uint64 seg =                eIdx / _valuesPerSegment;     //  Which segment are we in?
  uint64 pos = _valueWidth * (eIdx % _valuesPerSegment);    //  Which word in the segment?

  uint64 wrd = pos / 128;         //  The word we start in.
  uint64 bit = pos % 128;         //  Starting at this bit.

  uint64 lockW1 = 0;   //  Address of locks, computed inline with the
  uint64 lockW2 = 0;   //  setLock() function call below.

  //  Allocate more segment pointers and any missing segments.

  if (eIdx >= _numValuesAlloc) {
    setLock();
    allocate(eIdx);
    relLock();
  }

  //  Remember the largest element we've set.

  if (eIdx >= _validData) {
    setLock();
    _validData = std::max(_validData, eIdx+1);
    relLock();
  }

  //  Grab the locks for the two words we're going to be accessing.

  if (_wordsPerLock > 0)
    setLock(seg,
            lockW1 = (wrd + 0) / _wordsPerLock,
            lockW2 = (wrd + 1) / _wordsPerLock);

  //  Set the value in one word....
  //
  //          [--------------------]
  //                 [value]
  //           lSave           rSave
  //
  //  Or split the value across two words.
  //
  //            --lSave--               --rSave--
  //  [--word--][--first-word--][--second-word--][--word--]
  //                     [----value---=]
  //                      lSize  rSize

  value &= _valueMask;  //  Mask to width we can store, just in case.

  if (bit + _valueWidth <= 128) {
    uint32   lSave = bit;
    uint32   rSave = 128 - _valueWidth - bit;

    _segments[seg][wrd] = (saveLeftBits(_segments[seg][wrd], lSave) |
                           (value << rSave)                         |
                           saveRightBits(_segments[seg][wrd], rSave));
  }

  else {
    uint32   lSave =       bit,   rSave = 128 - _valueWidth - bit;
    uint32   lSize = 128 - bit,   rSize = _valueWidth - (128 - bit);

    _segments[seg][wrd+0] = saveLeftBits(_segments[seg][wrd+0], lSave) | (value >> rSize);
    _segments[seg][wrd+1] = (value << rSave) | saveRightBits(_segments[seg][wrd+1], rSave);
  }

  //  Release the locks.

  if (_wordsPerLock > 0)
    relLock(seg, lockW1, lockW2);
}

