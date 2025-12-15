
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

#ifndef MERYLUTIL_BITS_WORDARRAY_V1_H
#define MERYLUTIL_BITS_WORDARRAY_V1_H

#include <atomic>

#include "types.H"

//
//  wordArray - An array that efficiently stores non-machine-word size
//  integer words by packing the bits into machine-size words.  The array is
//  variable length but not sparse - accessing element 1,000,000 will
//  allocate elements 0 through 999,999.
//
//  The size, in bits, of each element is set at construction time.  All
//  elements must be the same size.
//
//  Words of at most 128 bits in size can be stored.
//
//  The elements are stored in a set of fixed-size blocks.  The block size
//  can also be set at construction time.  Note that this is specified IN
//  BITS.  The default size is 64 KB per block.  Decrease this if you know
//  you only need a few KB to store all values, or if you are storing several
//  GB of data.  There is no real performance loss/gain; it just adjusts the
//  number of blocks allocated.  There might be a slight degradation in
//  performance of the memory management system if millions of blocks are
//  allocated.
//

//
//  At the default segmentSize of 64 KB = 524288 bits, we'll allocate 4096
//  128-bit words per segment.  With _wordsPerLock = 64, we'll then have
//  4096 / 64 = 64+1 locks per segment.
//
//  Note that 'values' refers to the user-supplied data of some small size,
//  while 'words' are the 128-bit machine words used to store the data.
//

namespace merylutil::inline bits::inline v1 {

class wordArray {
public:
  wordArray(uint32 valueWidth, uint64 segmentsSizeInBits, bool useLocks);
  ~wordArray();

  void     erase(uint8 c, uint64 maxElt); //  Clear allocated space to c, set maxElement to maxElt.

  void     allocate(uint64 nElements);    //  Pre-allocate space for nElements.

  uint128  get(uint64 eIdx);              //  Get the value of element eIdx.
  void     set(uint64 eIdx, uint128 v);   //  Set the value of element eIdx to v.

public:
  void     show(void);                    //  Dump the wordArray to the screen; debugging.

private:
  void     setLock(void);
  void     relLock(void);
  void     setLock(uint64 seg, uint64 lockW1, uint64 lockW2);
  void     relLock(uint64 seg, uint64 lockW1, uint64 lockW2);
  void     setNval(uint32 eIdx);

private:
  uint64              _valueWidth       = 0;         //  Width of the values stored.
  uint128             _valueMask        = 0;         //  Mask the low _valueWidth bits
  uint64              _segmentSize      = 0;         //  Size, in bits, of each block of data.

  uint64              _valuesPerSegment = 0;         //  Number of values in each block.

  uint64              _wordsPerSegment  = 0;         //  Number of 128-bit words in each segment
  uint64              _wordsPerLock     = 0;         //  How many words are covered by each lock.
  uint64              _locksPerSegment  = 0;         //  Number of locks per segment

  uint64              _numValuesAlloc   = 0;
  uint64              _validData        = 0;

  std::atomic_flag    _lock;                         //  Global lock

  uint64              _segmentsLen      = 0;         //  Number of blocks in use.
  uint64              _segmentsMax      = 0;         //  Number of block pointers allocated.
  uint128           **_segments         = nullptr;   //  List of blocks allocated.

  std::atomic_flag  **_segLocks         = nullptr;   //  Locks on pieces of the segments.
};


inline
void
wordArray::setLock(void) {
  while (_lock.test_and_set(std::memory_order_acquire) == true)
    ;
}


inline
void
wordArray::relLock(void) {
  _lock.clear(std::memory_order_release);
}


inline
void
wordArray::setLock(uint64 seg, uint64 lockW1, uint64 lockW2) {

  if (lockW1 == lockW2) {
    while (_segLocks[seg][lockW1].test_and_set(std::memory_order_acquire) == true)
      ;
  }
  else {
    while (_segLocks[seg][lockW1].test_and_set(std::memory_order_acquire) == true)
      ;
    while (_segLocks[seg][lockW2].test_and_set(std::memory_order_acquire) == true)
      ;
  }
}


inline
void
wordArray::relLock(uint64 seg, uint64 lockW1, uint64 lockW2) {

  if (lockW1 == lockW2) {
    _segLocks[seg][lockW1].clear(std::memory_order_release);
  }
  else {
    _segLocks[seg][lockW2].clear(std::memory_order_release);
    _segLocks[seg][lockW1].clear(std::memory_order_release);
  }
}

}  //  namespace merylutil::bits::v1

#endif  //  MERYLUTIL_BITS_WORDARRAY_V1_H
