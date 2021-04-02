
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

#ifndef LIBBITS_H
#define LIBBITS_H

#include "types.H"
#include "arrays.H"
#include "files.H"

#include <algorithm>
#include <atomic>


//  Writing in the middle of data is toublesome.
//  This code currently will not split an object across two allocated blocks,
//  so if you want to rewrite in the middle, you need to make sure it
//  doesn't span a previously written block.  E.g., two writes
//  of 32 bits each could be in different blocks, and trying
//  to overwrite with a single 64-bit write could result in the first
//  block getting truncated (to the current position) and the write
//  replacing the first 64-bits in the second block, not just the 32-bits expected.
//
//  And don't even think of overwriting any of the variable length data.


inline
char *
displayWord(uint64 word, char *buffer=NULL) {
  static char b[65];

  if (buffer == NULL)
    buffer = b;

  memset(buffer, 'x', 64);

  for (uint32 ii=0; ii<64; ii++)
    buffer[ii] = (word & (uint64)1 << (63 - ii)) ? '1' : '0';

  buffer[64] = 0;

  return(buffer);
};



//  Generate a bit mask on the low (0x000fff) or high bits (0xfff000).
//
//  Algorithm:
//   - set the return value to all 1's
//   - shift left or right to keep the desired numBits in the word
//   - reset to all 0's if the numBits is zero
//     (if     zero, 'r & -0' == 'r & 0000..000)
//     (if not zero, 'r & -1' == 'r & 1111..111)
//   - reset to all 1's if the numBits is larger than the number of bits in the word
//
template<typename uintType>
uintType
buildLowBitMask(uint32 numBits) {
  uintType  r;

  r   = ~((uintType)0);
  r >>= 8 * sizeof(uintType) - numBits;
  r  &= -(uintType)(numBits != 0);
  r  |= -(uintType)(numBits  > 8 * sizeof(uintType));

  return(r);
}

template<typename uintType>
uintType
buildHighBitMask(uint32 numBits) {
  uintType  r;

  r   = ~((uintType)0);
  r <<= 8 * sizeof(uintType) - numBits;
  r  &= -(uintType)(numBits != 0);
  r  |= -(uintType)(numBits  > 8 * sizeof(uintType));

  return(r);
}



//  Return bits in a word:
//    Keeping the rightmost 64-n bits (mask out the leftmost  n bits)
//    Keeping the leftmost  64-n bits (mask out the rigthmost n bits)
//
inline uint64   clearLeftBits  (uint64 v,  uint32 l)  {  if (l >=  64) return(0);  return(v & (uint64max  >>        l));   };
inline uint64    saveLeftBits  (uint64 v,  uint32 l)  {  if (l ==   0) return(0);  return(v & (uint64max  << (64  - l)));  };
inline uint64   clearRightBits (uint64 v,  uint32 r)  {  if (r >=  64) return(0);  return(v & (uint64max  <<        r));   };
inline uint64    saveRightBits (uint64 v,  uint32 r)  {  if (r ==   0) return(0);  return(v & (uint64max  >> (64  - r)));  };

inline uint64   clearMiddleBits(uint64  v, uint32 l, uint32 r)  { return( saveRightBits(v, r) |  saveLeftBits(v, l)); };
inline uint64    saveMiddleBits(uint64  v, uint32 l, uint32 r)  { return(clearRightBits(v, r) & clearLeftBits(v, l)); };

inline uint128  clearLeftBits  (uint128 v, uint32 l)  {  if (l >= 128) return(0);  return(v & (uint128max >>        l));   };
inline uint128   saveLeftBits  (uint128 v, uint32 l)  {  if (l ==   0) return(0);  return(v & (uint128max << (128 - l)));  };
inline uint128  clearRightBits (uint128 v, uint32 r)  {  if (r >= 128) return(0);  return(v & (uint128max <<        r));   };
inline uint128   saveRightBits (uint128 v, uint32 r)  {  if (r ==   0) return(0);  return(v & (uint128max >> (128 - r)));  };

inline uint128  clearMiddleBits(uint128 v, uint32 l, uint32 r)  { return( saveRightBits(v, r) |  saveLeftBits(v, l)); };
inline uint128   saveMiddleBits(uint128 v, uint32 l, uint32 r)  { return(clearRightBits(v, r) & clearLeftBits(v, l)); };



//  Freed, Edwin E. 1983. "Binary Magic Number" Dr. Dobbs Journal Vol. 78 (April) pp. 24-37
//    Reverse the bits in a word,
//    Count the number of set bits in a words
//    and more.
//
inline
uint64
reverseBits64(uint64 x) {
  x = ((x >>  1) & 0x5555555555555555llu) | ((x <<  1) & 0xaaaaaaaaaaaaaaaallu);
  x = ((x >>  2) & 0x3333333333333333llu) | ((x <<  2) & 0xccccccccccccccccllu);
  x = ((x >>  4) & 0x0f0f0f0f0f0f0f0fllu) | ((x <<  4) & 0xf0f0f0f0f0f0f0f0llu);
  x = ((x >>  8) & 0x00ff00ff00ff00ffllu) | ((x <<  8) & 0xff00ff00ff00ff00llu);
  x = ((x >> 16) & 0x0000ffff0000ffffllu) | ((x << 16) & 0xffff0000ffff0000llu);
  x = ((x >> 32) & 0x00000000ffffffffllu) | ((x << 32) & 0xffffffff00000000llu);
  return(x);
}

inline
uint32
reverseBits32(uint32 x) {
  x = ((x >>  1) & 0x55555555lu) | ((x <<  1) & 0xaaaaaaaalu);
  x = ((x >>  2) & 0x33333333lu) | ((x <<  2) & 0xcccccccclu);
  x = ((x >>  4) & 0x0f0f0f0flu) | ((x <<  4) & 0xf0f0f0f0lu);
  x = ((x >>  8) & 0x00ff00fflu) | ((x <<  8) & 0xff00ff00lu);
  x = ((x >> 16) & 0x0000fffflu) | ((x << 16) & 0xffff0000lu);
  return(x);
}


inline
uint64
uint64Swap(uint64 x) {
  x = ((x >>  8) & 0x00ff00ff00ff00ffllu) | ((x <<  8) & 0xff00ff00ff00ff00llu);
  x = ((x >> 16) & 0x0000ffff0000ffffllu) | ((x << 16) & 0xffff0000ffff0000llu);
  x = ((x >> 32) & 0x00000000ffffffffllu) | ((x << 32) & 0xffffffff00000000llu);
  return(x);
}

inline
uint32
uint32Swap(uint32 x) {
  x = ((x >>  8) & 0x00ff00fflu) | ((x <<  8) & 0xff00ff00lu);
  x = ((x >> 16) & 0x0000fffflu) | ((x << 16) & 0xffff0000lu);
  return(x);
}

inline
uint16
uint16Swap(uint16 x) {
  x = ((x >>  8) & 0x00ff) | ((x <<  8) & 0xff00);
  return(x);
}


inline
uint32
countNumberOfSetBits32(uint32 x) {
  x = ((x >>  1) & 0x55555555lu) + (x & 0x55555555lu);
  x = ((x >>  2) & 0x33333333lu) + (x & 0x33333333lu);
  x = ((x >>  4) & 0x0f0f0f0flu) + (x & 0x0f0f0f0flu);
  x = ((x >>  8) & 0x00ff00fflu) + (x & 0x00ff00fflu);
  x = ((x >> 16) & 0x0000fffflu) + (x & 0x0000fffflu);
  return(x);
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  x = ((x >>  1) & 0x5555555555555555llu) + (x & 0x5555555555555555llu);
  x = ((x >>  2) & 0x3333333333333333llu) + (x & 0x3333333333333333llu);
  x = ((x >>  4) & 0x0f0f0f0f0f0f0f0fllu) + (x & 0x0f0f0f0f0f0f0f0fllu);
  x = ((x >>  8) & 0x00ff00ff00ff00ffllu) + (x & 0x00ff00ff00ff00ffllu);
  x = ((x >> 16) & 0x0000ffff0000ffffllu) + (x & 0x0000ffff0000ffffllu);
  x = ((x >> 32) & 0x00000000ffffffffllu) + (x & 0x00000000ffffffffllu);
  return(x);
}


//  Return the number of bits needed to represent 'x'.
//  It's really floor(log_2(x)) + 1.
//  Note that x=0 returns 0.
//
inline
uint32
countNumberOfBits32(uint32 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return(countNumberOfSetBits32(x));
}

inline
uint64
countNumberOfBits64(uint64 x) {
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;
  return(countNumberOfSetBits64(x));
}




//#if (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
//#define BUILTIN_POPCOUNT
//#endif

#ifdef BUILTIN_POPCOUNT

inline
uint32
countNumberOfSetBits32(uint32 x) {
  return(__builtin_popcount(x));
}

inline
uint64
countNumberOfSetBits64(uint64 x) {
  return(__builtin_popcountll(x));
}

#endif



//  Expand a 2-bit packed word into a 3-bit packed word.
//    input        aabbccdd
//    output   0aa0bb0cc0dd
//  Handy if you want to print such a number as octal.
//
inline
uint64
expandTo3(uint64 v) {
  uint64 o = 0;

  o  = (v & 0x0000000000000003llu) << 0;
  o |= (v & 0x000000000000000cllu) << 1;
  o |= (v & 0x0000000000000030llu) << 2;
  o |= (v & 0x00000000000000c0llu) << 3;
  o |= (v & 0x0000000000000300llu) << 4;
  o |= (v & 0x0000000000000c00llu) << 5;
  o |= (v & 0x0000000000003000llu) << 6;
  o |= (v & 0x000000000000c000llu) << 7;
  o |= (v & 0x0000000000030000llu) << 8;
  o |= (v & 0x00000000000c0000llu) << 9;
  o |= (v & 0x0000000000300000llu) << 10;
  o |= (v & 0x0000000000c00000llu) << 11;
  o |= (v & 0x0000000003000000llu) << 12;
  o |= (v & 0x000000000c000000llu) << 13;
  o |= (v & 0x0000000030000000llu) << 14;
  o |= (v & 0x00000000c0000000llu) << 15;
  o |= (v & 0x0000000300000000llu) << 16;
  o |= (v & 0x0000000c00000000llu) << 17;
  o |= (v & 0x0000003000000000llu) << 18;
  o |= (v & 0x000000c000000000llu) << 19;
  o |= (v & 0x0000030000000000llu) << 20;
  //   (v & 0x00000c0000000000llu) << 21;  //  This overflows.

  assert((v & 0xffffc0000000000llu) == 0);

  return(o);
}


//  Compress a 3-bit packed word into a 2-bit packed word, dropping the high bit.
inline
uint64
compressTo2(uint64 v) {
  uint64 o = 0;

  o  = (v & 0x0000000000000003llu) >> 0;
  o |= (v & 0x0000000000000018llu) >> 1;
  o |= (v & 0x00000000000000c0llu) >> 2;
  o |= (v & 0x0000000000000600llu) >> 3;
  o |= (v & 0x0000000000003000llu) >> 4;
  o |= (v & 0x0000000000018000llu) >> 5;
  o |= (v & 0x00000000000c0000llu) >> 6;
  o |= (v & 0x0000000000600000llu) >> 7;
  o |= (v & 0x0000000003000000llu) >> 8;
  o |= (v & 0x0000000018000000llu) >> 9;
  o |= (v & 0x00000000c0000000llu) >> 10;
  o |= (v & 0x0000000600000000llu) >> 11;
  o |= (v & 0x0000003000000000llu) >> 12;
  o |= (v & 0x0000018000000000llu) >> 13;
  o |= (v & 0x00000c0000000000llu) >> 14;
  o |= (v & 0x0000600000000000llu) >> 15;
  o |= (v & 0x0003000000000000llu) >> 16;
  o |= (v & 0x0018000000000000llu) >> 17;
  o |= (v & 0x00c0000000000000llu) >> 18;
  o |= (v & 0x0600000000000000llu) >> 19;
  o |= (v & 0x3000000000000000llu) >> 20;

  assert((o & 0xffffc0000000000llu) == 0);

  return(o);
}




class bitArray {
public:
  bitArray(uint64 maxNumBits=0) {
    _maxBitSet   = 0;
    _maxBitAvail = maxNumBits;
    _bits        = NULL;

    if (maxNumBits > 0)
      allocate(maxNumBits);
  };

  ~bitArray(void) {
    delete [] _bits;
  };

  bool     isAllocated(void) {
    return(_bits != NULL);
  }

  void     allocate(uint64 maxNumBits) {
    if (isAllocated() == true)
      return;

    _maxBitSet = 0;
    _maxBitAvail = maxNumBits;
    _bits        = new uint64 [_maxBitAvail / 64 + 1];

    clear();
  };

  void     clear(void) {
    memset(_bits, 0, sizeof(uint64) * (_maxBitAvail / 64 + 1));
  };

  bool     getBit(uint64 position) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);

    if (_maxBitAvail <= position)
      fprintf(stderr, "getBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
              position, _maxBitAvail);
    assert(position < _maxBitAvail);

    return((_bits[w] >> b) & 0x00000001);
  };

  void     setBit(uint64 position, bool value) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);
    uint64   m = ((uint64)1) << b;

    //fprintf(stderr, "SET pos %9" F_U64P " word %2" F_U64P " bit %2" F_U64P " value %c  0x%016" F_X64P " -> ",
    //        position, w, b, (value) ? '1' : '0', _bits[w]);

    if (_maxBitAvail <= position)
      fprintf(stderr, "setBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
              position, _maxBitAvail);
    assert(position < _maxBitAvail);

    _bits[w] &= ~m;
    _bits[w] |= ((uint64)value) << b;

    //fprintf(stderr, "0x%016" F_X64P "\n", _bits[w]);
  };

  //  Returns state of bit before flipping.
  bool     flipBit(uint64 position) {
    uint64   w =      (position / 64);
    uint64   b = 63 - (position % 64);
    uint64   m = ((uint64)1) << b;

    if (_maxBitAvail <= position)
      fprintf(stderr, "flipBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
              position, _maxBitAvail);
    assert(position < _maxBitAvail);

    uint64   v = _bits[w] & m;

    //fprintf(stderr, "FLIP w %lu b %lu m 0x%016lx v 0x%016lx FROM 0x%016lx", w, b, m, v, _bits[w]);

    _bits[w] ^= m;

    //fprintf(stderr, " TO 0x%016lx\n", _bits[w]);

    return(v >> b);
  };

private:
  uint64   _maxBitSet;
  uint64   _maxBitAvail;
  uint64  *_bits;
};



////////////////////////////////////////
//
//  wordArray - An array that efficiently stores non-machine-word size
//  integer words by packing the bits into machine-size words.  The array is
//  variable length but not sparse - accessing element 1,000,000 will
//  allocate elements 0 through 999,999.
//
//  The size, in bits, of each element is set at construction time.  All
//  elements must be the same size.
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
class wordArray {
public:
  wordArray(uint32 valueWidth, uint64 segmentsSizeInBits, bool useLocks);
  ~wordArray();

  void     clear(void);                   //  Reset the array to zero, doesn't deallocate space.

  void     allocate(uint64 nElements);    //  Pre-allocate space for nElements.

  uint128  get(uint64 eIdx);              //  Get the value of element eIdx.
  void     set(uint64 eIdx, uint128 v);   //  Set the value of element eIdx to v.

public:
  void     show(void);                    //  Dump the wordArray to the screen; debugging.

private:
  void     setLock(uint64 seg, uint64 lockW1, uint64 lockW2);
  void     relLock(uint64 seg, uint64 lockW1, uint64 lockW2);
  void     setNval(uint32 eIdx);

private:
  uint64              _valueWidth       = 0;         //  Width of the values stored.
  uint64              _valueMask        = 0;         //  Mask the low _valueWidth bits
  uint64              _segmentSize      = 0;         //  Size, in bits, of each block of data.

  uint64              _valuesPerSegment = 0;         //  Number of values in each block.

  uint64              _wordsPerSegment  = 0;         //  Number of 128-bit words in each segment
  uint64              _wordsPerLock     = 0;         //  How many words are covered by each lock.
  uint64              _locksPerSegment  = 0;         //  Number of locks per segment

  uint64              _numValues        = 0;         //  Number of values stored in the array.
  std::atomic_flag    _numValuesLock;                //  Lock on the above.

  uint64              _segmentsLen      = 0;         //  Number of blocks in use.
  uint64              _segmentsMax      = 0;         //  Number of block pointers allocated.
  uint128           **_segments         = nullptr;   //  List of blocks allocated.

  std::atomic_flag  **_segLocks         = nullptr;   //  Locks on pieces of the segments.
};



class stuffedBits {
public:
  stuffedBits(uint64 nBits=16 * 1024 * 1024 * 8);
  stuffedBits(const char *inputName);
  stuffedBits(FILE *inFile);
  stuffedBits(readBuffer *B);
  //stuffedBits(stuffedBits &that);   //  Untested.
  ~stuffedBits();

  //  Debugging.

  char    *displayWord(uint64 w) {
    return(::displayWord(_data[w]));
  };

  //  Files.

  void     dumpToBuffer(writeBuffer *B);
  bool     loadFromBuffer(readBuffer *B);

  void     dumpToFile(FILE *F);
  bool     loadFromFile(FILE *F);

  //  Management of the read/write head.

  void     setPosition(uint64 position, uint64 length = 0);
  uint64   getPosition(void);

  uint64   getLength(void);

  void     byteAlign(void);

  //  SINGLE BITS

  bool     getBit(void);           //  get a bit.
  bool     testBit(void);          //  get a bit, but don't move position.
  void     setBit(bool on=true);   //  set a bit.

  //  UNARY CODED DATA

  uint64   getUnary(void);
  uint64  *getUnary(uint64 number, uint64 *values);

  uint32   setUnary(uint64 value);
  uint32   setUnary(uint64 number, uint64 *values);

  //  BINARY CODED DATA

  uint64   getBinary(uint32 width);
  uint64  *getBinary(uint32 width, uint64 number, uint64 *values=NULL);

  uint32   setBinary(uint32 width, uint64 value);
  uint32   setBinary(uint32 width, uint64 number, uint64 *values);

  //  ELIAS GAMMA CODED DATA

  uint64   getEliasGamma(void);
  uint64  *getEliasGamma(uint64 number, uint64 *values=NULL);

  uint32   setEliasGamma(uint64 value);
  uint32   setEliasGamma(uint64 number, uint64 *values);

  //  ELIAS DELTA CODED DATA

  uint64   getEliasDelta(void);
  uint64  *getEliasDelta(uint64 number, uint64 *values=NULL);

  uint32   setEliasDelta(uint64 value);
  uint32   setEliasDelta(uint64 number, uint64 *values);

  //  ELIAS OMEGA CODED DATA - the omega code looks hard to implement - the
  //  encoding and decoding streams are backwards from each other.  The idea
  //  is:
  //
  //    push the binary representation of the value onto a stack.
  //    set value to one less than the number of bits emitted last.
  //    loop until the value is 1.
  //
  //  The stream is constructed by emitting the words on the stack, and
  //  terminating the stream with an extra '0'.
  //
#if 0
  uint64   getEliasOmega(void);
  uint64  *getEliasOmega(uint64 number, uint64 *values=NULL);

  uint32   setEliasOmega(uint64 value);
  uint32   setEliasOmega(uint64 number, uint64 *values);
#endif

  //  GOLOMB CODED DATA
  //
  //  Pick m.  For m == power_of_two, this is RICE CODED DATA.
  //
  //    q = floow(n/m).
  //    r = n-qm
  //    c = ceil(log_2 m)
  //
  //  Unary encode q, binary encode r.
  //
  //  The first 2^c-m values are encoded as c-1 bit values, starting with 00...00,
  //  The rest as c-bit numbers, ending with 11...11
  //


  //  FIBONACCI CODED DATA
  //
  //  A Fibonacci number is F(n) = F(n-1) + F(n-2), wher F(0) = F(1) = 1.
  //
  //  The Zeckendorf representation of a number encodes it such that no
  //  two consecurive Fibonacci numbers are used.  From the definition
  //  of a Fibonacci number, any pattern "100" can be replaced with "011".
  //  A number encoded after this transofmration is in the Fibonacci
  //  representation ("Zeckendorf representation" seems to be a real thing,
  //  I just made up "Fibonacci representation" - the two terms seem to
  //  be used synonymously in the real world).
  //
  //  Once encoded, it's added to the bit stream reversed.
  //
  //  For the Zeckendorf representation, a single 1-bit is added, terminating
  //  the number with the last '1' bit of data, followed immediately by
  //  another '1' bit.  (Because, by definition, there are no two adjacent
  //  set bits in the encoded number).
  //
  //  For the Fibonacci representation, we need to append two '0' bits.
  //  (Because, by definition, there are no two adjacent unset bits in the
  //  representation).  BUT, this representation saves at most one bit
  //  (replacing 100 at the start of the string by 011), but the savings
  //  is lost by the extra stop bit we need.
  //
  uint64   getZeckendorf(void);
  uint64  *getZeckendorf(uint64 number, uint64 *values=NULL);

  uint32   setZeckendorf(uint64 value);
  uint32   setZeckendorf(uint64 number, uint64 *values);

  //  Old meryl uses preDecrement() when using compressed bucket counting.
  //  Nothing else in canu uses these, and they're painful, so left unimplemented.
#if 0
  uint64   preIncrementBinary(uint64 width, uint64 position);
  uint64   postIncrementBinary(uint64 width, uint64 position);
  uint64   preDecrementBinary(uint64 width, uint64 position);
  uint64   postDecrementBinary(uint64 width, uint64 position);
#endif


private:

  //  For writing, update the length of the block to the maximum of where we're at now and the existing length.
  //
  void      updateLen(void) {
    _dataBlockLen[_dataBlk] = std::max(_dataPos, _dataBlockLen[_dataBlk]);
  };

  //  For both reading and writing, move to the next word if we're at the end of the current one.
  //
  void      updateBit(void) {
    if (_dataBit == 0) {
      _dataWrd += 1;
      _dataBit  = 64;
    }
  };

  //  For reading operations, move to the next block if we're at the end of the current one.
  //  For writing operations, this is done before the write, in ensureSpace().
  //
  //  Should be done before any reading operation.  It isn't (strictly) needed at the end
  //  of a read.  The next read will just do it at the start, and the next write
  //  only cares that there is enough space in the current block for the write.
  //
  void      updateBlk(uint64 readLength) {

    assert(_dataBit != 0);
    assert(_dataBit <= 64);

    //  If in this block, nothing to do.

    //fprintf(stderr, "updateBlk()--  at _dataPos " F_U64 "  _dataBlockLen[%u] = " F_U64 "\n", _dataPos, _dataBlk, _dataBlockLen[_dataBlk]);

    if (_dataPos + readLength <= _dataBlockLen[_dataBlk])
      return;

    //  Otherwise, we MUST be at the end of the current block.  If not,
    //  whatever we're trying to read is in the next block (or we're trying
    //  to read something longer than what is here).

    assert(_dataPos == _dataBlockLen[_dataBlk]);

    _dataBlk += 1;

    if (_dataBlk >= _dataBlocksLen)
      fprintf(stderr, "ERROR: _dataBlk = %lu  _dataBlocksLen = %u\n", _dataBlk, _dataBlocksLen);
    assert(_dataBlk < _dataBlocksLen);

    _dataPos  = 0;
    _data     = _dataBlocks[_dataBlk];

    _dataWrd  = 0;
    _dataBit  = 64;
  }

  void     clearBlock(void) {
    for (uint64 ii=0; ii<_dataBlockLenMaxW; ii++)
      _data[ii] = 0;
  };
  
  //  For writing operations, make sure there is enough space for the write in this block.
  //
  void     ensureSpace(uint64 spaceNeeded) {

    assert(_dataBit != 0);
    assert(_dataBit <= 64);

    //  If enough space in the current block, just return.

    if (_dataPos + spaceNeeded < _dataBlockLenMaxB)
      return;

    //  Othwerwise, terminate the current block.

    _dataBlockLen[_dataBlocksLen - 1] = _dataPos;

    //  Move to the new block.

    _dataBlocksLen++;

    if (_dataBlocksLen >= _dataBlocksMax) {
      setArraySize(_dataBlocks,   _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
      setArraySize(_dataBlockBgn, _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
      setArraySize(_dataBlockLen, _dataBlocksLen, _dataBlocksMax, _dataBlocksLen + 128);
    }

    assert(spaceNeeded <= _dataBlockLenMaxB);

    _dataPos = 0;
    _data    = _dataBlocks[_dataBlocksLen - 1] = new uint64 [_dataBlockLenMaxW];

    clearBlock();

    _dataBlockBgn[_dataBlocksLen - 1] = _dataBlockBgn[_dataBlocksLen - 2] + _dataBlockLen[_dataBlocksLen - 2];
    _dataBlockLen[_dataBlocksLen - 1] = 0;

    _dataBlk += 1;
    _dataWrd  = 0;
    _dataBit  = 64;
  };

  uint64    bitsToWords(uint64 bits) {
    return(bits / 64 + ((bits % 64) ? 1 : 0));
  };

  uint64   _dataBlockLenMaxB;  //  Allocated length of each block (in BITS).
  uint64   _dataBlockLenMaxW;  //  Allocated length of each block (in WORDS). 

  uint32   _dataBlocksLen;     //  Number of allocated data blocks.
  uint32   _dataBlocksMax;     //  Number of blocks we can allocate.

  uint64  *_dataBlockBgn;      //  Starting position, in the global file, of this block.
  uint64  *_dataBlockLen;      //  Length of this block.
  uint64 **_dataBlocks;        //  Just piles of bits.  Nothing interesting here.

  uint64   _dataPos;           //  Position in this block, in BITS.
  uint64  *_data;              //  Pointer to the currently active data block.

  uint64   _dataBlk;           //  Active data block.
  uint64   _dataWrd;           //  Active word in the active data block.
  uint64   _dataBit;           //  Active bit in the active word in the active data block (aka, number of bits left in this word)

  uint64   _fibData[93];       //  A pile of Fibonacci numbers.
};


//  Implementations.

#define BITS_IMPLEMENTATIONS

#include "bits-wordArray.H"

#undef BITS_IMPLEMENTATIONS


#endif  //  LIBBITS_H
