
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

#ifndef MERYLUTIL_BITS_BITS_V1_H
#define MERYLUTIL_BITS_BITS_V1_H

#include "types.H"

////////////////////////////////////////
//
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
namespace merylutil::inline bits::inline v1 {

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

}  //  namespace merylutil::bits::v1



////////////////////////////////////////
//
//  Return bits in a word:
//    Keeping the rightmost 64-n bits (mask out the leftmost  n bits)
//    Keeping the leftmost  64-n bits (mask out the rigthmost n bits)
//
namespace merylutil::inline bits::inline v1 {

inline uint64   clearLeftBits  (uint64 v,  uint32 l)  {  if (l >=  64) return(0);  return(v & (uint64max  >>        l));   }
inline uint64    saveLeftBits  (uint64 v,  uint32 l)  {  if (l ==   0) return(0);  return(v & (uint64max  << (64  - l)));  }
inline uint64   clearRightBits (uint64 v,  uint32 r)  {  if (r >=  64) return(0);  return(v & (uint64max  <<        r));   }
inline uint64    saveRightBits (uint64 v,  uint32 r)  {  if (r ==   0) return(0);  return(v & (uint64max  >> (64  - r)));  }

inline uint64   clearMiddleBits(uint64  v, uint32 l, uint32 r)  { return( saveRightBits(v, r) |  saveLeftBits(v, l)); }
inline uint64    saveMiddleBits(uint64  v, uint32 l, uint32 r)  { return(clearRightBits(v, r) & clearLeftBits(v, l)); }

inline uint128  clearLeftBits  (uint128 v, uint32 l)  {  if (l >= 128) return(0);  return(v & (uint128max >>        l));   }
inline uint128   saveLeftBits  (uint128 v, uint32 l)  {  if (l ==   0) return(0);  return(v & (uint128max << (128 - l)));  }
inline uint128  clearRightBits (uint128 v, uint32 r)  {  if (r >= 128) return(0);  return(v & (uint128max <<        r));   }
inline uint128   saveRightBits (uint128 v, uint32 r)  {  if (r ==   0) return(0);  return(v & (uint128max >> (128 - r)));  }

inline uint128  clearMiddleBits(uint128 v, uint32 l, uint32 r)  { return( saveRightBits(v, r) |  saveLeftBits(v, l)); }
inline uint128   saveMiddleBits(uint128 v, uint32 l, uint32 r)  { return(clearRightBits(v, r) & clearLeftBits(v, l)); }

}  //  namespace merylutil::bits::v1

////////////////////////////////////////
//
//  Rotate words
//
namespace merylutil::inline bits::inline v1 {

inline uint32   rotateBitsLeft (uint32  v, uint32 l)   { return (v << l) | (v >> (32  - l)); }
inline uint32   rotateBitsRight(uint32  v, uint32 r)   { return (v >> r) | (v << (32  - r)); }

inline uint64   rotateBitsLeft (uint64  v, uint32 l)   { return (v << l) | (v >> (64  - l)); }
inline uint64   rotateBitsRight(uint64  v, uint32 r)   { return (v >> r) | (v << (64  - r)); }

inline uint128  rotateBitsLeft (uint128 v, uint32 l)   { return (v << l) | (v >> (128 - l)); }
inline uint128  rotateBitsRight(uint128 v, uint32 r)   { return (v >> r) | (v << (128 - r)); }

}

////////////////////////////////////////
//
//  Freed, Edwin E. 1983. "Binary Magic Number" Dr. Dobbs Journal Vol. 78 (April) pp. 24-37
//    Reverse the bits in a word,
//    Count the number of set bits in a words
//    and more.
//
namespace merylutil::inline bits::inline v1 {

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

//#if (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
//#define BUILTIN_POPCOUNT
//#endif

//inline
//uint32
//countNumberOfSetBits32(uint32 x) {
//  return(__builtin_popcount(x));
//}

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

//inline
//uint64
//countNumberOfSetBits64(uint64 x) {
//  return(__builtin_popcountll(x));
//}

}  //  namespace merylutil::bits::v1



////////////////////////////////////////
//
//  Return the number of bits needed to represent 'x'.
//  It's really floor(log_2(x)) + 1.
//  Note that x=0 returns 0.
//
namespace merylutil::inline bits::inline v1 {

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

}  //  namespace merylutil::bits::v1



////////////////////////////////////////
//
//  Convert between 2- and 3-bit packed words:
//    2-bit        ddccbbaa
//    3-bit    0dd0cc0bb0aa
//
//  expandTo3() adds a zero-bit between each 2-bit word, while
//  compressTo2() drops that same bit.
//
//  At most 21 packed words can be accomodated - 21*3 = 63 bits.
//
//  Handy if you want to print 2-bit packed words as octal.
//
//  pdep/pext are probably more efficient.
//
namespace merylutil::inline bits::inline v1 {

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

}  //  namespace merylutil::bits::v1


#endif  //  MERYLUTIL_BITS_BITS_V1_H
