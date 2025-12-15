
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

#ifndef MERYLUTIL_BITS_BITARRAY_V1_H
#define MERYLUTIL_BITS_BITARRAY_V1_H

#include "types.H"

namespace merylutil::inline bits::inline v1 {

class bitArray {
public:
  bitArray(uint64 maxNumBits=0)   {  allocate(maxNumBits);      }
  ~bitArray(void)                 {  delete [] _bits;           }

  bool     isAllocated(void)      {  return(_bits != nullptr);  }

  void     allocate(uint64 maxNumBits);
  void     clear(void)            {  memset(_bits, 0, sizeof(uint64) * (_maxBitAvail / 64 + 1));  }

  bool     getBit(uint64 position);              //  Returns state of bit.
  void     setBit(uint64 position, bool value);  //  Sets bit to 'value'.
  bool     flipBit(uint64 position);             //  Returns state of bit before flipping.

private:
  uint64   _maxBitSet   = 0;
  uint64   _maxBitAvail = 0;
  uint64  *_bits        = nullptr;
};


inline
void
bitArray::allocate(uint64 maxNumBits) {

  if ((isAllocated() == true) || (maxNumBits == 0))
    return;

  _maxBitSet = 0;
  _maxBitAvail = maxNumBits;
  _bits        = new uint64 [_maxBitAvail / 64 + 1];

  clear();
}

inline
bool
bitArray::getBit(uint64 position) {
  uint64   w =      (position / 64);
  uint64   b = 63 - (position % 64);

  if (_maxBitAvail <= position)
    fprintf(stderr, "getBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
            position, _maxBitAvail);
  assert(position < _maxBitAvail);

  return((_bits[w] >> b) & 0x00000001);
}

inline
void
bitArray::setBit(uint64 position, bool value) {
  uint64   w =      (position / 64);
  uint64   b = 63 - (position % 64);
  uint64   m = ((uint64)1) << b;

  if (_maxBitAvail <= position)
    fprintf(stderr, "setBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
            position, _maxBitAvail);
  assert(position < _maxBitAvail);

  _bits[w] &= ~m;
  _bits[w] |= ((uint64)value) << b;
}

inline
bool
bitArray::flipBit(uint64 position) {
  uint64   w =      (position / 64);
  uint64   b = 63 - (position % 64);
  uint64   m = ((uint64)1) << b;

  if (_maxBitAvail <= position)
    fprintf(stderr, "flipBit()--  ERROR: position=" F_U64 " > maximum available=" F_U64 "\n",
            position, _maxBitAvail);
  assert(position < _maxBitAvail);

  uint64   v = _bits[w] & m;

  _bits[w] ^= m;

  return(v >> b);
}

}  //  namespace merylutil::bits::v1

#endif  //  MERYLUTIL_BITS_BITARRAY_V1_H
