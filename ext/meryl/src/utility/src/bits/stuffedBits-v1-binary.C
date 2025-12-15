
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
#include "stuffedBits-v1.H"

//  stuffedBits operations on binary coded data.

namespace merylutil::inline bits::inline v1 {

uint64
stuffedBits::getBinary(uint32 width) {
  uint64  value = 0;

  if (width == 0)
    return(0);

  assert(width < 65);

  //  Ensure we're in valid data.

  moveToNextBlock(width);

  //  If we're contained in a single word, special case.

  if (width < _dataBit) {
    value = saveRightBits(_data[_dataWrd] >> (_dataBit - width), width);

    _dataPos += width;
    _dataBit -= width;
  }

  //  If we're exactly in a single word, another special case.

  else if (width == _dataBit) {
    value = saveRightBits(_data[_dataWrd], width);

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we're spanning two words.

  else {
    uint64  w1 =         _dataBit;
    uint64  w2 = width - _dataBit;

    uint64  l  = saveRightBits(_data[_dataWrd + 0], w1) <<       w2;
    uint64  r  = saveLeftBits (_data[_dataWrd + 1], w2) >> (64 - w2);

    value = l | r;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  return(value);
}



uint64 *
stuffedBits::getBinary(uint32 width, uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getBinary(width);

  return(values);
}



uint32
stuffedBits::setBinary(uint32 width, uint64 value) {

  if (width == 0)
    return(0);

  assert(width < 65);

  ensureSpaceInCurrentBlock(width);

  //  Mask off any pieces we're not supposed to be seeing.

  value = saveRightBits(value, width);

  //  If we fit entirely within this word, handle it special.

  if (width < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + width) | (value << (_dataBit - width));

    _dataPos += width;
    _dataBit -= width;
  } else

  //  We fit _exactly_ in this word, special again!

  if (width == _dataBit) {
    _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit) | value;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we span two words.

  else {
    uint32  w1 =         _dataBit;
    uint32  w2 = width - _dataBit;

    _data[_dataWrd + 0] = clearRightBits(_data[_dataWrd + 0], w1) | (value >> (     w2));
    _data[_dataWrd + 1] = clearLeftBits (_data[_dataWrd + 1], w2) | (value << (64 - w2));

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  //  updateBit() isn't needed; it's handled in the special cases.

  updateLen();

  return(width);
}



uint32
stuffedBits::setBinary(uint32 width, uint64 number, uint64 *values) {
  uint32 size = 0;

  if (width == 0)
    return(0);

  for (uint64 ii=0; ii<number; ii++)
    size += setBinary(width, values[ii]);

  return(size);
}

}  //  namespace merylutil::bits::v1
