
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

//  stuffedBits operations on unary encoded data.

namespace merylutil::inline bits::inline v1 {

uint64
stuffedBits::getUnary(void) {
  uint64  value = 0;
  uint64  wrd;

  //  Ensure we're in valid data.

  moveToNextBlock(1);

  //  Word align us first.

  wrd = _data[_dataWrd] << (64 - _dataBit);

  //  Skip entire words.  For the first word, if we're left with only zeros
  //  after the shifting, we increase the 'value' by the number of bits
  //  we could read in the word, not the full 64 bits that are zero.

  while (wrd == 0) {
    value    += _dataBit;

    _dataPos += _dataBit;
    _dataWrd += 1;

    _dataBit  = 64;
    wrd       = _data[_dataWrd];
  }

  //  Decode the last partial word.

  wrd       = 64 - countNumberOfBits64(wrd);

  value    += wrd;

  _dataPos += wrd + 1;
  _dataBit -= wrd + 1;

  updateBit();

  return(value);
}



uint64 *
stuffedBits::getUnary(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getUnary();

  return(values);
}



uint32
stuffedBits::setUnary(uint64 value) {

  ensureSpaceInCurrentBlock(value+1);

  //  If we fit entirely within this word, handle it special.

  if (value + 1 < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + value + 1);

    _dataPos += value + 1;
    _dataBit -= value + 1;

    _data[_dataWrd] |= (uint64)1 << _dataBit;

    updateLen();

    return(value + 1);
  }

  //  We fit _exactly_ in this word, special again!

  if (value + 1 == _dataBit) {
    _data[_dataWrd]  = clearRightBits(_data[_dataWrd], _dataBit);
    _data[_dataWrd] |= 1;   //  ALWAYS the last bit.

    _dataPos += value + 1;  //  ALWAYS move to the next word.
    _dataWrd += 1;
    _dataBit  = 64;

    updateLen();

    return(value + 1);
  }

  //  Otherwise, we span at least two words.  First, get us word aligned,
  //  by clearing the rest of this word.

  assert(value + 1 > _dataBit);

  _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit);

  value    -= _dataBit;

  _dataPos += _dataBit;
  _dataWrd += 1;
  _dataBit  = 64;

  assert((_dataPos % 64) == 0);

  //  Then, set entire words to zero.

  while (value >= 64) {
    _data[_dataWrd] = 0;

    value    -= 64;

    _dataPos += 64;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Finally, set the last partial word.  This is similar to the cases
  //  at the start, but we know the bits will always be on the left of the word.

  _data[_dataWrd] = clearLeftBits(_data[_dataWrd], value + 1);

  _dataPos += value + 1;                      //  Skip the zero bits.
  _dataBit -= value + 1;                      //  (and the sentinel)

  _data[_dataWrd] |= (uint64)1 << _dataBit;   //  Add the sentinel.

  //  Update the pointers.

  updateLen();
  updateBit();

  return(value + 1);
}



uint32
stuffedBits::setUnary(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setUnary(values[ii]);

  return(size);
}

}  //  namespace merylutil::bits::v1
