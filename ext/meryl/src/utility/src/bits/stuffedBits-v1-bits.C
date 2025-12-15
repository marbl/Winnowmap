
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

//  stuffedBits operations on single bits.  getBit() and setBit() are
//  special cases on getBinary() and setBinary().

namespace merylutil::inline bits::inline v1 {

bool
stuffedBits::getBit(void) {
  moveToNextBlock(1);

  bool value = saveRightBits(_data[_dataWrd] >> (_dataBit - 1), 1);

  _dataPos++;
  _dataBit--;

  updateBit();

  return(value);
}



bool
stuffedBits::testBit(void) {
  moveToNextBlock(1);

  bool value = saveRightBits(_data[_dataWrd] >> (_dataBit - 1), 1);

  return(value);
}



void
stuffedBits::setBit(bool bit) {
  ensureSpaceInCurrentBlock(1);

  uint64  mask = ((uint64)1 << (_dataBit - 1));

  if (bit)
    _data[_dataWrd] |=  mask;
  else
    _data[_dataWrd] &= ~mask;

  _dataPos++;
  _dataBit--;

  updateBit();
  updateLen();
}

}  //  namespace merylutil::bits::v1
