
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

//  stuffedBits operations on Elias Gamma coded data.

namespace merylutil::inline bits::inline v1 {

uint64
stuffedBits::getEliasGamma(void) {
  uint32  N = getUnary();
  uint64  V = getBinary(N);

  V |= (uint64)1 << N;

  return(V);
}



uint64 *
stuffedBits::getEliasGamma(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getEliasGamma();

  return(values);
}



uint32
stuffedBits::setEliasGamma(uint64 value) {
  uint32 size = 0;
  uint64 N    = countNumberOfBits64(value) - 1;

  assert(value > 0);

  size += setUnary(N);
  size += setBinary(N, value);

  return(size);
}



uint32
stuffedBits::setEliasGamma(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setEliasGamma(values[ii]);

  return(size);
}

}  //  namespace merylutil::bits::v1
