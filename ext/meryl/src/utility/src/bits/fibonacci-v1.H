
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

#ifndef MERYLUTIL_BITS_FIBONACCI_V1_H
#define MERYLUTIL_BITS_FIBONACCI_V1_H

#include "types.H"

namespace merylutil::inline bits::inline v1 {

extern uint64 __fibNumV[93];

inline
uint64
fibonacciNumber(uint32 f) {
  assert(f < 93);
  return(__fibNumV[f]);
}

}  //  namespace merylutil::bits::v1

#endif  //  MERYLUTIL_BITS_FIBONACCI_V1_H
