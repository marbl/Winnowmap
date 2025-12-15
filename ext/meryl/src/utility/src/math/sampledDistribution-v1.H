
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

#ifndef MERYLUTIL_MATH_SAMPLED_DISTRIBUTION_V1_H
#define MERYLUTIL_MATH_SAMPLED_DISTRIBUTION_V1_H

#include "types.H"

namespace merylutil::inline math::inline v1 {

class sampledDistribution {
public:
  sampledDistribution()   {                  }
  ~sampledDistribution()  { delete [] _data; }

  void     loadDistribution(char const *path);

  bool     empty(void)    { return(_dataLen == 0);  }
  uint32   getValue(double d);   //  The value of the distribution at 0 <= d <= 1.0.

public:
  uint32   _dataLen = 0;         //  Length of the valid data in _data: _data[_dataLen-1] is the last valid data point.
  uint32   _dataMax = 0;         //  Allocated size of _data.
  uint32  *_data    = nullptr;   //  Number of data points of value x:  _data[x] == v <-> 'x' was present 'v' times in the input.

  uint64   _dataSum = 0;         //  Number of points in the input.  It's the sum of all _data[x].
};

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_MATH_SAMPLED_DISTRIBUTION_V1_H
