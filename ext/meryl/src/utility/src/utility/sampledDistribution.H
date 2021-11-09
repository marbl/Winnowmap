
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

#ifndef SAMPLED_DISTRIBUTION_H
#define SAMPLED_DISTRIBUTION_H

#include "types.H"

class sampledDistribution {
public:
  sampledDistribution();
  ~sampledDistribution();

  void     loadDistribution(char const *path);
  bool     empty(void) {
    return(_dataLen == 0);
  };

  //  Imagine our input histogram (value occurences) is expanded into an
  //  array of (_dataSum) values, where each 'value' in the input histogram
  //  is listed 'occurences' times (and that it's sorted).
  //
  //  For input 0 <= 'd' < 1, we want to return the 'value' that is at that
  //  position in the array.
  //
  //  Scan the _data, looking for the _data element that contains the 'lim'th
  //  entry in the (imagined) array.
  //
  uint32   getValue(double d) {
    if (d < 0.0)  d = 0.0;          //  Limit our input (random) number
    if (d > 1.0)  d = 1.0;          //  to valid scaling range.

    uint64  lim = (uint64)floor(_dataSum * d);
    uint64  val = 0;

    while (_data[val] <= lim) {     //  If _data[value] is more than the current
      lim -= _data[val];            //  limit, we're found the value, otherwise,
      val += 1;                     //  decrement the limit by the occurrences
    }                               //  for this value and move to the next.

    assert(val < _dataLen);

    return(val);
  };


public:
  uint32   _dataLen = 0;         //  Length of the valid data in _data: _data[_dataLen-1] is the last valid data point.
  uint32   _dataMax = 0;         //  Allocated size of _data.
  uint32  *_data    = nullptr;   //  Number of data points of value x:  _data[x] == v <-> 'x' was present 'v' times in the input.

  uint64   _dataSum = 0;         //  Number of points in the input.  It's the sum of all _data[x].
};


#endif  //  SAMPLED_DISTRIBUTION_H
