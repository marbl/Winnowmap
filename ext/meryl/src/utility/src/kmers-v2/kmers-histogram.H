
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

#ifndef MERYLUTIL_KMERS_HISTOGRAM_V2_H
#define MERYLUTIL_KMERS_HISTOGRAM_V2_H

#ifndef MERYLUTIL_KMERS_H
#error "include kmers.H, not this."
#endif

#include <map>

#include "types.H"
#include "bits.H"

namespace merylutil::inline kmers::v2 {

//  Stores a histogram of kmer count values.
//
//  The representation allows updates at any time (v1, v2 and v3 did not)
//  but has a slight penalty of using a map<> for large values.

class merylHistogram {
public:
  merylHistogram(uint32 size = 32 * 1048576);
  ~merylHistogram();

  void      addValue(kmvalu value, uint64 occur=1);

  void      clear(void);

  void      dump(stuffedBits *bits);
  void      dump(FILE        *outFile);

  void      load_v01(stuffedBits *bits);
  void      load_v03(stuffedBits *bits);
  void      load(stuffedBits *bits,    uint32 version);
  void      load(FILE        *inFile,  uint32 version);

  void      insert(merylHistogram *that);
  void      insert(merylHistogram &that)          { insert(&that);        };

  uint64    numUnique(void)                       { return(_numUnique);   };
  uint64    numDistinct(void)                     { return(_numDistinct); };
  uint64    numTotal(void)                        { return(_numTotal);    };

  void      reportHistogram(FILE *F);
  void      reportStatistics(FILE *F);

private:
  uint64                   _numUnique;
  uint64                   _numDistinct;
  uint64                   _numTotal;

  uint32                   _histMax;    //  Max value that can be stored in _histSml.
  uint64                  *_histSml;    //  Values smaller than _histMax.
  std::map<uint64, uint64> _histBig;    //  Values bigger than _histMax; <value,occurrances>

  friend class merylHistogramIterator;
};


//  Converts a merylHistogram into something that allows iteration over the
//  histogram and excludes the values with no occurrences.

class merylHistogramIterator {
public:
  merylHistogramIterator(merylHistogram &that)   { construct( that); }
  merylHistogramIterator(merylHistogram *that)   { construct(*that); }

  ~merylHistogramIterator() {
    delete [] _val;
    delete [] _occ;
  }

  uint64    maxValue(void) {
    if (_len == 0)
      return(0);
    else
      return(_val[_len-1]);
  };

  uint32    histogramLength(void)                 { return(_len);         };
  uint64    histogramValue(uint32 i)              { return(_val[i]);      };
  uint64    histogramOccurrences(uint32 i)        { return(_occ[i]);      };

private:
  void      construct(merylHistogram &that);

private:
  uint64   _len = 0;

  uint64  *_val = nullptr;
  uint64  *_occ = nullptr;
};



inline
void
merylHistogram::addValue(kmvalu value, uint64 occur) {

  if (value == 0)
    return;

  if (value == 1)
    _numUnique += occur;

  _numDistinct += occur;
  _numTotal    += occur * value;

  if (value < _histMax)
    _histSml[value] += occur;
  else
    _histBig[value] += occur;
}

}  //  namespace merylutil::kmers::v2

#endif  //  MERYLUTIL_KMERS_HISTOGRAM_V2_H
