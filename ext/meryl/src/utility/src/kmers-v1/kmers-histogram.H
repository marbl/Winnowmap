
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

#ifndef MERYLUTIL_KMERS_HISTOGRAM_V1_H
#define MERYLUTIL_KMERS_HISTOGRAM_V1_H

#ifndef MERYLUTIL_KMERS_H
#error "include kmers.H, not this."
#endif

#include <map>

#include "types.H"
#include "bits.H"

namespace merylutil::inline kmers::v1 {

//  Stores a histogram of kmer count values.

class merylHistogram {
public:
  merylHistogram();
  ~merylHistogram();

  void      addValue(kmvalu value) {

    if (value == 0)
      return;

    if (value == 1)
      _numUnique++;

    _numDistinct += 1;
    _numTotal    += value;

    if (value < _histMax)
      _hist[value]++;
    else
      _histBig[value]++;
  };

  void      clear(void);

  void      dump(stuffedBits *bits);
  void      dump(FILE        *outFile);

  void      load_v01(stuffedBits *bits);
  void      load_v03(stuffedBits *bits);
  void      load(stuffedBits *bits,    uint32 version);
  void      load(FILE        *inFile,  uint32 version);
  void      load(const char  *filename);

  uint64    numUnique(void)                       { return(_numUnique);   }
  uint64    numDistinct(void)                     { return(_numDistinct); }
  uint64    numTotal(void)                        { return(_numTotal);    }

  uint32    histogramLength(void)                 { return(_histLen);     }
  uint64    histogramValue(uint32 i)              { return(_histVs[i]);   }
  uint64    histogramOccurrences(uint32 i)        { return(_histOs[i]);   }

  uint64    maxValue(void)                        { return(_histVs[_histLen-1]); }

private:
  void      dumpSmoothed(uint32 n, double *h, double *s);
  void      dumpDerivs  (uint32 n, double *h);
  void      dumpMirror  (uint32 n, double *h);

public:
  void      computePloidyPeaks(bool debugPloidy=false);

  double    getNoiseTrough(void)    { return           _coveragePeaks[0];      }
  double    getCoverage(uint32 x)   { return (x < 9) ? _coveragePloidy[x] : 0; }
  double    getDepth   (uint32 x)   { return (x < 9) ? _coveragePeaks[x]  : 0; }


private:
  uint64                   _numUnique   = 0;
  uint64                   _numDistinct = 0;
  uint64                   _numTotal    = 0;

  double                   _coveragePloidy[9] = { 0 };   //  0.5x, 1x, 2x, etc.
  double                   _coveragePeaks[9]  = { 0 };   //  depth for the above.

  uint32                   _histMax = 0;         //  Max value that can be stored in _hist.
  uint64                  *_hist    = nullptr;   //
  std::map<uint64, uint64> _histBig;             //  Values bigger than _histMax; <value,occurrances>

  uint64                   _histLen = 0;         //  If loaded from disk, this is the unpacked histogram.
  //nt64                   _histMax = 0;         //
  uint64                  *_histVs  = nullptr;   //  The value this histogram entry is counting.
  uint64                  *_histOs  = nullptr;   //  The number of occurrences of that value.
};

}  //  namespace merylutil::kmers::v1

#endif  //  MERYLUTIL_KMERS_HISTOGRAM_V1_H
