
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

#include "kmers.H"

namespace merylutil::inline kmers::v2 {

merylHistogram::merylHistogram(uint32 size) {
  _numUnique     = 0;
  _numDistinct   = 0;
  _numTotal      = 0;

  _histMax       = size;
   _histSml       = new uint64 [_histMax];

  for (uint64 ii=0; ii<_histMax; ii++)
    _histSml[ii] = 0;
}


merylHistogram::~merylHistogram() {
  delete [] _histSml;
}



void
merylHistogram::clear(void) {
  _numUnique     = 0;
  _numDistinct   = 0;
  _numTotal      = 0;

  for (uint64 ii=0; ii<_histMax; ii++)
    _histSml[ii] = 0;

  _histBig.clear();
}



void
merylHistogram::dump(stuffedBits *bits) {

  //  This only writes the latest version.

  bits->setBinary(64, _numUnique);
  bits->setBinary(64, _numDistinct);
  bits->setBinary(64, _numTotal);

  //  Find the maximum value, and count how many values we have
  //  in the histogram.

  uint64   numValues = _histBig.size();

  for (uint32 ii=0; ii<_histMax; ii++)
    if (_histSml[ii] > 0)
      numValues++;

  bits->setBinary(64, numValues);

  //  Now the data!

  for (uint32 ii=0; ii<_histMax; ii++) {
    if (_histSml[ii] > 0) {
      bits->setBinary(64,          ii);     //  Value
      bits->setBinary(64, _histSml[ii]);    //  Number of occurrences
    }
  }

  for (auto it=_histBig.begin(); it != _histBig.end(); it++) {
    bits->setBinary(64, it->first);      //  Value
    bits->setBinary(64, it->second);     //  Number of occurrences
  }
}



void
merylHistogram::dump(FILE        *outFile) {
  stuffedBits  *bits = new stuffedBits;

  dump(bits);

  bits->dumpToFile(outFile);

  delete bits;
}



void
merylHistogram::load_v01(stuffedBits *bits) {

  _numUnique   = bits->getBinary(64);
  _numDistinct = bits->getBinary(64);
  _numTotal    = bits->getBinary(64);

  _histMax     = bits->getBinary(32);   //  Number of small values
  uint64 big   = bits->getBinary(32);   //  Number of big values

  //  Versions 1 and 2 failed to store the big histogram values.

  delete [] _histSml;

  _histSml = bits->getBinary(64, _histMax);
}



void
merylHistogram::load_v03(stuffedBits *bits) {

  _numUnique   = bits->getBinary(64);
  _numDistinct = bits->getBinary(64);
  _numTotal    = bits->getBinary(64);

  uint64 hl    = bits->getBinary(64);

  //  Version 3 stores the histogram as 'hl' pairs of value,occurrence.  To
  //  load, we read each pair, then insert into either the small or big
  //  space.

  for (uint64 ii=0; ii<hl; ii++) {
    uint64  v = bits->getBinary(64);
    uint64  o = bits->getBinary(64);

    if (v < _histMax)
      _histSml[v] = o;
    else
      _histBig[v] = o;
  }
}



void
merylHistogram::load(stuffedBits *bits,
                     uint32       version) {

  switch (version) {
    case 1:
    case 2:
      load_v01(bits);
      break;
    case 3:
      load_v03(bits);
      break;
    default:
      fprintf(stderr, "merylHistogram::load()-- Unknown version %u\n", version), exit(1);
      break;
  }
}



void
merylHistogram::load(FILE        *inFile,
                     uint32       version) {
  stuffedBits  *bits = new stuffedBits;

  bits->loadFromFile(inFile);

  load(bits, version);

  delete bits;
}



//  The insert() function does not need to explicitly add in _numUnique,
//  _numDistinct or _numTotal - they're implicilty updated by addValue().
//
void
merylHistogram::insert(merylHistogram *that) {

  if (that == nullptr)
    return;

  for (uint64 val=0; val<that->_histMax; val++)
    addValue(val, that->_histSml[val]);

  for (auto it=that->_histBig.begin(); it != that->_histBig.end(); it++)
    addValue(it->first, it->second);
}



void
merylHistogram::reportHistogram(FILE *F) {

  for (uint64 val=0; val<_histMax; val++)
    if (_histSml[val] > 0)
      fprintf(F, F_U64 "\t" F_U64 "\n", val, _histSml[val]);

  for (auto it=_histBig.begin(); it != _histBig.end(); it++)
    fprintf(F, F_U64 "\t" F_U64 "\n", it->first, it->second);
}



void
merylHistogram::reportStatistics(FILE *F) {

  uint64  nUniverse = buildLowBitMask<uint64>(kmer::merSize() * 2) + 1;
  uint64  sDistinct = 0;
  uint64  sTotal    = 0;

  fprintf(F, "Number of %u-mers that are:\n", kmer::merSize());
  fprintf(F, "  unique   %20" F_U64P "  (exactly one instance of the kmer is in the input)\n", _numUnique);
  fprintf(F, "  distinct %20" F_U64P "  (non-redundant kmer sequences in the input)\n", _numDistinct);
  fprintf(F, "  present  %20" F_U64P "  (...)\n", _numTotal);
  fprintf(F, "  missing  %20" F_U64P "  (non-redundant kmer sequences not in the input)\n", nUniverse - _numDistinct);
  fprintf(F, "\n");
  fprintf(F, "             number of   cumulative   cumulative     presence\n");
  fprintf(F, "              distinct     fraction     fraction   in dataset\n");
  fprintf(F, "frequency        kmers     distinct        total       (1e-6)\n");
  fprintf(F, "--------- ------------ ------------ ------------ ------------\n");

  auto emitLine = [&] (uint64 value, uint64 occur) {
                    sDistinct  += occur;
                    sTotal     += occur * value;

                    fprintf(F, "%9" F_U64P " %12" F_U64P " %12.4f %12.4f %12.6f\n",
                            value,
                            occur,
                            (double)sDistinct / _numDistinct,
                            (double)sTotal    / _numTotal,
                            (double)value     / _numTotal * 1000000.0);
                  };

  for (uint64 val=0; val<_histMax; val++)
    if (_histSml[val] > 0)
      emitLine(val, _histSml[val]);

  for (auto it=_histBig.begin(); it != _histBig.end(); it++)
    emitLine(it->first, it->second);
}



void
merylHistogramIterator::construct(merylHistogram &that) {
  uint32  nV = that._histBig.size();

  for (uint64 val=0; val<that._histMax; val++)
    if (that._histSml[val] > 0)
      nV++;

  _val = new uint64 [nV];
  _occ = new uint64 [nV];

  for (uint64 val=0; val<that._histMax; val++) {
    if (that._histSml[val] > 0) {
      _val[_len] = val;
      _occ[_len] = that._histSml[val];
      _len++;
    }
  }

  for (auto it=that._histBig.begin(); it != that._histBig.end(); it++) {
    _val[_len] = it->first;
    _occ[_len] = it->second;
    _len++;
  }

  assert(_len == nV);
}

}  //  namespace merylutil::kmers::v2
