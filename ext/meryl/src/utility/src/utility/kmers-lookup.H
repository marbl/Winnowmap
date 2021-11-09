
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

#ifndef MERYL_UTIL_KMER_LOOKUP_H
#define MERYL_UTIL_KMER_LOOKUP_H

#ifndef MERYL_UTIL_KMER_H
#error "include kmers.H, not this."
#endif


class merylExactLookupProgress {
public:
  merylExactLookupProgress(char const *format, bool enabled=false) {
    if (enabled)
      strncpy(_formatString, format, 128);
  }

public:
  void     tick(uint32 ff) {
    _progressCounter[ff]++;
  }

  void     show(void) {
    if (_formatString[0] == 0)
      return;

    for (uint32 bb=0; bb<64;bb++) {
      if      (_progressCounter[bb] == 0)
        _progressString[bb] = '.';
      else if (_progressCounter[bb] == uint64max)
        _progressString[bb] = '*';
      else
        _progressString[bb] = (char [4]){ '|', '/', '-', '\\' } [ _progressCounter[bb] & 0x03 ];
    }

    _progressString[64] = 0;

    fprintf(stderr, _formatString, _progressString);
  }

  void     stop(uint32 ff) {
    _progressCounter[ff] = uint64max;

    for (uint32 ii=0; ii<64; ii++)
      if (_progressCounter[ii] != uint64max)
        return;

    if (_formatString[0] == 0) {
      show();
      fprintf(stderr, "\n");
    }
  }

private:
  uint64  _progressCounter[64] = {0};
  char    _progressString[65]  = {0};
  char    _formatString[128]   = {0};
};



class merylExactLookup {
public:
  merylExactLookup() {
  };
  ~merylExactLookup() {
    delete    _sufPointer;
    delete    _sufData;
    delete    _valData;
  };

public:
  //  Optional.  Quickly analyze the input kmers and compute the minimum
  //  memory needed for the lookup tables.  The size returned is in GB.
  //
  //  maxMemInGB is used as an upper limit.
  //
  double   estimateMemoryUsage(merylFileReader *input_,
                               double           maxMemInGB_,
                               uint32           prefixSize_,
                               kmvalu           minValue_ = 0,
                               kmvalu           maxValue_ = kmvalumax) {
    initialize(input_, minValue_, maxValue_);
    return(configure(maxMemInGB_, prefixSize_, true, true));
  };

public:
  //  Load a new meryl database into the lookup table.
  //
  //  maxMemInGB is used as an upper limit on the size of the lookup table.
  //  The actual size used is determined from useMinimalMemory or
  //  useOptimalMemory, as returned from estinmateMemoryUsage().
  //
  //  The difference between 'minimal' and 'optimal' is one of speed; lookups
  //  with 'minimal' memory will be slower than with 'optimal' memory;
  //  however, it isn't known how significant this is.
  //
  //  The return value is the actual memory used, in GB, or -1 if loading
  //  failed.
  //
  double   load(merylFileReader *input_,
                double           maxMemInGB_,
                uint32           prefixSize,
                kmvalu           minValue_      = 0,
                kmvalu           maxValue_      = kmvalumax);

public:
  //  For describing what we've loaded.
  //
  uint64   nKmers(void)  {  return(_nKmersLoaded);  };

  //  The accessors.
  //
  //  Return true/false if the kmer exists/does not.
  //  Return true/false if the kmer exists/does not, and populate 'value' with the value.
  //  Return the value of the kmer, or zero if it doesn't exist.
  //
  bool     exists(kmer k);
  bool     exists(kmer k, kmvalu &value);
  kmvalu   value(kmer k);

  //  For testing the implementation.
  //
  bool     exists_test(kmer k);

private:
  //  Used internally for construction.  As tempting is it seems to call
  //  initialize() or configure() directly, you can't.
  //
  void     initialize(merylFileReader *input_, kmvalu minValue_, kmvalu maxValue_);
  void     computeSpace(bool reportMemory,
                        bool reportSizes,
                        bool compute,
                        uint32 &pbMin, uint64 &minSpace);
  double   configure(double  memInGB,
                     uint32  prefixSize,
                     bool    reportMemory,
                     bool    reportSizes);
  void     count(void);
  double   allocate(void);
  void     load(void);

  kmvalu   value_value(kmvalu value);

private:
  merylFileReader  *_input         = nullptr;

  uint64            _maxMemory     = 0;
  bool              _verbose       = true;

  kmvalu            _minValue      = 0;    //  Minimum value stored in the table -| both of these filter the
  kmvalu            _maxValue      = 0;    //  Maximum value stored in the table -| input kmers.
  kmvalu            _valueOffset   = 0;    //  Offset of values stored in the table.

  uint64            _nKmersLoaded  = 0;
  uint64            _nKmersTooLow  = 0;
  uint64            _nKmersTooHigh = 0;

  uint32            _Kbits         = uint32max;

  uint32            _prefixBits    = 0;    //  How many high-end bits of the kmer are an index into _sufPointer.
  uint32            _suffixBits    = 0;    //  How many low-end bits of the kmer are in the suffix table.
  uint32            _valueBits     = 0;    //  How many bits of the suffix entry are data.

  kmdata            _suffixMask    = 0;
  kmvalu            _valueMask     = 0;

  uint64            _nPrefix       = 0;    //  How many entries in _sufPointer:  2 ^ _prefixBits.
  uint64            _nSuffix       = 0;    //  How many entries in _suffixData:  nDistinct in the input database.

  uint32            _blkPtrBits    = 0;    //  How many bits wide is the begin field.
  uint64            _blkPtrMask    = 0;
  uint32            _blkLenBits    = 0;    //  How many bits wide is the length field.
  uint64            _blkLenMask    = 0;

  void               setPointers(uint64 prefix, uint64  bgn, uint64  len);
  void               getPointers(uint64 prefix, uint64 &bgn, uint64 &end);

  wordArray        *_sufPointer= nullptr;  //  Pointer to start and length of each block.

  wordArray        *_sufData   = nullptr;  //  Finally, kmer suffix data!
  wordArray        *_valData   = nullptr;  //  Finally, value data!
};





inline
void
merylExactLookup::setPointers(uint64 prefix, uint64  bgn, uint64  len) {
  _sufPointer->set(prefix, (bgn << _blkLenBits) | len);
}


inline
void
merylExactLookup::getPointers(uint64 prefix, uint64 &bgn, uint64 &end) {
  uint64  dat = _sufPointer->get(prefix);

  bgn = (dat >> _blkLenBits);
  end = (dat >> _blkLenBits) + (dat & _blkLenMask);
}


inline
kmvalu
merylExactLookup::value_value(kmvalu value) {
  if (_valueBits == 0)               //  Return 'true' if no value
    return(1);                       //  is stored.

  value &= _valueMask;

  return(value + _valueOffset);      //  Otherwise, return the value.
};



//  Return true/false if the kmer exists/does not.
inline
bool
merylExactLookup::exists(kmer k) {
  kmdata  kmer   = (kmdata)k;
  uint64  prefix = kmer >> _suffixBits;
  kmdata  suffix = kmer  & _suffixMask;
  kmdata  tag;
  uint64  bgn, mid, end;

  getPointers(prefix, bgn, end);

  //  Binary search for the matching tag.

  while (bgn + 8 < end) {
    mid = bgn + (end - bgn) / 2;

    tag = _sufData->get(mid);

    if (tag == suffix)
      return(true);

    if (suffix < tag)
      end = mid;

    else
      bgn = mid + 1;
  }

  //  Switch to linear search when we're down to just a few candidates.

  for (mid=bgn; mid < end; mid++) {
    tag = _sufData->get(mid);

    if (tag == suffix)
      return(true);
  }

  return(false);
}



//  Return true/false if the kmer exists/does not.
//  And populate 'value' with the value of the kmer.
inline
bool
merylExactLookup::exists(kmer k, kmvalu &value) {
  kmdata  kmer   = (kmdata)k;
  kmdata  prefix = kmer >> _suffixBits;
  kmdata  suffix = kmer  & _suffixMask;
  kmdata  tag;
  uint64  bgn, mid, end;

  getPointers(prefix, bgn, end);

  //  Binary search for the matching tag.

  while (bgn + 8 < end) {
    mid = bgn + (end - bgn) / 2;

    tag = _sufData->get(mid);

    if (tag == suffix) {
      if (_valueBits == 0)
        value = 1;
      else
        value = _valData->get(mid);
      return(true);
    }

    if (suffix < tag)
      end = mid;

    else
      bgn = mid + 1;
  }

  //  Switch to linear search when we're down to just a few candidates.

  for (mid=bgn; mid < end; mid++) {
    tag = _sufData->get(mid);

    if (tag == suffix) {
      if (_valueBits == 0)
        value = 1;
      else
        value = _valData->get(mid);
      return(true);
    }
  }

  value = 0;
  return(false);
}


//  Returns the value of the kmer, '0' if it doesn't exist.
inline
kmvalu
merylExactLookup::value(kmer k) {
  kmdata  kmer   = (kmdata)k;
  kmdata  prefix = kmer >> _suffixBits;
  kmdata  suffix = kmer  & _suffixMask;
  kmdata  tag;
  uint64  bgn, mid, end;

  getPointers(prefix, bgn, end);

  //  Binary search for the matching tag.

  while (bgn + 8 < end) {
    mid = bgn + (end - bgn) / 2;

    tag = _sufData->get(mid);

    if (tag == suffix) {
      if (_valueBits == 0)
        return(1);
      else
        return(_valData->get(mid));
    }

    if (suffix < tag)
      end = mid;

    else
      bgn = mid + 1;
  }

  //  Switch to linear search when we're down to just a few candidates.

  for (mid=bgn; mid < end; mid++) {
    tag = _sufData->get(mid);

    if (tag == suffix) {
      if (_valueBits == 0)
        return(1);
      else
        return(_valData->get(mid));
    }
  }

  return(0);
};


#endif  //  MERYL_UTIL_KMER_LOOKUP_H
