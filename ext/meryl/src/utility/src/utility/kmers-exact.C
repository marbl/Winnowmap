
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

#include <vector>
#include <algorithm>


//  Set some basic boring stuff.
//
void
merylExactLookup::initialize(merylFileReader *input_, kmvalu minValue_, kmvalu maxValue_) {

  //  Save a pointer to the input data.

  _input = input_;

  //  Silently make minValue and maxValue be valid values.

  if (minValue_ == 0)
    minValue_ = 1;

  if (maxValue_ == kmvalumax) {
    uint32  nV = _input->stats()->histogramLength();

    maxValue_ = _input->stats()->histogramValue(nV - 1);
  }

  //  Now initialize filtering!

  _minValue       = minValue_;
  _maxValue       = maxValue_;
  _valueOffset    = minValue_ - 1;                   //  "1" stored in the data is really "minValue" to the user.

  _nKmersLoaded   = 0;
  _nKmersTooLow   = 0;
  _nKmersTooHigh  = 0;

  //  Now initialize table parameters!

  _Kbits          = kmer::merSize() * 2;

  _prefixBits     = 0;                               //  Bits of the kmer used as an index into the table.
  _suffixBits     = 0;                               //  Width of an entry in the suffix table.
  _valueBits      = 0;                               //  (also in the suffix table)

  if (_maxValue >= _minValue)
    _valueBits = countNumberOfBits64(_maxValue + 1 - _minValue);

  _suffixMask     = 0;

  _nPrefix        = 0;                               //  Number of entries in pointer table.
  _nSuffix        = 0;                               //  Number of entries in suffix dable.

  //  Scan the histogram to count the number of kmers in range.

  for (uint32 ii=0; ii<_input->stats()->histogramLength(); ii++) {
    kmvalu  v = _input->stats()->histogramValue(ii);

    if ((_minValue <= v) &&
        (v <= _maxValue))
      _nSuffix += _input->stats()->histogramOccurrences(ii);
  }

  _prePtrBits     = countNumberOfBits64(_nSuffix);   //  Width of an entry in the prefix table.
  _prePtrBits     = 64;

  _suffixBgn      = nullptr;
  _suffixLen      = nullptr;
  _suffixEnd      = nullptr;
  _sufData        = nullptr;
  _valData        = nullptr;
}



//  Analyze the number of kmers to store in the table, to decide on
//  various parameters for allocating the table - how many bits to
//  use for indexing (prefixSize), and how many bits of data we need
//  to store explicitly (suffixBits and valueBits).
//
void
merylExactLookup::configure(double  memInGB,
                            double &memInGBmin,
                            double &memInGBmax,
                            bool    useMinimalMemory,
                            bool    useOptimalMemory,
                            bool    reportMemory,
                            bool    reportSizes) {

  //  Convert the memory in GB to memory in BITS.  If no memory
  //  size is supplied, as the OS how big we can get.

  if (memInGB == 0.0)
    _maxMemory = getMaxMemoryAllowed() * 8;
  else
    _maxMemory = (uint64)(memInGB * 1024.0 * 1024.0 * 1024.0 * 8);

  //  Find the prefixBits that results in the smallest allocated memory size.
  //  Due to threading over the files, we cannot use a prefix smaller than 6
  //  bits.
  //
  //  While it's nice to find the smallest memory size possible, that's also
  //  about the slowest possible.  Instead, empirically determined on a small
  //  test, allow a very sparse table of 16 to 32 prefixes per kmer (if
  //  possible).

  uint64  minSpace   = uint64max;
  uint64  optSpace   = uint64max;
  uint64  usdSpace   = uint64max;

  //  _nSuffix here is just the number of distinct kmers in the input.  We'll
  //  search for prefix sizes up to that size plus a bit more to show that
  //  what we pick really is the best size.
  //
  //  We save the smallest size, and the 'optimal' size, defined as something
  //  at least as big as the smallest, but not more than 8 times larger.

  uint32  pbMin      = 0;
  uint32  pbOpt      = 0;
  uint32  pbMax      = countNumberOfBits64(_nSuffix) + 1;

  if (pbMax > kmer::merSize() * 2)
    pbMax = kmer::merSize() * 2;

  for (uint32 pb=0; pb<pbMax; pb++) {
    uint64  nprefix = (uint64)1 << pb;
    uint64  space   = nprefix * _prePtrBits + _nSuffix * (_Kbits - pb) + _nSuffix * _valueBits;

    if (space < minSpace) {
      pbMin        = pb;
      minSpace     = space;
    }

    if ((space < _maxMemory) && (pb < pbMin + 4)) {
      pbOpt        = pb;
      optSpace     = space;
    }
  }

  //  Set parameters.  For logging, we need these set even if
  //  useMinimalMemory and useOptimalMemory are false -- this happens when
  //  we're called from estimateMemoryUsage.

  if (useMinimalMemory == true) {
    usdSpace = minSpace;

    _prefixBits  =          pbMin;
    _suffixBits  = _Kbits - pbMin;

    _suffixMask  = buildLowBitMask<kmdata>(_suffixBits);

    _nPrefix     = (uint64)1 << pbMin;
  }

  if (useOptimalMemory == true) {
    usdSpace = optSpace;

    _prefixBits  =          pbOpt;
    _suffixBits  = _Kbits - pbOpt;

    _suffixMask  = buildLowBitMask<kmdata>(_suffixBits);

    _nPrefix     = (uint64)1 << pbOpt;
  }

  //  And do it all again to keep the users entertained.

  if (reportMemory) {
    fprintf(stderr, "\n");
    fprintf(stderr, " p       prefixes             bits gigabytes (allowed: %lu GB)\n", _maxMemory >> 33);
    fprintf(stderr, "-- -------------- ---------------- ---------\n");

    uint32  minpb = (pbMin < 4)          ? 1      : pbMin - 4;  //  Show four values before and
    uint32  maxpb = (_Kbits < pbOpt + 5) ? _Kbits : pbOpt + 5;  //  four after the smallest.

    if (pbOpt == 0)
      maxpb = minpb + 10;

    for (uint32 pb=minpb; pb < maxpb; pb++) {
      uint64  nprefix = (uint64)1 << pb;
      uint64  space   = nprefix * _prePtrBits + _nSuffix * (_Kbits - pb) + _nSuffix * _valueBits;

      if     ((pb == pbMin) &&
              (pb == pbOpt))
        fprintf(stderr, "%2u %14lu %16lu %9.3f (smallest)\n", pb, nprefix, space, bitsToGB(space));
      else if (pb == pbMin)
        fprintf(stderr, "%2u %14lu %16lu %9.3f (smallest)\n", pb, nprefix, space, bitsToGB(space));
      else if (pb == pbOpt)
        fprintf(stderr, "%2u %14lu %16lu %9.3f (faster)\n",   pb, nprefix, space, bitsToGB(space));
      else
        fprintf(stderr, "%2u %14lu %16lu %9.3f\n",            pb, nprefix, space, bitsToGB(space));
    }

    fprintf(stderr, "-- -------------- ---------------- ---------\n");
    fprintf(stderr, "   %14lu total kmers\n", _nSuffix);
    fprintf(stderr, "\n");
  }

  if (reportSizes) {
    fprintf(stderr, "\n");
    fprintf(stderr, "For %lu distinct %u-mers (with %u bits used for indexing and %u bits for tags):\n", _nSuffix, _Kbits / 2, _prefixBits, _suffixBits);
    fprintf(stderr, "  %7.3f GB memory for kmer indices - %12lu elements %2u bits wide)\n", bitsToGB(_nPrefix * _prePtrBits), _nPrefix, _prePtrBits);
    fprintf(stderr, "  %7.3f GB memory for kmer tags    - %12lu elements %2u bits wide)\n", bitsToGB(_nSuffix * _suffixBits), _nSuffix, _suffixBits);
    fprintf(stderr, "  %7.3f GB memory for kmer values  - %12lu elements %2u bits wide)\n", bitsToGB(_nSuffix * _valueBits),  _nSuffix, _valueBits);
    fprintf(stderr, "  %7.3f GB memory\n",                                                  bitsToGB(usdSpace));
    fprintf(stderr, "\n");
  }

  //  Copy the min and optimal memory sizes to the output variables.

  memInGBmin = bitsToGB(minSpace);
  memInGBmax = bitsToGB(optSpace);
}



//  Make one pass through the file to count how many kmers per prefix we will end
//  up with.  This is needed only if kmers are filtered, but does
//  make the rest of the loading a little easier.
//
//  The loop control and kmer loading is the same in the two loops.
void
merylExactLookup::count(void) {

  _suffixBgn = new uint64 [_nPrefix];
  _suffixLen = new uint64 [_nPrefix];
  _suffixEnd = new uint64 [_nPrefix];

  for (uint64 ii=0; ii<_nPrefix; ii++)
    _suffixBgn[ii] = _suffixLen[ii] = _suffixEnd[ii] = uint64zero;

  //  Scan all kmer files, counting the number of kmers per prefix.
  //  This is thread safe when _prefixBits is more than 6 (the number of files).

  uint32   nf = _input->numFiles();

  assert(nf == 64);

  uint64   minp[nf];
  uint64   maxp[nf];

  for (uint32 ii=0; ii<nf; ii++) {
    minp[ii] = uint64max;
    maxp[ii] = uint64min;
  }

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<nf; ff++) {
    FILE                  *blockFile = _input->blockFile(ff);
    merylFileBlockReader  *block     = new merylFileBlockReader;

    //  Keep local counters, otherwise, we collide when updating the global counts.

    uint64  tooLow  = 0;
    uint64  tooHigh = 0;
    uint64  loaded  = 0;

    //  Load blocks until there are no more.

    while (block->loadBlock(blockFile, ff) == true) {
      block->decodeBlock();

      for (uint32 ss=0; ss<block->nKmers(); ss++) {
        kmdata   kbits  = 0;
        kmdata   prefix = 0;
        kmvalu   value  = block->values()[ss];

        if (value < _minValue) {
          tooLow++;
          continue;
        }

        if (_maxValue < value) {
          tooHigh++;
          continue;
        }

        loaded++;

        kbits   = block->prefix();         //  Combine the file prefix and
        kbits <<= _input->suffixSize();    //  suffix data to reconstruct
        kbits  |= block->suffixes()[ss];   //  the kmer bits.

        prefix = kbits >> _suffixBits;     //  Then extract the prefix

        minp[ff] = std::min(minp[ff], (uint64)prefix);
        maxp[ff] = std::max(maxp[ff], (uint64)prefix);

        assert(prefix < _nPrefix);

        _suffixLen[prefix]++;              //  Count the number of kmers per prefix.
      }
    }

#pragma omp critical (count_stats)
    {
      _nKmersTooLow  += tooLow;
      _nKmersTooHigh += tooHigh;
      _nKmersLoaded  += loaded;
    }

    delete block;

    AS_UTL_closeFile(blockFile);
  }

  //  If the min/max intersect, we've got a problem somewhere.  Each 'prefix'
  //  will map to exactly one file, and they're supposed to map
  //  consecutively.  Good luck figuring out what broke if this triggers.

  for (uint32 ii=1; ii<nf; ii++)
    assert(maxp[ii-1] < minp[ii]);

  //  Now that we know the length of each block, we can set _suffixBgn to the
  //  address of the first element.  _suffixEnd is set to that too; we'll use
  //  it to load data into the table.
  //
  //  To allow threads without locks, we need to pad the end of each block so
  //  that two blocks don't share a wordArray word.  Instead, we just pad the
  //  last block that each thread (one thread per file) will access.  A
  //  little bit harder to figure out, but less memory used.
  //
  //  For a prefix of [ffffffppp..pppp] a single thread will process all
  //  kmers [ffffff......].  Thus, when the prefix ends in 1111...111, we can
  //  just bump up 'bgn' a bit, just enough to get to the next 128-bit word.
  //  But since this index is used both in storing suffixes and values,
  //  that's impossible and we just add 256 bits.

  uint64 mask = (_nPrefix - 1) >> 6;

  for (uint64 bgn=0, ii=0; ii<_nPrefix; ii++) {
    _suffixBgn[ii] = bgn;
    _suffixEnd[ii] = bgn;

    bgn += _suffixLen[ii];

    if ((ii & mask) == mask)
      bgn += 256;
  }

  //  Log.

  if (_verbose)
    fprintf(stderr, "Will load " F_U64 " kmers.  Skipping " F_U64 " (too low) and " F_U64 " (too high) kmers.\n",
            _nKmersLoaded, _nKmersTooLow, _nKmersTooHigh);
}



//  With all parameters known, just grab and clear memory.
//
//  The block size used in the wordArray _sufData is chosen so that large
//  arrays have not-that-many allocations.  The array is pre-allocated, to
//  prevent the need for any locking or coordination when filling out the
//  array.
//
double
merylExactLookup::allocate(void) {
  uint64  arraySize;
  uint64  arrayBlockMin;
  double  memInGBused = 0.0;

  uint64  ns = _suffixEnd[_nPrefix-1];   //  The largest word we access in wordArray.

  if (_suffixBits > 0) {
    arraySize      = ns * _suffixBits;
    arrayBlockMin  = std::max(arraySize / 1024llu, 268435456llu);   //  In bits, so 32MB per block.
    memInGBused   += bitsToGB(arraySize);

    if (_verbose)
      fprintf(stderr, "Allocating space for %lu suffixes of %u bits each -> %lu bits (%.3f GB) in blocks of %.3f MB\n",
              ns, _suffixBits, arraySize, bitsToGB(arraySize), bitsToMB(arrayBlockMin));

    assert(_suffixBits <= 128);

    _sufData = new wordArray(_suffixBits, arrayBlockMin, false);
    _sufData->allocate(ns);
  }

  if (_valueBits > 0) {
    arraySize     = ns * _valueBits;
    arrayBlockMin = std::max(arraySize / 1024llu, 268435456llu);   //  In bits, so 32MB per block.
    memInGBused   += bitsToGB(arraySize);

    if (_verbose)
      fprintf(stderr, "                     %lu values   of %u bits each -> %lu bits (%.3f GB) in blocks of %.3f MB\n",
              ns, _valueBits,  arraySize, bitsToGB(arraySize), bitsToMB(arrayBlockMin));

    assert(_valueBits <= 64);

    _valData = new wordArray(_valueBits, arrayBlockMin, false);
    _valData->allocate(ns);
  }

  return(memInGBused);
}



//  Each file can be processed independently IF we know how many kmers are in
//  each prefix.  For that, we need to load the merylFileReader index.
//  We don't, actually, know that if we're filtering out low/high count kmers.
//  In this case, we overallocate, but cannot cleanup at the end.
void
merylExactLookup::load(void) {
  uint32   nf      = _input->numFiles();
  uint64   sufMask = buildLowBitMask<kmdata>(_suffixBits);
  uint64   valMask = buildLowBitMask<kmvalu>(_valueBits);

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<nf; ff++) {
    FILE                  *blockFile = _input->blockFile(ff);
    merylFileBlockReader  *block     = new merylFileBlockReader;

    //  Load blocks until there are no more.

    while (block->loadBlock(blockFile, ff) == true) {
      block->decodeBlock();

      for (uint32 ss=0; ss<block->nKmers(); ss++) {
        kmdata   kbits  = 0;
        kmdata   prefix = 0;
        kmdata   suffix = 0;
        kmvalu   value  = block->values()[ss];

        if ((value < _minValue) ||         //  Sanity checking and counting done
            (_maxValue < value))           //  in count() above.
          continue;

        kbits   = block->prefix();         //  Combine the file prefix and
        kbits <<= _input->suffixSize();    //  suffix data to reconstruct
        kbits  |= block->suffixes()[ss];   //  the kmer bits.

        suffix = kbits  & sufMask;         //  Then extract the prefix
        prefix = kbits >> _suffixBits;     //  and suffix to use in the table

        _sufData->set(_suffixEnd[prefix], suffix);

        //  Compute and store the value, if requested.

        if (_valueBits > 0) {
          value -= _valueOffset;

          if (value > _maxValue + 1 - _minValue)
            fprintf(stderr, "minValue " F_U32 " maxValue " F_U32 " value " F_U32 " bits " F_U32 "\n",
                    _minValue, _maxValue, value, _valueBits);
          assert(value <= valMask);

          _valData->set(_suffixEnd[prefix], value);
        }

        //  Move to the next item.

        _suffixEnd[prefix]++;
      }
    }

    delete block;

    AS_UTL_closeFile(blockFile);
  }

  //  Check that we loaded the expected number of kmers into each space

  for (uint64 ii=0; ii<_nPrefix; ii++)
    assert(_suffixBgn[ii] + _suffixLen[ii] == _suffixEnd[ii]);
  
  //  Now just log.

  if (_verbose)
    fprintf(stderr, "Loaded " F_U64 " kmers.  Skipped " F_U64 " (too low) and " F_U64 " (too high) kmers.\n",
            _nKmersLoaded, _nKmersTooLow, _nKmersTooHigh);
}



void
merylExactLookup::estimateMemoryUsage(merylFileReader *input_,
                                      double           maxMemInGB_,
                                      double          &minMemInGB_,
                                      double          &optMemInGB_,
                                      kmvalu           minValue_,
                                      kmvalu           maxValue_) {
  initialize(input_, minValue_, maxValue_);
  configure(maxMemInGB_, minMemInGB_, optMemInGB_, false, false, true, false);
}



double
merylExactLookup::load(merylFileReader *input_,
                       double           maxMemInGB_,
                       bool             useMinimalMemory,
                       bool             useOptimalMemory,
                       kmvalu           minValue_,
                       kmvalu           maxValue_) {
  double  minMem  = 0.0;
  double  maxMem  = 0.0;
  double  memInGBused = 0.0;

  initialize(input_, minValue_, maxValue_);            //  Initialize ourself.

  configure(maxMemInGB_,                               //  Find parameters.
            minMem,
            maxMem,
            useMinimalMemory,
            useOptimalMemory,
            false,
            true);

  if (_prefixBits == 0)                                //  Fail if needed.
    return(0.0);

  count();                                             //  Count kmers/prefix.
  memInGBused = allocate();                            //  Allocate space.
  load();                                              //  Load data.

  return(memInGBused);
}






bool
merylExactLookup::exists_test(kmer k) {
  char    kmerString[65];
  kmdata  kmer   = (kmdata)k;
  kmdata  prefix = kmer >> _suffixBits;
  kmdata  suffix = kmer  & _suffixMask;

  fprintf(stderr, "\n");
  fprintf(stderr, "kmer        %s  %s\n", toHex(kmer, 2 * k.merSize()), k.toString(kmerString));
  fprintf(stderr, "suffixBits  %s  %3u bits\n", toHex(_suffixMask, _suffixBits), _suffixBits);
  fprintf(stderr, "prefix      %s  %3u bits\n", toHex(prefix, 2 * k.merSize() - _suffixBits), 2 * k.merSize() - _suffixBits);
  fprintf(stderr, "suffix      %s\n", toHex(suffix, _suffixBits));

  uint64  bgn = _suffixBgn[prefix];
  uint64  mid;
  uint64  end = _suffixEnd[prefix];

  kmdata  tag;

  //  Binary search for the matching tag.

  fprintf(stderr, "BINARY SEARCH the bucket %lu-%lu for suffix %s.\n", bgn, end, toHex(suffix));

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

  fprintf(stderr, "\n");
  fprintf(stderr, "FAILED kmer   0x%s\n", toHex(kmer));
  fprintf(stderr, "FAILED prefix 0x%s\n", toHex(prefix));
  fprintf(stderr, "FAILED suffix 0x%s\n", toHex(suffix));
  fprintf(stderr, "\n");
  fprintf(stderr, "original  %9lu %9lu\n", _suffixBgn[prefix], _suffixEnd[prefix]);
  fprintf(stderr, "final     %9lu %9lu\n", bgn, end);
  fprintf(stderr, "\n");

  bgn = _suffixBgn[prefix];
  end = _suffixEnd[prefix];

  fprintf(stderr, "BINARY SEARCH the bucket %lu-%lu for suffix %s.\n", bgn, end, toHex(suffix));

  while (bgn + 8 < end) {
    mid = bgn + (end - bgn) / 2;

    tag = _sufData->get(mid);

    fprintf(stderr, "TEST bgn %8lu %8lu %8lu end -- dat %s =?= %s suffix\n",
            bgn, mid, end, toHex(tag), toHex(suffix));

    if (tag == suffix)
      return(true);

    if (suffix < tag)
      end = mid;

    else
      bgn = mid + 1;
  }

  //  Exhaustively search the bucket.

  fprintf(stderr, "LINEAR SEARCH the bucket %lu-%lu for suffix %s.\n", bgn, end, toHex(suffix));

  for (mid=bgn; mid < end; mid++) {
    tag = _sufData->get(mid);

    fprintf(stderr, "ITER bgn %8lu %8lu %8lu end -- dat %s\n",
            bgn, mid, end, toHex(tag));

    if (tag == suffix)
      return(true);
  }

  //  Exhaustively search all buckets.
  //
  //  THIS IS WRONG - it needs to skip the empty buckets in the middle, so needs to
  //  iterate over each suffixBgn/suffixEnd pair individually.

  bgn = _suffixBgn[0];
  end = _suffixEnd[_nPrefix - 1];

  fprintf(stderr, "LINEAR SEARCH the entire table %lu-%lu for suffix %s.\n", bgn, end, toHex(suffix));

  for (mid=bgn; mid < end; mid++) {
    tag = _sufData->get(mid);

    fprintf(stderr, "ITER bgn %8lu %8lu %8lu end -- dat %s\n",
            bgn, mid, end, toHex(tag));

    if (tag == suffix)
      return(true);
  }

  assert(0);
};
