
/******************************************************************************
 *
 *  This file is part of meryl, a genomic k-kmer counter with nice features.
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

#include "meryl.H"
#include "strings.H"
#include "system.H"




//  Memory used by the simple algorithm depends only on the kmer size and the count
//  of the most frequent kmer (or highest count allowed).
//
void
findExpectedSimpleSize(uint64  nKmerEstimate,
                       char   *countSuffixString,
                       uint32  countSuffixLength,
                       uint64 &memoryUsed_) {
  uint32   lowBitsSize     = sizeof(lowBits_t) * 8;
  uint64   nEntries        = (uint64)1 << (2 * kmerTiny::merSize() - 2 * countSuffixLength);

  uint64   expMaxCount     = 0.004 * nKmerEstimate;
  uint64   expMaxCountBits = countNumberOfBits64(expMaxCount) + 1;
  uint64   extraBits       = (expMaxCountBits < lowBitsSize) ? (0) : (expMaxCountBits - lowBitsSize);

  uint64   lowMem          = nEntries * lowBitsSize;
  uint64   highMem         = nEntries * extraBits;
  uint64   totMem          = (lowMem + highMem) / 8;

  memoryUsed_ = UINT64_MAX;

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "SIMPLE MODE\n");
  fprintf(stderr, "-----------\n");
  fprintf(stderr, "\n");

  if (2 * kmerTiny::merSize() - 2 * countSuffixLength > 42) {    //  21-mers need 8 TB memory for simple mode.
    fprintf(stderr, "  Not possible.\n");                        //  Just give up and say we can't do bigger.
    return;                                                      //  Of note, 31-mers (possibly 30) overflow
  }                                                              //  lowMem.

  if (countSuffixLength == 0)
    fprintf(stderr, "  %u-mers\n", kmerTiny::merSize());
  else
    fprintf(stderr, "  %u-mers with constant %u-mer suffix '%s'\n", kmerTiny::merSize(), countSuffixLength, countSuffixString);

  fprintf(stderr, "    -> %lu entries for counts up to %u.\n", nEntries, ((uint32)1 << lowBitsSize) - 1);
  fprintf(stderr, "    -> %lu %cbits memory used\n", scaledNumber(lowMem),  scaledUnit(lowMem));
  fprintf(stderr, "\n");
  fprintf(stderr, "  %lu input bases\n", nKmerEstimate);
  fprintf(stderr, "    -> expected max count of %lu, needing %lu extra bits.\n", expMaxCount, extraBits);
  if (extraBits > 0)
    fprintf(stderr, "    -> %lu %cbits memory used\n", scaledNumber(highMem), scaledUnit(highMem));
  else
    fprintf(stderr, "    -> no memory used\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  %lu %cB memory needed\n", scaledNumber(totMem),  scaledUnit(totMem));

  memoryUsed_ = totMem;
}



//  Returns bestPrefix_ and memoryUsed_ corresponding to the minimal memory
//  estimate for the supplied nKmerEstimate.
//
//  If no estimate is below memoryAllowed, 0 and UINT64_MAX, respectively,
//  are returned.
//
//  IMPORTANT!  Prefixes must be at least six bits - the number of bits
//  we use for deciding on a file - and probably need to then leave one
//  bit for a block id.  So only save a bestPrefix if it is at least 7.
//
//  IMPORTANT!  Smaller prefixes mean bigger blocks in the output file,
//  and any operation on those outputs needs to load all kmers in
//  the block into memory, as full kmers.  Only save a bestPrefix
//  if it is at least 10 - giving us at least 16 blocks per file.
//
//  IMPORTANT!  Start searching at 1 and stop at merSize-1, otherwise
//  we end up with a prefix or a suffix of size zero.
//
void
findBestPrefixSize(uint64  nKmerEstimate,
                   uint64  memoryAllowed,
                   uint32 &bestPrefix_,
                   uint64 &memoryUsed_) {
  uint32  merSize      = kmerTiny::merSize();
  uint32  segSizeBits  = merylCountArray::pagesPerSegment() * getPageSize() * 8;
  uint32  segSizeBytes = merylCountArray::pagesPerSegment() * getPageSize();

  bestPrefix_  = 0;
  memoryUsed_  = UINT64_MAX;

  for (uint32 wp=1; wp < 2 * merSize - 1; wp++) {
    uint64  nPrefix          = (uint64)1 << wp;                    //  Number of prefix == number of blocks of data
    uint64  kmersPerPrefix   = nKmerEstimate / nPrefix + 1;        //  Expected number of kmers we need to store per prefix
    uint64  kmersPerSeg      = segSizeBits / (2 * merSize - wp);   //  Kmers per segment
    uint64  segsPerPrefix    = kmersPerPrefix / kmersPerSeg + 1;   //

    if (wp + countNumberOfBits64(segsPerPrefix) + countNumberOfBits64(segSizeBytes) >= 64)
      break;   //  Otherwise, dataMemory overflows.

    uint64  structMemory     = ((sizeof(merylCountArray) * nPrefix) +                  //  Basic structs
                                (sizeof(uint64 *)        * nPrefix * segsPerPrefix));  //  Pointers to segments
    uint64  dataMemoryMin    = nPrefix *                 segSizeBytes;                 //  Minimum memory needed for this size.
    uint64  dataMemory       = nPrefix * segsPerPrefix * segSizeBytes;                 //  Expected memory for full batch.
    uint64  totalMemory      = structMemory + dataMemory;

    //  If we're too big, stop, we'll only get bigger.
    if (structMemory + dataMemoryMin > memoryAllowed)
      break;

    //  Pick a larger prefix if it is dramatically smaller than what we have.
    //  More prefixes seem to run a bit slower, but also have smaller buckets
    //  for sorting at the end.
    if ((wp > 9) && (totalMemory + 16 * wp * 1024 * 1024 < memoryUsed_)) {
      memoryUsed_ = totalMemory;
      bestPrefix_ = wp;
    }

    //  If we're vastly bigger than the 'best' size, stop.
    if (totalMemory > 16 * memoryUsed_)
      break;
  }
}



//  Find the values that match wPrefix == bestPrefix == _wPrefix.  It does
//  the same compute as findBestPrefix, but reports it and saves the porper
//  values.
//
void
merylOpCounting::findBestValues(uint64 nKmers, uint32  bestPrefix, uint64  memoryUsed) {
  uint32  merSize      = kmerTiny::merSize();
  uint32  segSizeBits  = merylCountArray::pagesPerSegment() * getPageSize() * 8;
  uint32  segSizeBytes = merylCountArray::pagesPerSegment() * getPageSize();

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "COMPLEX MODE\n");
  fprintf(stderr, "------------\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "prefix     # of   struct   kmers/    segs/      min     data    total\n");
  fprintf(stderr, "  bits   prefix   memory   prefix   prefix   memory   memory   memory\n");
  fprintf(stderr, "------  -------  -------  -------  -------  -------  -------  -------\n");

  //  Iterate over all reasonable prefix sizes emitting the table of values.
  //  When we hit the 'bestPrefix' save a copy of the values.

  for (uint32 wp=1; wp < 2 * merSize - 1; wp++) {
    uint64  nPrefix          = (uint64)1 << wp;                     //  Number of prefix == number of blocks of data
    uint64  kmersPerPrefix   = nKmers / nPrefix + 1;                //  Expected number of kmers we need to store per prefix
    uint64  kmersPerSeg      = segSizeBits / (2 * merSize - wp);    //  Kmers per segment
    uint64  segsPerPrefix    = kmersPerPrefix / kmersPerSeg + 1;    //

    if (wp + countNumberOfBits64(segsPerPrefix) + countNumberOfBits64(segSizeBytes) >= 64)
      break;   //  Otherwise, dataMemory overflows.

    uint64  structMemory     = ((sizeof(merylCountArray) * nPrefix) +                  //  Basic structs
                                (sizeof(uint64 *)        * nPrefix * segsPerPrefix));  //  Pointers to segments
    uint64  dataMemoryMin    = nPrefix *                 segSizeBytes;                 //  Minimum memory needed for this size.
    uint64  dataMemory       = nPrefix * segsPerPrefix * segSizeBytes;                 //  Expected memory for full batch.
    uint64  totalMemory      = structMemory + dataMemory;

    fprintf(stderr, "%6" F_U32P "  %4" F_U64P " %cP  %4" F_U64P " %cB  %4" F_U64P " %cM  %4" F_U64P " %cS  %4" F_U64P " %cB  %4" F_U64P " %cB  %4" F_U64P " %cB",
            wp,
            scaledNumber(nPrefix),        scaledUnit(nPrefix),
            scaledNumber(structMemory),   scaledUnit(structMemory),
            scaledNumber(kmersPerPrefix), scaledUnit(kmersPerPrefix),
            scaledNumber(segsPerPrefix),  scaledUnit(segsPerPrefix),
            scaledNumber(dataMemoryMin),  scaledUnit(dataMemoryMin),
            scaledNumber(dataMemory),     scaledUnit(dataMemory),
            scaledNumber(totalMemory),    scaledUnit(totalMemory));

    //  If this is the magic prefix, save all the values.
    if (wp == bestPrefix) {
      fprintf(stderr, "  Best Value!\n");

      _wPrefix       = wp;
      _wSuffix       = 2 * merSize - wp;

      _wSuffixMask   = 0;                               //  Build a mask by setting all bits
      _wSuffixMask   = ~_wSuffixMask;                   //  to one, then shifting in the
      _wSuffixMask >>= 8 * sizeof(kmdata) - _wSuffix;   //  correct number of zeros.

      _nPrefix       = nPrefix;
    }
    else {
      fprintf(stderr, "\n");
    }

    //  If we're vastly bigger than the 'best' size, stop.
    if (totalMemory > 16 * memoryUsed)
      break;
  }
}






//  Return a complete guess at the number of kmers in the input files.  No
//  rigorous went into the multipliers, just looked at a few sets of lambda reads.
uint64
merylOpCounting::guesstimateNumberOfkmersInInput_dnaSeqFile(dnaSeqFile *sequence) {
  uint64       numMers = 0;
  char const  *name    = sequence->filename();
  uint32       len     = strlen(name);

  if ((name[0] == '-') && (len == 1))
    return(0);

  uint64  size = merylutil::sizeOfFile(name);

  if      ((len > 3) && (name[len-3] == '.') && (name[len-2] == 'g') && (name[len-1] == 'z'))
    numMers += size * 3;

  else if ((len > 4) && (name[len-4] == '.') && (name[len-3] == 'b') && (name[len-2] == 'z') && (name[len-1] == '2'))
    numMers += size * 3.5;

  else if ((len > 3) && (name[len-3] == '.') && (name[len-2] == 'x') && (name[len-1] == 'z'))
    numMers += size * 4.0;

  else
    numMers += size;

  return(numMers);
}


uint64
merylOpCounting::guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID) {
  uint64  numMers = 0;

#ifdef CANU
  for (uint32 ii=bgnID; ii<endID; ii++)
    numMers += store->sqStore_getReadLength(ii);
#endif

  return(numMers);
}


uint64
merylOpCounting::guesstimateNumberOfkmersInInput(std::vector<merylInput *> &inputs) {
  uint64  guess = 0;

  for (uint32 ii=0; ii<inputs.size(); ii++) {
    if (inputs[ii]->isFromSequence())
      guess += guesstimateNumberOfkmersInInput_dnaSeqFile(inputs[ii]->_sequence);

    if (inputs[ii]->isFromStore())
      guess += guesstimateNumberOfkmersInInput_sqStore(inputs[ii]->_store, inputs[ii]->_sqBgn, inputs[ii]->_sqEnd);
  }

  return(guess);
}




//  Perform the counting operation, then close the output.
//
void
merylOpCounting::doCounting(std::vector<merylInput *> &inputs,
                            uint64                     memoryAllowed,
                            uint32                     threadsAllowed,
                            merylFileWriter           *output) {
  char    name[FILENAME_MAX + 1] = { 0 };
  bool    doSimple   = false;   //  Algorithm to use
  bool    doThreaded = true;    //  Always true

  //  Get this out of the way.

  omp_set_num_threads(threadsAllowed);

  //  Unless the user told us, make a guess on how many kmers the input will
  //  contain, then report what we're attempting to count.

  if (_expNumKmers == 0)
    _expNumKmers = guesstimateNumberOfkmersInInput(inputs);

  fprintf(stderr, "\n");
  fprintf(stderr, "Counting %lu (estimated)%s %s%s%s " F_U32 "-mers from " F_SIZE_T " input file%s:\n",
          scaledNumber(_expNumKmers), scaledName(_expNumKmers),
          (_countCanonical == true) ? "canonical" : "",
          (_countForward   == true) ? "forward" : "",
          (_countReverse   == true) ? "reverse" : "",
          kmerTiny::merSize(), inputs.size(), (inputs.size() == 1) ? "" : "s");

  for (uint32 ii=0; ii<inputs.size(); ii++)
    fprintf(stderr, "  %15s: %s\n", inputs[ii]->inputType(), inputs[ii]->inputName());

  //
  //  Set up to use the simple algorithm.
  //

  uint64  memoryUsedSimple = UINT64_MAX;

  findExpectedSimpleSize(_expNumKmers, _countSuffixString, _countSuffixLength, memoryUsedSimple);

  //
  //  Set up to use the complex algorithm.
  //
  uint64   memoryUsedComplex = UINT64_MAX;
  uint32   bestPrefix        = 0;
  uint32   nBatches          = 0;

#warning best prefix and memory size iterate over batches but this looks wrong
  for (nBatches=1; memoryUsedComplex > memoryAllowed; nBatches++)
    findBestPrefixSize(_expNumKmers / nBatches, memoryAllowed, bestPrefix, memoryUsedComplex);

  findBestValues(_expNumKmers / nBatches, bestPrefix, memoryUsedComplex);

  //
  //  Decide simple or complex.  useSimple_ is an output.
  //  If a count-suffix, we must use simple mode.
  //

  uint64  memoryUsed = 0;

  if ((memoryUsedSimple < memoryUsedComplex) &&
      (memoryUsedSimple < memoryAllowed)) {
    doSimple   = true;
    memoryUsed = memoryUsedSimple;
  }

  else {
    doSimple   = false;
    memoryUsed = memoryUsedComplex;
  }

  if (_countSuffixLength > 0) {
    doSimple    = true;
    memoryUsed  = memoryUsedSimple;
  }

  //
  //  Output the configuration.
  //

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "FINAL CONFIGURATION\n");
  fprintf(stderr, "-------------------\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Estimated to require %lu %cB memory out of %lu %cB allowed.\n",
          scaledNumber(memoryUsed),    scaledUnit(memoryUsed),
          scaledNumber(memoryAllowed), scaledUnit(memoryAllowed));
  fprintf(stderr, "Estimated to require %u batch%s.\n", nBatches, (nBatches == 1) ? "" : "es");
  fprintf(stderr, "\n");
  fprintf(stderr, "Configured %s mode for %.3f GB memory per batch, and up to %u batch%s.\n",      //  This is parsed
          (doSimple == true) ? "simple" : "complex",                                               //  by Canu.
          ((memoryUsed < memoryAllowed) ? memoryUsed : memoryAllowed) / 1024.0 / 1024.0 / 1024.0,  //  DO NOT CHANGE!
          nBatches, (nBatches == 1) ? "" : "es");
  fprintf(stderr, "\n");

  //
  //  Now, do some counting.
  //

  if (doSimple) {
    fprintf(stderr, "Start counting with SIMPLE method.\n");
    countSimple(inputs, memoryAllowed, threadsAllowed, output);
  }

  else if (doThreaded) {
    fprintf(stderr, "Start counting with THREADED method.\n");
    countThreads(inputs, memoryAllowed, threadsAllowed, output);
  }

  else {
    fprintf(stderr, "Start counting with SEQUENTIAL method.\n");
    countSequential(inputs, memoryAllowed, threadsAllowed, output);
  }
}
