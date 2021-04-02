
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


//
//  mcaSize       = sizeof(merylCountArray)  == 80
//  ptrSize       = sizeof(uint64 *)         == 8
//  segSize       = (as above)               ~= 64 KB
//
//  mersPerPrefix = nKmers / nPrefix+1
//  mersPerSeg    = 8 * segSize / 2*merSize-prefixSize
//
//  nSegPerPrefix = mersPerPrefix / mersPerSeg
//
//             basic structure             pointers to data                         data
//           |-----------------|   |------------------------------|   |------------------------------|
//           |                 |   |                              |   |                              |
//  memory = (mcaSize * nPrefix) + (ptrSize * nPrefix * nSegPrefix) + (segSize * nPrefix * nSegPrefix)
//
//  Solving for nKmers
//
//  memory = (mcaSize * nPrefix) + nSegPerPrefix * (ptrSize * nPrefix + segSize * nPrefix)
//
//  memory - (mcaSize * nPrefix) = nSegPerPrefix * (ptrSize * nPrefix + segSize * nPrefix)
//
//  nSegPerPrefix = (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + nPrefix * segSize)
//
//  (nKmers / nPrefix+1) / mersPerSeg = (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + segSize * nPrefix)
//   nKmers / nPrefix+1  = mersPerSeg * (memory - mcaSize * nPrefix) / (ptrSize * nPrefix + segSize * nPrefix)
#if 0
uint64
findMaxInputSizeForMemorySize(uint32 merSize, uint64 memSize) {
  uint64  mcaSize = sizeof(merylCountArray);
  uint64  ptrSize = sizeof(uint64 *);
  uint64  segSize = SEGMENT_SIZE * 1024;

  //  Free variable - prefixSize (wp)

  fprintf(stderr, "For memory size %lu bytes\n", memSize);
  fprintf(stderr, "                                                  |---------memory-breakdown--------\n");
  fprintf(stderr, "pBits    nKmers     nPrefix   nSegment   mers/seg  structure   pointers         data\n");
  fprintf(stderr, "----- ---------- ---------- ---------- ---------- ---------- ---------- ------------\n");

  for (uint32 wp=1; wp < 2*merSize; wp++) {
    uint64  nPrefix       = (uint64)1 << wp;

    if (mcaSize * nPrefix > memSize)
      break;

    //  Compute how many mer suffixes we can fit in a segment of memory.

    uint64  mersPerSeg = (segSize * 8) / (2 * merSize - wp);

    //  Compute the number of segments we can fit in memeory.  Each prefix must have a segment.  We
    //  don't actually know how many pointers we'll need -- it depends on the number of prefixes and
    //  number of segments -- so we assume there are 64 per prefix.
    //
    uint64  nSeg = (memSize - mcaSize * nPrefix - ptrSize * 64 * nPrefix) / segSize;

    if (nSeg < nPrefix)
      nSeg = nPrefix;

    //  Assuming some segment loading average, compute the number of kmers we can fit.  Each prefix
    //  has a 2/3 full segment, then any left over segments are 100% full.

    uint64  nKmers = nPrefix * mersPerSeg * 0.666;

    if (nPrefix < nSeg)
      nKmers += (nSeg - nPrefix) * mersPerSeg;

    //  For flavoring, show the memory breakdown.  nSegPrefix underflows, so needs to be included directly.

    uint64  ptrPerPrefix = 64; //32 * (nSeg / nPrefix + 1);  //  Somewhat incorrect.  The basic allocation is 32 pointers, then it doubles.

    uint64  basicMem     = mcaSize * nPrefix;                                    //  Size of the merCountArray structure itself.
    uint64  ptrMem       = ptrSize * nPrefix * ptrPerPrefix;                     //  Additional allocation for pointers to segments.
    uint64  dataMem      = segSize * nSeg;                                       //  Additional allocation of segments.

    if (basicMem + ptrMem + dataMem > 4 * memSize)
      break;

    fprintf(stderr, "%5u %10lu %10lu %10lu %10lu %10lu %10lu %12lu%s\n",
            wp, nKmers,
            nPrefix, nSeg,
            mersPerSeg,
            basicMem, ptrMem, dataMem,
            (basicMem + ptrMem + dataMem < memSize) ? "" : " INVALID");
  }

  exit(0);
}
#endif



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
//  estimate for the supplied nKmerEstimate.  If no estimate is below
//  memoryAllowed, 0 and UINT64_MAX, respectively, are returned.
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

  //  IMPORTANT!  Prefixes must be at least six bits - the number of bits
  //  we use for deciding on a file - and probably need to then leave one
  //  bit for a block id.  So only save a bestPrefix if it is at least 7.

  //  IMPORTANT!  Smaller prefixes mean bigger blocks in the output file,
  //  and any operation on those outputs needs to load all kmers in
  //  the block into memory, as full kmers.  Only save a bestPrefix
  //  if it is at least 10 - giving us at least 16 blocks per file.

  //  IMPORTANT!  Start searching at 1 and stop at merSize-1, otherwise
  //  we end up with a prefix or a suffix of size zero.

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

    //  Pick a larger prefix if it is dramatically smaller than what we have.
    //  More prefixes seem to run a bit slower, but also have smaller buckets
    //  for sorting at the end.

    if (structMemory + dataMemoryMin > memoryAllowed)
      break;

    if ((wp > 9) && (totalMemory + 16 * wp * 1024 * 1024 < memoryUsed_)) {
      memoryUsed_ = totalMemory;
      bestPrefix_ = wp;
    }

    if (totalMemory > 16 * memoryUsed_)
      break;
  }
}



void
findBestValues(uint64  nKmerEstimate,
               uint32  bestPrefix,
               uint64  memoryUsed,
               uint32 &wPrefix_,
               uint64 &nPrefix_,
               uint32 &wData_,
               kmdata &wDataMask_) {
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

  for (uint32 wp=1; wp < 2 * merSize - 1; wp++) {
    uint64  nPrefix          = (uint64)1 << wp;                          //  Number of prefix == number of blocks of data
    uint64  kmersPerPrefix   = nKmerEstimate / nPrefix + 1;              //  Expected number of kmers we need to store per prefix
    uint64  kmersPerSeg      = segSizeBits / (2 * merSize - wp);         //  Kmers per segment
    uint64  segsPerPrefix    = kmersPerPrefix / kmersPerSeg + 1;         //

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

    if (wp == bestPrefix) {
      fprintf(stderr, "  Best Value!\n");

      wPrefix_     = wp;
      nPrefix_     = nPrefix;
      wData_       = 2 * merSize - wp;

      wDataMask_   = 0;                              //  Build a mask by setting all bits
      wDataMask_   = ~wDataMask_;                    //  to one, then shifting in the
      wDataMask_ >>= 8 * sizeof(kmdata) - wData_;    //  correct number of zeros.

    } else {
      fprintf(stderr, "\n");
    }

    if (totalMemory > 16 * memoryUsed)
      break;
  }
}



void
merylOperation::configureCounting(uint64   memoryAllowed,      //  Input:  Maximum allowed memory in bytes
                                  bool    &useSimple_,         //  Output: algorithm to use
                                  uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
                                  uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
                                  uint32  &wData_,             //  Output: Number of bits in kmer data
                                  kmdata  &wDataMask_) {       //  Output: A mask to return just the data of the mer

  //
  //  Check kmer size, presence of output, and guess how many bases are in the inputs.
  //

  if (kmerTiny::merSize() == 0)
    fprintf(stderr, "ERROR: Kmer size not supplied with modifier k=<kmer-size>.\n"), exit(1);

  if ((_outputO == NULL) && (_onlyConfig == false))
    fprintf(stderr, "ERROR: No output specified for count operation.\n"), exit(1);

  if (_expNumKmers == 0)
    _expNumKmers = guesstimateNumberOfkmersInInput();

  //
  //  Report what we're trying.
  //

  fprintf(stderr, "\n");
  fprintf(stderr, "Counting %lu (estimated)%s %s%s%s " F_U32 "-mers from " F_SIZE_T " input file%s:\n",
          scaledNumber(_expNumKmers), scaledName(_expNumKmers),
          (_operation == opCount)        ? "canonical" : "",
          (_operation == opCountForward) ? "forward" : "",
          (_operation == opCountReverse) ? "reverse" : "",
          kmerTiny::merSize(), _inputs.size(), (_inputs.size() == 1) ? "" : "s");

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    fprintf(stderr, "  %15s: %s\n", _inputs[ii]->inputType(), _inputs[ii]->_name);

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

  for (nBatches=1; memoryUsedComplex > memoryAllowed; nBatches++)
    findBestPrefixSize(_expNumKmers / nBatches, memoryAllowed, bestPrefix, memoryUsedComplex);

  findBestValues(_expNumKmers / nBatches, bestPrefix, memoryUsedComplex, wPrefix_, nPrefix_, wData_, wDataMask_);

  //
  //  Decide simple or complex.  useSimple_ is an output.
  //  If a count-suffix, we must use simple mode.
  //

  uint64  memoryUsed = 0;

  if ((memoryUsedSimple < memoryUsedComplex) &&
      (memoryUsedSimple < memoryAllowed)) {
    useSimple_ = true;
    memoryUsed = memoryUsedSimple;
  }

  else {
    useSimple_ = false;
    memoryUsed = memoryUsedComplex;
  }

  if (_countSuffixLength > 0) {
    useSimple_  = true;
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
          (useSimple_ == true) ? "simple" : "complex",                                             //  by Canu.
          ((memoryUsed < memoryAllowed) ? memoryUsed : memoryAllowed) / 1024.0 / 1024.0 / 1024.0,  //  DO NOT CHANGE!
          nBatches, (nBatches == 1) ? "" : "es");
  fprintf(stderr, "\n");
}



//  Return a complete guess at the number of kmers in the input files.  No
//  rigorous went into the multipliers, just looked at a few sets of lambda reads.
uint64
merylOperation::guesstimateNumberOfkmersInInput_dnaSeqFile(dnaSeqFile *sequence) {
  uint64       numMers = 0;
  char const  *name    = sequence->filename();
  uint32       len     = strlen(name);

  if ((name[0] == '-') && (len == 1))
    return(0);

  uint64  size = AS_UTL_sizeOfFile(name);

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
merylOperation::guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID) {
  uint64  numMers = 0;

#ifdef CANU
  for (uint32 ii=bgnID; ii<endID; ii++)
    numMers += store->sqStore_getReadLength(ii);
#endif

  return(numMers);
}


uint64
merylOperation::guesstimateNumberOfkmersInInput(void) {
  uint64  guess = 0;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    if (_inputs[ii]->isFromSequence())
      guess += guesstimateNumberOfkmersInInput_dnaSeqFile(_inputs[ii]->_sequence);

    if (_inputs[ii]->isFromStore())
      guess += guesstimateNumberOfkmersInInput_sqStore(_inputs[ii]->_store, _inputs[ii]->_sqBgn, _inputs[ii]->_sqEnd);
  }

  return(guess);
}



void
merylOperation::count(uint32  wPrefix,
                      uint64  nPrefix,
                      uint32  wData,
                      kmdata  wDataMask) {

  //configureCounting(_maxMemory, useSimple, wPrefix, nPrefix, wData, wDataMask);

  //  If we're only configuring, stop now.

  if (_onlyConfig)
    return;

  //  Configure the writer for the prefix bits we're counting with.
  //
  //  We split the kmer into wPrefix and wData (bits) pieces.
  //  The prefix is used by the filewriter to decide which file to write to, by simply
  //  shifting it to the right to keep the correct number of bits in the file.
  //
  //           kmer -- [ wPrefix (18) = prefixSize               | wData (36) ]
  //           file -- [ numFileBits  | prefixSize - numFileBits ]

  _outputO->initialize(wPrefix);

  merylBlockWriter  *_writer = _outputP->getBlockWriter();

  //  Allocate buckets.  The buckets don't allocate space for mers until they're added,
  //  and allocate space for these mers in blocks of 8192 * 64 bits.
  //
  //  Need someway of balancing the number of prefixes we have and the size of each
  //  initial allocation.

  merylCountArray  *data = new merylCountArray [nPrefix];

  //  Load bases, count!

  uint64          bufferMax  = 1024 * 1024;
  uint64          bufferLen  = 0;
  char           *buffer     = new char [bufferMax];
  bool            endOfSeq   = false;

  memset(buffer, 0, sizeof(char) * bufferMax);

  //char            fstr[65];
  //char            rstr[65];

  uint64          memBase     = getProcessSize();   //  Overhead memory.
  uint64          memUsed     = 0;                  //  Sum of actual memory used.
  uint64          memReported = 0;                  //  Memory usage at last report.

  memUsed = memBase;

  for (uint32 pp=0; pp<nPrefix; pp++)
    memUsed += data[pp].initialize(pp, wData);

  uint64          kmersAdded  = 0;

  kmerIterator    kiter;

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    fprintf(stderr, "Loading kmers from '%s' into buckets.\n", _inputs[ii]->_name);

    while (_inputs[ii]->loadBases(buffer, bufferMax, bufferLen, endOfSeq)) {
      if (bufferLen == 0)
        continue;

      //fprintf(stderr, "read " F_U64 " bases from '%s'\n", bufferLen, _inputs[ii]->_name);

      kiter.addSequence(buffer, bufferLen);

      while (kiter.nextMer()) {
        bool    useF = (_operation == opCountForward);
        kmdata  pp   = 0;
        kmdata  mm   = 0;

        if (_operation == opCount)
          useF = (kiter.fmer() < kiter.rmer());

        if (useF == true) {
          pp = (kmdata)kiter.fmer() >> wData;
          mm = (kmdata)kiter.fmer()  & wDataMask;
          //fprintf(stderr, "useF F=%s R=%s ms=%u pp %lu mm %lu\n", kiter.fmer().toString(fstr), kiter.rmer().toString(rstr), kiter.fmer().merSize(), pp, mm);
        }

        else {
          pp = (kmdata)kiter.rmer() >> wData;
          mm = (kmdata)kiter.rmer()  & wDataMask;
          //fprintf(stderr, "useR F=%s R=%s ms=%u pp %lu mm %lu\n", kiter.fmer().toString(fstr), kiter.rmer().toString(rstr), kiter.rmer().merSize(), pp, mm);
        }

        assert(pp < nPrefix);

        memUsed += data[pp].add(mm);

        kmersAdded++;
      }

      if (endOfSeq)      //  If the end of the sequence, clear
        kiter.reset();   //  the running kmer.

      //  Report that we're actually doing something.

      if (memUsed - memReported > (uint64)128 * 1024 * 1024) {
        memReported = memUsed;

        fprintf(stderr, "Used %3.3f GB out of %3.3f GB to store %12lu kmers.\n",
                memUsed    / 1024.0 / 1024.0 / 1024.0,
                _maxMemory / 1024.0 / 1024.0 / 1024.0,
                kmersAdded);
      }

      //  If we're out of space, process the data and dump.

      if (memUsed > _maxMemory) {
        fprintf(stderr, "Memory full.  Writing results to '%s', using " F_S32 " threads.\n",
                _outputO->filename(), getMaxThreadsAllowed());
        fprintf(stderr, "\n");

#pragma omp parallel for schedule(dynamic, 1)
        for (uint32 ff=0; ff<_outputO->numberOfFiles(); ff++) {
          //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
          //        omp_get_thread_num(), ff, _outputO->firstPrefixInFile(ff), _outputO->lastPrefixInFile(ff));

          for (uint64 pp=_outputO->firstPrefixInFile(ff); pp <= _outputO->lastPrefixInFile(ff); pp++) {
            data[pp].countKmers();                //  Convert the list of kmers into a list of (kmer, count).
            data[pp].dumpCountedKmers(_writer);   //  Write that list to disk.
            data[pp].removeCountedKmers();        //  And remove the in-core data.
          }
        }

        _writer->finishBatch();

        kmersAdded = 0;

        memUsed = memBase;                        //  Reinitialize or memory used.
        for (uint32 pp=0; pp<nPrefix; pp++)
          memUsed += data[pp].usedSize();
      }

    }

    //  Would like some kind of report here on the kmers loaded from this file.

    delete _inputs[ii]->_sequence;
    _inputs[ii]->_sequence = NULL;
  }

  //  Finished loading kmers.  Free up some space.

  //delete [] kmers;
  delete [] buffer;

  //  Sort, dump and erase each block.
  //
  //  A minor complication is that within each output file, the blocks must be in order.

  fprintf(stderr, "\n");
  fprintf(stderr, "Writing results to '%s', using " F_S32 " threads.\n",
          _outputO->filename(), getMaxThreadsAllowed());

  //for (uint64 pp=0; pp<nPrefix; pp++)
  //  fprintf(stderr, "Prefix 0x%016lx writes to file %u\n", pp, _outputO->fileNumber(pp));

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<_outputO->numberOfFiles(); ff++) {
    //fprintf(stderr, "thread %2u writes file %2u with prefixes 0x%016lx to 0x%016lx\n",
    //        omp_get_thread_num(), ff, _outputO->firstPrefixInFile(ff), _outputO->lastPrefixInFile(ff));

    for (uint64 pp=_outputO->firstPrefixInFile(ff); pp <= _outputO->lastPrefixInFile(ff); pp++) {
      data[pp].countKmers();                //  Convert the list of kmers into a list of (kmer, count).
      data[pp].dumpCountedKmers(_writer);   //  Write that list to disk.
      data[pp].removeCountedKmers();        //  And remove the in-core data.
    }
  }

  //  Merge any iterations into a single file, or just rename
  //  the single file to the final name.

  _writer->finish();

  delete _writer;
  _writer = NULL;

  //  Cleanup.

  delete [] data;

  fprintf(stderr, "\n");
  fprintf(stderr, "Finished counting.\n");
}
