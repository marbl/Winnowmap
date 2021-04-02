
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


void
merylFileWriter::initialize(uint32 prefixSize, bool isMultiSet) {

  //  Fail if we're already initialized and asked to change the prefix size.
  //  But just ignore the re-init request if the prefix size is the same.

  if ((_initialized == true) &&
      (prefixSize != _prefixSize))
    fprintf(stderr, "merylFileWriter::initialize()-- asked to initialize with different prefixSize (new %u existing %u).\n", prefixSize, _prefixSize), exit(1);

  if (_initialized == true)
    return;

  //  If the global mersize isn't set, we're hosed.

  if (kmer::merSize() == 0)
    fprintf(stderr, "merylFileWriter::initialize()-- asked to initialize, but kmer::merSize() is zero!\n"), exit(1);

  //  The count operations call initialize() exactly once, but nextMer() calls
  //  it once per file and so we need some kind of concurrency control here.

#pragma omp critical (merylFileWriterInit)
  if (_initialized == false) {

    //  If the prefixSize is zero, set it to (arbitrary) 1/4 the kmer size.
    //  This happens in the streaming writer (which is used when meryl does
    //  any non-count operation).  The prefixSize here just controls how
    //  often we dump blocks to the file.

    if (_prefixSize == 0)
      _prefixSize = prefixSize;

#warning how to set prefix size for streaming operations?
    if (_prefixSize == 0)
      _prefixSize = 12;  //max((uint32)8, 2 * kmer::merSize() / 3);

    _suffixSize         = 2 * kmer::merSize() - _prefixSize;
    _suffixMask         = buildLowBitMask<kmdata>(_suffixSize);

    //  Decide how many files to write.  We can make up to 2^32 files, but will
    //  run out of file handles _well_ before that.  For now, limit to 2^6 = 64 files.

    _numFilesBits       = 6;
    _numBlocksBits      = _prefixSize - _numFilesBits;

    _numFiles           = (uint64)1 << _numFilesBits;
    _numBlocks          = (uint64)1 << _numBlocksBits;

    _isMultiSet         = isMultiSet;

    //  Now we're initialized!

    //fprintf(stderr, "merylFileWriter()-- Creating '%s' for %u-mers, with prefixSize %u suffixSize %u numFiles %lu\n",
    //        _outName, (_prefixSize + _suffixSize) / 2, _prefixSize, _suffixSize, _numFiles);

    _initialized = true;
  }
}



merylFileWriter::merylFileWriter(const char *outputName,
                                 uint32      prefixSize) {

  //  Note that we're not really initialized yet.  We could call initialize() in some cases,
  //  but the interesting one can't initialized() until the first meryl input file is opened,
  //  so we don't initialize any of them.

  _initialized   = false;

  //  Save the output directory name, and try to make it.  If we can't we'll fail quickly.

  strncpy(_outName, outputName, FILENAME_MAX);

  AS_UTL_mkdir(_outName);

  //  Parameters on how the suffixes/values are encoded are set once we know
  //  the kmer size.  See initialize().

  _prefixSize    = prefixSize;

  _suffixSize    = 0;
  _suffixMask    = 0;

  _numFilesBits  = 0;
  _numBlocksBits = 0;
  _numFiles      = 0;
  _numBlocks     = 0;

  _isMultiSet    = false;
}



merylFileWriter::~merylFileWriter() {
  uint32   flags = (uint32)0x0000;

  //  Set flags.

  if (_isMultiSet)
    flags |= (uint32)0x0001;

  //  Create a master index with the parameters.

  stuffedBits  *masterIndex = new stuffedBits(32 * 1024);

  masterIndex->setBinary(64, 0x646e496c7972656dllu);  //  HEX: ........  ONDISK: merylInd
  masterIndex->setBinary(64, 0x33302e765f5f7865llu);  //       30.v__xe          ex__v.03
  masterIndex->setBinary(32, _prefixSize);
  masterIndex->setBinary(32, _suffixSize);
  masterIndex->setBinary(32, _numFilesBits);
  masterIndex->setBinary(32, _numBlocksBits);
  masterIndex->setBinary(32, flags);

  _stats.dump(masterIndex);

  //  Store the master index (and stats) to disk.

  char     N[FILENAME_MAX+1];
  FILE    *F;

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _outName);

  F = AS_UTL_openOutputFile(N);
  masterIndex->dumpToFile(F);
  AS_UTL_closeFile(F);

  delete masterIndex;
}



uint32
merylFileWriter::fileNumber(uint64  prefix) {

  assert(_initialized);

  //  Based on the prefix, decide what output file to write to.
  //  The prefix has _prefixSize = _numFilesBits + _numBlocksBits.
  //  We want to save the highest _numFilesBits.

  uint64  oi  = prefix >> _numBlocksBits;

  if (oi >= _numFiles) {
    fprintf(stderr, "merylFileWriter()-- Formed invalid file number %lu >= number of files %lu:\n", oi, _numFiles);
    fprintf(stderr, "merylFileWriter()--   prefix          0x%016lx\n", prefix);
    fprintf(stderr, "merylFileWriter()--   prefixSize      %u\n", _prefixSize);
    fprintf(stderr, "merylFileWriter()--   suffixSize      %u\n", _suffixSize);
    fprintf(stderr, "merylFileWriter()--   numFilesBits    %u\n", _numFilesBits);
    fprintf(stderr, "merylFileWriter()--   numBlocksBits   %u\n", _numBlocksBits);
  }
  assert(oi < _numFiles);

  return((uint32)oi);
}



void
merylFileWriter::writeBlockToFile(FILE            *datFile,
                                  merylFileIndex  *datFileIndex,
                                  kmpref           blockPrefix,
                                  uint64           nKmers,
                                  kmdata          *suffixes,
                                  kmvalu          *values) {

  //  Figure out the optimal size of the Elias-Fano prefix.  It's just log2(N)-1.

  uint32  unaryBits = 0;
  uint64  unarySum  = 1;
  while (unarySum < nKmers) {
    unaryBits  += 1;
    unarySum  <<= 1;
  }

  uint32  binaryBits = _suffixSize - unaryBits;      //  Only _suffixSize is used from the class.

  //  Decide how to encode the data.
  //
  //    kmer coding type 1 == Elias Fano
  //
  //    valu coding type 1 == 32-bit binary data
  //    valu coding type 2 == 64-bit binary data

  uint32  kct = 1;
  uint32  vct = sizeof(kmvalu) / 4;

  //  Dump data.
  //
  //  Unary coding requires that the whole value fit in one block.  The
  //  largest unary number we'll see is 2^unaryBits -- and that's already
  //  computed in unarySum.  We'll set the stuffedBits block size
  //  to 

  uint64         blockSize;

  blockSize  = 10 * 64;                    //  For the header.
  blockSize += 2 * unarySum;               //  For the unary encoded prefix bits
  blockSize += nKmers * binaryBits / 16;   //  For the binary encoded suffix bits
  blockSize += nKmers * 32 / 16;           //  For the value bits

  blockSize = (blockSize & 0xfffffffffffffc00llu) + 1024;   //  Make it a multiple of 1024.

  stuffedBits   *dumpData = new stuffedBits(blockSize);

  dumpData->setBinary(64, 0x7461446c7972656dllu);    //  Magic number, part 1.
  dumpData->setBinary(64, 0x0a3030656c694661llu);    //  Magic number, part 2.

  dumpData->setBinary(64, blockPrefix);
  dumpData->setBinary(64, nKmers);

  dumpData->setBinary(8,  kct);                      //  Kmer coding type
  dumpData->setBinary(32, unaryBits);                //  Kmer coding parameters
  dumpData->setBinary(32, binaryBits);
  dumpData->setBinary(64, 0);

  dumpData->setBinary(8,  vct);                      //  Value coding type
  dumpData->setBinary(64, 0);                        //  Value coding parameters
  dumpData->setBinary(64, 0);

  //  Split the kmer suffix into two pieces, one unary encoded offsets and one binary encoded.

  uint64  lastPrefix = 0;
  uint64  thisPrefix = 0;

  assert(kct == 1);  //  Eventually could add more...

  for (uint32 kk=0; kk<nKmers; kk++) {
    thisPrefix = suffixes[kk] >> binaryBits;

    uint64  l = suffixes[kk] >> 64;
    uint64  r = suffixes[kk];

    uint32 ls = (binaryBits <= 64) ? (0)          : (binaryBits - 64);
    uint32 rs = (binaryBits <= 64) ? (binaryBits) : (64);

    dumpData->setUnary(thisPrefix - lastPrefix);
    dumpData->setBinary(ls, l);
    dumpData->setBinary(rs, r);

    lastPrefix = thisPrefix;
  }

  //  Save the values, too.  Eventually these will be cleverly encoded.  Really.

  uint64  lastValue = 0;
  uint64  thisValue = 0;

  assert((vct == 1) || (vct == 2));

  for (uint32 kk=0; kk<nKmers; kk++) {
    dumpData->setBinary(32 * vct, values[kk]);
  }

  //  Save the index entry.

  uint64  block = blockPrefix & buildLowBitMask<uint64>(_numBlocksBits);

  datFileIndex[block].set(blockPrefix, datFile, nKmers);

  //  Dump data to disk, cleanup, and done!

  dumpData->dumpToFile(datFile);

  delete dumpData;
}
