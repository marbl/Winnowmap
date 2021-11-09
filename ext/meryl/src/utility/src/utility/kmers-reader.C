
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


//  Clear all members and allocate buffers.
void
merylFileReader::initializeFromMasterI_v00(void) {

  _prefixSize    = 0;
  _suffixSize    = 0;

  _numFilesBits  = 0;
  _numBlocksBits = 0;

  _numFiles      = 0;
  _numBlocks     = 0;

  _stats         = NULL;

  _datFile       = NULL;

  _block         = new merylFileBlockReader();
  _blockIndex    = NULL;

  _kmer          = kmer();
  _value         = 0;

  _prefix        = 0;

  _activeMer     = 0;
  _activeFile    = 0;

  _threadFile    = UINT32_MAX;

  _nKmers        = 0;
  _nKmersMax     = 1024;
  _suffixes      = new kmdata [_nKmersMax];
  _values        = new kmvalu [_nKmersMax];
}



//  Initialize for the original.
void
merylFileReader::initializeFromMasterI_v01(stuffedBits  *masterIndex,
                                           bool          doInitialize) {

  if (doInitialize == true) {
    initializeFromMasterI_v00();

    _prefixSize    = masterIndex->getBinary(32);
    _suffixSize    = masterIndex->getBinary(32);

    _numFilesBits  = masterIndex->getBinary(32);
    _numBlocksBits = masterIndex->getBinary(32);

    _numFiles      = (uint64)1 << _numFilesBits;    //  The same for all formats, but
    _numBlocks     = (uint64)1 << _numBlocksBits;   //  awkward to do outside of here.
  }

  //  If we didn't initialize, set the file position to the start
  //  of the statistics.
  else {
    masterIndex->setPosition(64 + 64 + 32 + 32 + 32 + 32);
  }
}



//  Initialize for the format that includes multi sets.
void
merylFileReader::initializeFromMasterI_v02(stuffedBits  *masterIndex,
                                           bool          doInitialize) {

  if (doInitialize == true) {
    initializeFromMasterI_v00();

    _prefixSize    = masterIndex->getBinary(32);
    _suffixSize    = masterIndex->getBinary(32);

    _numFilesBits  = masterIndex->getBinary(32);
    _numBlocksBits = masterIndex->getBinary(32);

    uint32 flags   = masterIndex->getBinary(32);

    _isMultiSet    = flags & (uint32)0x0001;        //  This is new in v02.

    _numFiles      = (uint64)1 << _numFilesBits;    //  The same for all formats, but
    _numBlocks     = (uint64)1 << _numBlocksBits;   //  awkward to do outside of here.
  }

  //  If we didn't initialize, set the file position to the start
  //  of the statistics.
  else {
    masterIndex->setPosition(64 + 64 + 32 + 32 + 32 + 32 + 32);
  }
}



void
merylFileReader::initializeFromMasterI_v03(stuffedBits  *masterIndex,
                                           bool          doInitialize) {
  initializeFromMasterI_v02(masterIndex, doInitialize);
}



void
merylFileReader::initializeFromMasterIndex(bool  doInitialize,
                                           bool  loadStatistics,
                                           bool  beVerbose) {
  char   N[FILENAME_MAX+1];

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _inName);

  if (fileExists(N) == false)
    fprintf(stderr, "ERROR: '%s' doesn't appear to be a meryl input; file '%s' doesn't exist.\n",
            _inName, N), exit(1);

  //  Open the master index.

  stuffedBits  *masterIndex = new stuffedBits(N);

  //  Based on the magic number, initialzie.

  uint64  m1 = masterIndex->getBinary(64);
  uint64  m2 = masterIndex->getBinary(64);
  uint32  vv = 1;

  if        ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
             (m2 == 0x31302e765f5f7865llu)) {   //  ex__v.01
    initializeFromMasterI_v01(masterIndex, doInitialize);
    vv = 1;

  } else if ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
             (m2 == 0x32302e765f5f7865llu)) {   //  ex__v.02
    initializeFromMasterI_v02(masterIndex, doInitialize);
    vv = 2;

  } else if ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
             (m2 == 0x33302e765f5f7865llu)) {   //  ex__v.03
    initializeFromMasterI_v03(masterIndex, doInitialize);
    vv = 3;

  } else {
    fprintf(stderr, "ERROR: '%s' doesn't look like a meryl input; file '%s' fails magic number check.\n",
            _inName, N), exit(1);
  }

  //  Check that the mersize is set and valid.

  uint32  merSize = (_prefixSize + _suffixSize) / 2;

  if (kmer::merSize() == 0)         //  If the global kmer size isn't set yet,
    kmer::setSize(merSize);         //  set it.

  if (kmer::merSize() != merSize)   //  And if set, make sure we're compatible.
    fprintf(stderr, "mer size mismatch, can't process this set of files.\n"), exit(1);

  //  If loading statistics is enabled, load the stats assuming the file is in
  //  the proper position.

  if (loadStatistics == true) {
    _stats = new merylHistogram;
    _stats->load(masterIndex, vv);
  }

  //  And report some logging.

  if (beVerbose) {
    char    m[17] = { 0 };

    for (uint32 i=0, s=0; i<8; i++, s+=8) {
      m[i + 0] = (m1 >> s) & 0xff;
      m[i + 8] = (m2 >> s) & 0xff;
    }

    fprintf(stderr, "Opened '%s'.\n", _inName);
    fprintf(stderr, "  magic          0x%016lx%016lx '%s'\n", m1, m2, m);
    fprintf(stderr, "  prefixSize     %u\n", _prefixSize);
    fprintf(stderr, "  suffixSize     %u\n", _suffixSize);
    fprintf(stderr, "  numFilesBits   %u (%u files)\n", _numFilesBits, _numFiles);
    fprintf(stderr, "  numBlocksBits  %u (%u blocks)\n", _numBlocksBits, _numBlocks);
  }

  delete masterIndex;
}



merylFileReader::merylFileReader(const char *inputName,
                                         bool        beVerbose) {
  strncpy(_inName, inputName, FILENAME_MAX);
  initializeFromMasterIndex(true, false, beVerbose);
}



merylFileReader::merylFileReader(const char *inputName,
                                         uint32      threadFile,
                                         bool        beVerbose) {
  strncpy(_inName, inputName, FILENAME_MAX);
  initializeFromMasterIndex(true, false, beVerbose);
  enableThreads(threadFile);
}



merylFileReader::~merylFileReader() {

  delete [] _blockIndex;

  delete [] _suffixes;
  delete [] _values;

  delete    _stats;

  AS_UTL_closeFile(_datFile);

  delete    _block;
}



void
merylFileReader::loadStatistics(void) {
  if (_stats == NULL)
    initializeFromMasterIndex(false, true, false);
}



void
merylFileReader::dropStatistics(void) {
  delete _stats;
  _stats = NULL;
}



void
merylFileReader::enableThreads(uint32 threadFile) {
  _activeFile = threadFile;
  _threadFile = threadFile;
}



void
merylFileReader::loadBlockIndex(void) {

  if (_blockIndex != NULL)
    return;

  _blockIndex = new merylFileIndex [_numFiles * _numBlocks];

  for (uint32 ii=0; ii<_numFiles; ii++) {
    char  *idxname = constructBlockName(_inName, ii, _numFiles, 0, true);
    FILE  *idxfile = AS_UTL_openInputFile(idxname);

    loadFromFile(_blockIndex + _numBlocks * ii, "merylFileReader::blockIndex", _numBlocks, idxfile);

    AS_UTL_closeFile(idxfile, idxname);

    delete [] idxname;
  }
}



//  Like loadBlock, but just reports all blocks in the file, ignoring
//  the kmer data.
//
void
dumpMerylDataFile(char *name) {
  FILE            *F = NULL;
  merylFileIndex   I;
  stuffedBits     *D = NULL;

  //  Dump the merylIndex for this block.

  if (fileExists(name, '.', "merylIndex") == false)
    fprintf(stderr, "ERROR: '%s.merylIndex' doesn't exist.  Can't dump it.\n",
            name), exit(1);

  F = AS_UTL_openInputFile(name, '.', "merylIndex");

  fprintf(stdout, "\n");
  fprintf(stdout, "    prefix    blkPos    nKmers\n");
  fprintf(stdout, "---------- --------- ---------\n");

  while (loadFromFile(I, "merylFileIndex", F, false) != 0) {
    fprintf(stdout, "0x%08x %9lu %9lu\n", I.blockPrefix(), I.blockPosition(), I.numKmers());
  }

  AS_UTL_closeFile(F);

  //  Read each block, sequentially, and report the header.

  if (fileExists(name, '.', "merylData") == false)
    fprintf(stderr, "ERROR: '%s.merylData' doesn't exist.  Can't dump it.\n",
            name), exit(1);

  F = AS_UTL_openInputFile(name, '.', "merylData");
  D = new stuffedBits;

  fprintf(stdout, "\n");
  fprintf(stdout, "            prefix   nKmers kCode uBits bBits                 k1 cCode                 c1                 c2\n");
  fprintf(stdout, "------------------ -------- ----- ----- ----- ------------------ ----- ------------------ ------------------\n");

  while (D->loadFromFile(F)) {
    uint64 position   = D->getPosition();

    uint64 m1         = D->getBinary(64);
    uint64 m2         = D->getBinary(64);

    uint64 prefix     = D->getBinary(64);
    uint64 nKmers     = D->getBinary(64);

    uint8  kCode      = D->getBinary(8);
    uint32 unaryBits  = D->getBinary(32);
    uint32 binaryBits = D->getBinary(32);
    uint64 k1         = D->getBinary(64);

    uint8  cCode      = D->getBinary(8);
    uint64 c1         = D->getBinary(64);
    uint64 c2         = D->getBinary(64);

    if ((m1 != 0x7461446c7972656dllu) ||
        (m2 != 0x0a3030656c694661llu)) {
      fprintf(stderr, "merylFileReader::nextMer()-- Magic number mismatch at position " F_U64 ".\n", position);
      fprintf(stderr, "merylFileReader::nextMer()-- Expected 0x7461446c7972656d got 0x%016" F_X64P "\n", m1);
      fprintf(stderr, "merylFileReader::nextMer()-- Expected 0x0a3030656c694661 got 0x%016" F_X64P "\n", m2);
      exit(1);
    }

    fprintf(stdout, "0x%016lx %8lu %5u %5u %5u 0x%016lx %5u 0x%016lx 0x%016lx\n",
            prefix, nKmers, kCode, unaryBits, binaryBits, k1, cCode, c1, c2);
  }

  delete D;

  AS_UTL_closeFile(F);

  //  Read each block again, dump the kmers in the block.

  F = AS_UTL_openInputFile(name, '.', "merylData");
  D = new stuffedBits;

  while (D->loadFromFile(F)) {
    uint64 position   = D->getPosition();

    uint64 m1         = D->getBinary(64);
    uint64 m2         = D->getBinary(64);

    uint64 prefix     = D->getBinary(64);
    uint64 nKmers     = D->getBinary(64);

    uint8  kCode      = D->getBinary(8);
    uint32 unaryBits  = D->getBinary(32);
    uint32 binaryBits = D->getBinary(32);
    uint64 k1         = D->getBinary(64);

    uint8  cCode      = D->getBinary(8);
    uint64 c1         = D->getBinary(64);
    uint64 c2         = D->getBinary(64);

    fprintf(stdout, "\n");
    fprintf(stdout, " kmerIdx prefixDelta      prefix |--- suffix-size and both suffixes ---|    value\n");
    fprintf(stdout, "-------- ----------- ----------- -- ---------------- -- ---------------- --------\n");

    uint64   *pd = new uint64 [nKmers];
    uint64   *s1 = new uint64 [nKmers];
    uint64   *s2 = new uint64 [nKmers];
    uint64   *va = new uint64 [nKmers];

    uint32    ls = (binaryBits <= 64) ? (0)          : (binaryBits - 64);
    uint32    rs = (binaryBits <= 64) ? (binaryBits) : (64);

    uint64    tp = 0;

    //  Get all the kmers.
    for (uint32 kk=0; kk<nKmers; kk++) {
      if (kCode == 1) {
        pd[kk] = D->getUnary();
        s1[kk] = D->getBinary(ls);
        s2[kk] = D->getBinary(rs);
      }

      else {
        fprintf(stderr, "ERROR: unknown kCode %u\n", kCode), exit(1);
      }
    }

    //  Get all the values.
    for (uint32 kk=0; kk<nKmers; kk++) {
      if      (cCode == 1) {
        va[kk] = D->getBinary(32);
      }

      else if (cCode == 2) {
        va[kk] = D->getBinary(64);
      }

      else {
        fprintf(stderr, "ERROR: unknown cCode %u\n", cCode), exit(1);
      }
    }

    //  Dump.
    for (uint32 kk=0; kk<nKmers; kk++) {
      tp += pd[kk];

      fprintf(stdout, "%8u %11lu %011lx %2u %016lx %2u %016lx %8lx\n",
              kk, pd[kk], tp, ls, s1[kk], rs, s2[kk], va[kk]);
    }
  }

  delete D;

  AS_UTL_closeFile(F);
}



bool
merylFileReader::nextMer(void) {

  _activeMer++;

  //  If we've still got data, just update and get outta here.
  //  Otherwise, we need to load another block.

  if (_activeMer < _nKmers) {
    _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
    _value = _values[_activeMer];
    return(true);
  }

  //  If no file, open whatever is 'active'.  In thread mode, the first file
  //  we open is the 'threadFile'; in normal mode, the first file we open is
  //  the first file in the database.

 loadAgain:
  if (_datFile == NULL)
    _datFile = openInputBlock(_inName, _activeFile, _numFiles);

  //  Load blocks.

  bool loaded = _block->loadBlock(_datFile, _activeFile);

  //  If nothing loaded. open a new file and try again.

  if (loaded == false) {
    AS_UTL_closeFile(_datFile);

    if (_activeFile == _threadFile)   //  Thread mode, if no block was loaded,
      return(false);                  //  we're done.

    _activeFile++;

    if (_numFiles <= _activeFile)
      return(false);

    goto loadAgain;
  }

  //  Got a block!  Stash what we loaded.

  _prefix = _block->prefix();
  _nKmers = _block->nKmers();

#ifdef SHOW_LOAD
  fprintf(stdout, "LOADED prefix %016lx nKmers %lu\n", _prefix, _nKmers);
#endif

  //  Make sure we have space for the decoded data

  resizeArrayPair(_suffixes, _values, 0, _nKmersMax, _nKmers, _raAct::doNothing);

  //  Decode the block into _OUR_ space.
  //
  //  decodeBlock() marks the block as having no data, so the next time we loadBlock() it will
  //  read more data from disk.  For blocks that don't get decoded, they retain whatever was
  //  loaded, and do not load another block in loadBlock().

  _block->decodeBlock(_suffixes, _values);

  //  But if no kmers in this block, load another block.  Sadly, the block must always
  //  be decoded, otherwise, the load will not load a new block.

  if (_nKmers == 0)
    goto loadAgain;

  //  Reset iteration, and load the first kmer.

  _activeMer = 0;

  _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
  _value = _values[_activeMer];

  return(true);
}
