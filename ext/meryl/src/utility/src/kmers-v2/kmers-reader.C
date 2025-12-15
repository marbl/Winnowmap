
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


stuffedBits *
merylFileReader::openMasterIndex(void) {
  char   N[FILENAME_MAX+1];

  snprintf(N, FILENAME_MAX, "%s/merylIndex", _inName);

  if (fileExists(N) == false)
    return nullptr;
  else
    return new stuffedBits(N);
}

//  Clear all members and allocate buffers.
void
merylFileReader::initializeFromMasterI_v00(void) {

  _block         = new merylFileBlockReader();
  _blockIndex    = nullptr;

  _nKmers        = 0;
  _nKmersMax     = 1024;
  _suffixes      = new kmdata [_nKmersMax];
  _values        = new kmvalu [_nKmersMax];
  _labels        = new kmlabl [_nKmersMax];
}



//  Initialize for the original.
void
merylFileReader::initializeFromMasterI_v01(stuffedBits  *masterIndex) {

  initializeFromMasterI_v00();

  _prefixSize    = masterIndex->getBinary(32);
  _suffixSize    = masterIndex->getBinary(32);

  _numFilesBits  = masterIndex->getBinary(32);
  _numBlocksBits = masterIndex->getBinary(32);

  _numFiles      = (uint64)1 << _numFilesBits;    //  The same for all formats, but
  _numBlocks     = (uint64)1 << _numBlocksBits;   //  awkward to do outside of here.

  _statsVersion  = 1;
  _statsOffset   = masterIndex->getPosition();   assert(_statsOffset == 64 + 64 + 32 + 32 + 32 + 32);
}



//  Initialize for the format that includes multi sets.
void
merylFileReader::initializeFromMasterI_v02(stuffedBits  *masterIndex) {

  initializeFromMasterI_v00();

  _prefixSize    = masterIndex->getBinary(32);
  _suffixSize    = masterIndex->getBinary(32);

  _numFilesBits  = masterIndex->getBinary(32);
  _numBlocksBits = masterIndex->getBinary(32);

  uint32 flags   = masterIndex->getBinary(32);

  _isMultiSet    = flags & (uint32)0x0001;        //  This is new in v02.

  _numFiles      = (uint64)1 << _numFilesBits;    //  The same for all formats, but
  _numBlocks     = (uint64)1 << _numBlocksBits;   //  awkward to do outside of here.

  _statsVersion  = 2;
  _statsOffset   = masterIndex->getPosition();   assert(_statsOffset == 64 + 64 + 32 + 32 + 32 + 32 + 32);
}



void
merylFileReader::initializeFromMasterI_v03(stuffedBits  *masterIndex) {
  initializeFromMasterI_v02(masterIndex);

  _statsVersion  = 3;
}



void
merylFileReader::initializeFromMasterI_v04(stuffedBits  *masterIndex) {
  initializeFromMasterI_v02(masterIndex);

  _statsVersion  = 4;
}



bool
merylFileReader::initializeFromMasterIndex(std::vector<char const *> *errors) {

  stuffedBits  *masterIndex = openMasterIndex();

  if (masterIndex == nullptr)
    return fatalError(errors, "Input '%s' isn't a meryl database; master index file not found.", _inName);

  uint64  m1 = masterIndex->getBinary(64);    //  Based on the magic number in the index,
  uint64  m2 = masterIndex->getBinary(64);    //  initialize!

  if      ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
           (m2 == 0x31302e765f5f7865llu))     //  ex__v.01
    initializeFromMasterI_v01(masterIndex);

  else if ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
           (m2 == 0x32302e765f5f7865llu))     //  ex__v.02
    initializeFromMasterI_v02(masterIndex);

  else if ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
           (m2 == 0x33302e765f5f7865llu))     //  ex__v.03
    initializeFromMasterI_v03(masterIndex);

  else if ((m1 == 0x646e496c7972656dllu) &&   //  merylInd
           (m2 == 0x34302e765f5f7865llu))     //  ex__v.04
    initializeFromMasterI_v04(masterIndex);

  delete masterIndex;

  if (_statsVersion == 0)                     //  Failed to initialize.
    return fatalError(errors, "Input '%s' isn't a meryl database; magic number check failed.", _inName);

  uint32  merSize = (_prefixSize + _suffixSize) / 2;

  if (kmer::merSize() == 0)                   //  If the global kmer size isn't set yet,
    kmer::setSize(merSize);                   //  set it.

  if (kmer::merSize() != merSize)             //  And if set, make sure we're compatible.
    return fatalError(errors, "Input '%s' contains %u-mers, but expecting %u-mers.",
                      _inName, merSize, kmer::merSize());

  return true;
}



merylFileReader::merylFileReader(const char *inputName, std::vector<char const *> *errors) {
  strncpy(_inName, inputName, FILENAME_MAX);
  initializeFromMasterIndex(errors);
}



merylFileReader::merylFileReader(const char *inputName,
                                 uint32      threadFile, std::vector<char const *> *errors) {
  strncpy(_inName, inputName, FILENAME_MAX);
  initializeFromMasterIndex(errors);
  enableThreads(threadFile);
}



merylFileReader::~merylFileReader() {

  delete [] _blockIndex;

  delete [] _suffixes;
  delete [] _values;
  delete [] _labels;

  delete    _stats;

  merylutil::closeFile(_datFile);

  delete    _block;
}



void
merylFileReader::loadStatistics(std::vector<char const *> *errors) {

  if (_stats)    //  Already have stats, don't load again.
    return;

  stuffedBits  *masterIndex = openMasterIndex();

  if (masterIndex) {
    masterIndex->setPosition(_statsOffset);

    _stats = new merylHistogram;
    _stats->load(masterIndex, _statsVersion);
  }
  else {
    fatalError(errors, "Failed to load statistics for input '%s'; master index file not found.", _inName);
  }

  delete masterIndex;
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
    FILE  *idxfile = merylutil::openInputFile(idxname);

    loadFromFile(_blockIndex + _numBlocks * ii, "merylFileReader::blockIndex", _numBlocks, idxfile);

    merylutil::closeFile(idxfile, idxname);

    delete [] idxname;
  }
}



bool
merylFileReader::nextMer(void) {

  _activeMer++;

  //  If we've still got data, just update and get outta here.
  //  Otherwise, we need to load another block.

  if (_activeMer < _nKmers) {
    _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
    _kmer._val = _values[_activeMer];
    return(true);
  }

  //  If no file, open whatever is 'active'.  In thread mode, the first file
  //  we open is the 'threadFile'; in normal mode, the first file we open is
  //  the first file in the database.

 loadAgain:
  if (_datFile == NULL)
    _datFile = openInputBlock(_inName, _activeFile, _numFiles);

  //  Load blocks.

  bool loaded = _block->loadKmerFileBlock(_datFile, _activeFile);

  //  If nothing loaded. open a new file and try again.

  if (loaded == false) {
    merylutil::closeFile(_datFile);

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

  //  Make sure we have space for the decoded data

  resizeArray(_suffixes, _values, _labels, 0, _nKmersMax, _nKmers, _raAct::doNothing);

  //  Decode the block into _OUR_ space.
  //
  //  decodeKmerFileBlock() marks the block as having no data, so the next
  //  time we loadBlock() it will read more data from disk.  For blocks that
  //  don't get decoded, they retain whatever was loaded, and do not load
  //  another block in loadBlock().

  _block->decodeKmerFileBlock(_suffixes, _values, _labels);

  //  But if no kmers in this block, load another block.  Sadly, the block must always
  //  be decoded, otherwise, the load will not load a new block.

  if (_nKmers == 0)
    goto loadAgain;

  //  Reset iteration, and load the first kmer.

  _activeMer = 0;

  _kmer.setPrefixSuffix(_prefix, _suffixes[_activeMer], _suffixSize);
  _kmer._val = _values[_activeMer];
  _kmer._lab = _labels[_activeMer];

  return(true);
}

}  //  namespace merylutil::kmers::v2
