
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

#include "files.H"
#include "bits.H"
#include "stuffedBits-v1.H"


inline
uint64
bitsToWords(uint64 bits) {
  return(bits / 64 + ((bits % 64) ? 1 : 0));
}




namespace merylutil::inline bits::inline v1 {

stuffedBits::stuffedBits(uint64 nBits) {
  _maxBits = roundMaxSizeUp(nBits);
  allocateBlock();
}

stuffedBits::stuffedBits(char const *inputName) {
  FILE *inFile = openInputFile(inputName);
  load(inFile, nullptr);
  closeFile(inFile);
}

stuffedBits::stuffedBits(FILE *inFile) {
  load(inFile, nullptr);
}

stuffedBits::stuffedBits(readBuffer *B) {
  load(nullptr, B);
}

stuffedBits::~stuffedBits() {
  for (uint32 ii=0; ii<_blocksMax; ii++)
    delete [] _blocks[ii]._dat;
  delete [] _blocks;
}




void
stuffedBits::dump(FILE *F, writeBuffer *B) {
  uint32   outLen = 0;

  //  Count the number of blocks we're going to write.

  for (outLen=0; ((outLen < _blocksMax) &&
                  (_blocks[outLen]._len > 0)); outLen++)
    ;

  if (F) {
    writeToFile(_maxBits,   "dataBlockLenMax", F);
    writeToFile(outLen,     "blocksLen",       F);
    writeToFile(_blocksMax, "blocksMax",       F);
  }

  if (B) {
    B->write(&_maxBits,   sizeof(uint64));
    B->write(&outLen,     sizeof(uint32));
    B->write(&_blocksMax, sizeof(uint32));
  }

  //  Copy the begin positions and lengths to arrays, then dump them as such.
  //  This is because the original format wrote this way, and I wasn't smart
  //  enough to put a version number in the stream.

  uint64  *ua = new uint64 [outLen];

  for (uint32 ii=0; ii<outLen; ii++)
    ua[ii] = _blocks[ii]._bgn;
  if (F)   writeToFile(ua, "dataBlockBgn", outLen, F);
  if (B)   B->write(ua, sizeof(uint64) * outLen);

  for (uint32 ii=0; ii<outLen; ii++)
    ua[ii] = _blocks[ii]._len;
  if (F)   writeToFile(ua, "dataBlockLen", outLen, F);
  if (B)   B->write(ua, sizeof(uint64) * outLen);

  delete [] ua;

  //  Write any block with data.

  for (outLen=0; ((outLen < _blocksMax) &&
                  (_blocks[outLen]._len > 0)); outLen++) {
    if (F)   writeToFile(_blocks[outLen]._dat, "dataBlocks", bitsToWords(_blocks[outLen]._len), F);
    if (B)   B->write(_blocks[outLen]._dat, sizeof(uint64) * bitsToWords(_blocks[outLen]._len));
  }
}



bool
stuffedBits::load(FILE *F, readBuffer *B) {
  uint32   nLoad    = 0;
  uint32   inLen    = 0;   //  Number of blocks we need to load.
  uint32   inMax    = 0;   //  Maximum number of blocks to allocate, not used here.

  eraseBlocks();

  //  Try to load the parameters of the block.  If any fail to read, we've
  //  hit the end-of-buffer and there is no valid stuffedBits to load.

  if ((F == nullptr) &&
      (B == nullptr))
    return(false);

  if ((F != nullptr) &&
      ((1 != merylutil::loadFromFile(_maxBits, "dataBlockLenMax", F, false)) ||
       (1 != merylutil::loadFromFile(inLen,    "blocksLen",       F, false)) ||
       (1 != merylutil::loadFromFile(inMax,    "blocksMax",       F, false))))
    return(false);

  if ((B != nullptr) &&
      ((8 != B->read(&_maxBits, sizeof(uint64))) ||
       (4 != B->read(&inLen,    sizeof(uint32))) ||
       (4 != B->read(&inMax,    sizeof(uint32)))))
    return(false);

  //  Make sure that maxBits is a multiple of 64.

  _maxBits = roundMaxSizeUp(_maxBits);

  //  Increase the number of blocks we have, if needed.

  resizeArray(_blocks, _blocksMax, _blocksMax, inLen, _raAct::copyData | _raAct::clearNew);

  //  Load the data.  The in-memory structure changed after the disk format
  //  was fixed, hence the convoluted load here.

  uint64  *ua = new uint64 [inLen];

  if (F)   merylutil::loadFromFile(ua,  "dataBlockBgn", inLen, F);
  if (B)   B->read(ua, sizeof(uint64) * inLen);
  for (uint32 ii=0; ii<inLen; ii++)
    _blocks[ii]._bgn = ua[ii];

  if (F)   merylutil::loadFromFile(ua,  "dataBlockLen", inLen, F);
  if (B)   B->read(ua, sizeof(uint64) * inLen);
  for (uint32 ii=0; ii<inLen; ii++)
    _blocks[ii]._len = ua[ii];

  delete [] ua;

  //  Load the data, for real.  Any words in the block that aren't read are cleared.

  for (uint32 ii=0; ii<inLen; ii++) {
    uint64  nWordsToRead  = bitsToWords(_blocks[ii]._len);
    uint64  nWordsToClear = _blocks[ii]._max / 64 - nWordsToRead;

    assert(_blocks[ii]._max % 64 == 0);

    if (_blocks[ii]._len > _blocks[ii]._max) {
      delete [] _blocks[ii]._dat;

      _blocks[ii]._max = nWordsToRead * 64;
      _blocks[ii]._dat = new uint64 [_blocks[ii]._max / 64];

      nWordsToClear = 0;
    }

    if (F)   merylutil::loadFromFile(_blocks[ii]._dat, "dataBlocks", nWordsToRead, F);
    if (B)   B->read(_blocks[ii]._dat, sizeof(uint64) * nWordsToRead);

    memset(_blocks[ii]._dat + nWordsToRead, 0, nWordsToClear);
  }

  if (inLen == 0)        //  If the input is empty, no blocks are loaded.
    allocateBlock();     //  Initialize one just so we have something.

  setPosition(0);

  return(true);
}




//  Set the position of stuffedBits to 'position'.  If that position
//  doesn't exist, position is set to the end of the data.
//
//  Terminates when
//    position isn't in the array
void
stuffedBits::setPosition(uint64 position) {

  _dataBlk = 0;

  //  While
  //    Still a valid block,
  //    Block ends before position,
  //    Block is not empty:
  //      Move to the next block.
  //
  while ((_dataBlk < _blocksMax) &&
         (_blocks[_dataBlk]._bgn + _blocks[_dataBlk]._len < position) &&
         (_blocks[_dataBlk]._len > 0))
    _dataBlk++;

  //  If there are no blocks, clear our position.  This really shouldn't
  //  occur, but did before load() above allocated a block for completely
  //  empty inputs.
  if      (_blocksMax == 0) {
    fprintf(stderr, "WARNING: setPosition() called on empty stuffedBits.\n");
    _dataPos = 0;
    _data    = nullptr;
    _dataWrd = 0;
    _dataBit = 64;
    return;
  }

  //  If we've run off the end of the array, all we can do is set the
  //  position to the actual end.
  //
  if (_dataBlk == _blocksMax) {
    assert(_dataBlk > 0);
    _dataPos = _blocks[--_dataBlk]._len;
  }

  //  Otherwise, we're in valid data.  If we're within the contents of the
  //  block, set to that position.
  else if (position <= _blocks[_dataBlk]._bgn + _blocks[_dataBlk]._len) {
    _dataPos = position - _blocks[_dataBlk]._bgn;
  }

  //  Otherwise, the position doesn't exist in stuffedBits and we also
  //  set it to the end of the data.
  else {
    if (_dataBlk > 0)  //  An empty stuffedBits will terminate the loop at
      _dataBlk--;      //  _dataBlk == 0 because _len of that block is zero.

    _dataPos = _blocks[_dataBlk]._len;
  }

  //  We've found the correct block and position.  Update
  //  pointers to that spot.

  assert(_dataBlk <  _blocksMax);
  assert(_dataPos <= _blocks[_dataBlk]._len);

  _data    = _blocks[_dataBlk]._dat;

  _dataWrd =      _dataPos / 64;
  _dataBit = 64 - _dataPos % 64;
}


uint64
stuffedBits::getPosition(void) {
  if (_dataBlk < _blocksMax)
    return(_blocks[_dataBlk]._bgn + _dataPos);
  return(0);
}


uint64
stuffedBits::getLength(void) {
  uint64 nBits = 0;

  for (uint32 ii=0; ((ii < _blocksMax) &&
                     (_blocks[ii]._len > 0)); ii++)
    nBits += _blocks[ii]._len;

  return(nBits);
}

}  //  namespace merylutil::bits::v1
