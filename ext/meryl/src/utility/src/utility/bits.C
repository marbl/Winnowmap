
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

#include "bits.H"
#include "files.H"
#include "system.H"


//  Set up the Fibonacci encoding table.
//
//  It takes 46 values to saturate a uint32 (fib[47] > uint32).
//  It takes 92 values to saturate a uint64 (fib[93] > uint64).
//
//  This is NOT an efficient method; it should be used only to get a
//  compile-time evaluatable function for populating fibonaciNumber[].
//
constexpr
uint64
__fibNum(uint32 f) {
  uint64  fs[93] = { 1, 1, 0 };

  for (uint32 ii=2; ii<=f; ii++)
    fs[ii] = fs[ii-1] + fs[ii-2];

  return(fs[f]);
}  //  C++14 required!

uint64
__fibNumV[93] = { __fibNum( 0), __fibNum( 1), __fibNum( 2), __fibNum( 3), __fibNum( 4), __fibNum( 5), __fibNum( 6), __fibNum( 7), __fibNum( 8), __fibNum( 9),
                  __fibNum(10), __fibNum(11), __fibNum(12), __fibNum(13), __fibNum(14), __fibNum(15), __fibNum(16), __fibNum(17), __fibNum(18), __fibNum(19),
                  __fibNum(20), __fibNum(21), __fibNum(22), __fibNum(23), __fibNum(24), __fibNum(25), __fibNum(26), __fibNum(27), __fibNum(28), __fibNum(29),
                  __fibNum(30), __fibNum(31), __fibNum(32), __fibNum(33), __fibNum(34), __fibNum(35), __fibNum(36), __fibNum(37), __fibNum(38), __fibNum(39),
                  __fibNum(40), __fibNum(41), __fibNum(42), __fibNum(43), __fibNum(44), __fibNum(45), __fibNum(46), __fibNum(47), __fibNum(48), __fibNum(49),
                  __fibNum(50), __fibNum(51), __fibNum(52), __fibNum(53), __fibNum(54), __fibNum(55), __fibNum(56), __fibNum(57), __fibNum(58), __fibNum(59),
                  __fibNum(60), __fibNum(61), __fibNum(62), __fibNum(63), __fibNum(64), __fibNum(65), __fibNum(66), __fibNum(67), __fibNum(68), __fibNum(69),
                  __fibNum(70), __fibNum(71), __fibNum(72), __fibNum(73), __fibNum(74), __fibNum(75), __fibNum(76), __fibNum(77), __fibNum(78), __fibNum(79),
                  __fibNum(80), __fibNum(81), __fibNum(82), __fibNum(83), __fibNum(84), __fibNum(85), __fibNum(86), __fibNum(87), __fibNum(88), __fibNum(89),
                  __fibNum(90), __fibNum(91), __fibNum(92) };

static_assert(__fibNum(45) < __fibNum(46));  //  Fails if 32-bit signed
static_assert(__fibNum(46) < __fibNum(47));  //  Fails if 32-bit unsigned
static_assert(__fibNum(91) < __fibNum(92));  //  Fails if 64-bit signed
//atic_assert(__fibNum(92) < __fibNum(93));  //  Always fails.






inline
uint64
bitsToWords(uint64 bits) {
  return(bits / 64 + ((bits % 64) ? 1 : 0));
}




stuffedBits::stuffedBits(uint64 nBits) {

  _maxBits = roundMaxSizeUp(nBits);

  allocateBlock();
}



stuffedBits::stuffedBits(char const *inputName) {
  FILE *inFile = AS_UTL_openInputFile(inputName);
  load(inFile, nullptr);
  AS_UTL_closeFile(inFile);
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
      ((1 != ::loadFromFile(_maxBits, "dataBlockLenMax", F, false)) ||
       (1 != ::loadFromFile(inLen,    "blocksLen",       F, false)) ||
       (1 != ::loadFromFile(inMax,    "blocksMax",       F, false))))
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

  if (F)   ::loadFromFile(ua,  "dataBlockBgn", inLen, F);
  if (B)   B->read(ua, sizeof(uint64) * inLen);
  for (uint32 ii=0; ii<inLen; ii++)
    _blocks[ii]._bgn = ua[ii];

  if (F)   ::loadFromFile(ua,  "dataBlockLen", inLen, F);
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

    if (F)   ::loadFromFile(_blocks[ii]._dat, "dataBlocks", nWordsToRead, F);
    if (B)   B->read(_blocks[ii]._dat, sizeof(uint64) * nWordsToRead);

    memset(_blocks[ii]._dat + nWordsToRead, 0, nWordsToClear);
  }

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

  //  If we've run off the end of the array, all we can do is set the
  //  position to the actual end.
  //
  if (_dataBlk == _blocksMax) {
    assert(_dataBlk > 0);
    //_dataBlk = _blocksMax-1;
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



//  A special case of getBinary().
bool
stuffedBits::getBit(void) {

  moveToNextBlock(1);

  bool value = saveRightBits(_data[_dataWrd] >> (_dataBit - 1), 1);

  _dataPos++;
  _dataBit--;

  updateBit();

  return(value);
}



bool
stuffedBits::testBit(void) {

  moveToNextBlock(1);

  bool value = saveRightBits(_data[_dataWrd] >> (_dataBit - 1), 1);

  return(value);
}



//  A special case of setBinary().
void
stuffedBits::setBit(bool bit) {

  ensureSpaceInCurrentBlock(1);

  uint64  mask = ((uint64)1 << (_dataBit - 1));

  if (bit)
    _data[_dataWrd] |=  mask;
  else
    _data[_dataWrd] &= ~mask;

  _dataPos++;
  _dataBit--;

  updateBit();
  updateLen();
}



uint64
stuffedBits::getUnary(void) {
  uint64  value = 0;
  uint64  wrd;

  //  Ensure we're in valid data.

  moveToNextBlock(1);

  //  Word align us first.

  wrd = _data[_dataWrd] << (64 - _dataBit);

  //  Skip entire words.  For the first word, if we're left with only zeros
  //  after the shifting, we increase the 'value' by the number of bits
  //  we could read in the word, not the full 64 bits that are zero.

  while (wrd == 0) {
    value    += _dataBit;

    _dataPos += _dataBit;
    _dataWrd += 1;

    _dataBit  = 64;
    wrd       = _data[_dataWrd];
  }

  //  Decode the last partial word.

  wrd       = 64 - countNumberOfBits64(wrd);

  value    += wrd;

  _dataPos += wrd + 1;
  _dataBit -= wrd + 1;

  updateBit();

  return(value);
}

uint64 *
stuffedBits::getUnary(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getUnary();

  return(values);
}




uint32
stuffedBits::setUnary(uint64 value) {

  ensureSpaceInCurrentBlock(value+1);

  //  If we fit entirely within this word, handle it special.

  if (value + 1 < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + value + 1);

    _dataPos += value + 1;
    _dataBit -= value + 1;

    _data[_dataWrd] |= (uint64)1 << _dataBit;

    updateLen();

    return(value + 1);
  }

  //  We fit _exactly_ in this word, special again!

  if (value + 1 == _dataBit) {
    _data[_dataWrd]  = clearRightBits(_data[_dataWrd], _dataBit);
    _data[_dataWrd] |= 1;   //  ALWAYS the last bit.

    _dataPos += value + 1;  //  ALWAYS move to the next word.
    _dataWrd += 1;
    _dataBit  = 64;

    updateLen();

    return(value + 1);
  }

  //  Otherwise, we span at least two words.  First, get us word aligned,
  //  by clearing the rest of this word.

  assert(value + 1 > _dataBit);

  _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit);

  value    -= _dataBit;

  _dataPos += _dataBit;
  _dataWrd += 1;
  _dataBit  = 64;

  assert((_dataPos % 64) == 0);

  //  Then, set entire words to zero.

  while (value >= 64) {
    _data[_dataWrd] = 0;

    value    -= 64;

    _dataPos += 64;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Finally, set the last partial word.  This is similar to the cases
  //  at the start, but we know the bits will always be on the left of the word.

  _data[_dataWrd] = clearLeftBits(_data[_dataWrd], value + 1);

  _dataPos += value + 1;                      //  Skip the zero bits.
  _dataBit -= value + 1;                      //  (and the sentinel)

  _data[_dataWrd] |= (uint64)1 << _dataBit;   //  Add the sentinel.

  //  Update the pointers.

  updateLen();
  updateBit();

  return(value + 1);
}



uint32
stuffedBits::setUnary(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setUnary(values[ii]);

  return(size);
}





uint64
stuffedBits::getBinary(uint32 width) {
  uint64  value = 0;

  if (width == 0)
    return(0);

  assert(width < 65);

  //  Ensure we're in valid data.

  moveToNextBlock(width);

  //  If we're contained in a single word, special case.

  if (width < _dataBit) {
    value = saveRightBits(_data[_dataWrd] >> (_dataBit - width), width);

    _dataPos += width;
    _dataBit -= width;
  }

  //  If we're exactly in a single word, another special case.

  else if (width == _dataBit) {
    value = saveRightBits(_data[_dataWrd], width);

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we're spanning two words.

  else {
    uint64  w1 =         _dataBit;
    uint64  w2 = width - _dataBit;

    uint64  l  = saveRightBits(_data[_dataWrd + 0], w1) <<       w2;
    uint64  r  = saveLeftBits (_data[_dataWrd + 1], w2) >> (64 - w2);

    value = l | r;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  return(value);
}



uint64 *
stuffedBits::getBinary(uint32 width, uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getBinary(width);

  return(values);
}




uint32
stuffedBits::setBinary(uint32 width, uint64 value) {

  if (width == 0)
    return(0);

  assert(width < 65);

  ensureSpaceInCurrentBlock(width);

  //  Mask off any pieces we're not supposed to be seeing.

  value = saveRightBits(value, width);

  //  If we fit entirely within this word, handle it special.

  if (width < _dataBit) {
    _data[_dataWrd] = clearMiddleBits(_data[_dataWrd], 64 - _dataBit, _dataBit + width) | (value << (_dataBit - width));

    _dataPos += width;
    _dataBit -= width;
  } else

  //  We fit _exactly_ in this word, special again!

  if (width == _dataBit) {
    _data[_dataWrd] = clearRightBits(_data[_dataWrd], _dataBit) | value;

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64;
  }

  //  Otherwise, we span two words.

  else {
    uint32  w1 =         _dataBit;
    uint32  w2 = width - _dataBit;

    _data[_dataWrd + 0] = clearRightBits(_data[_dataWrd + 0], w1) | (value >> (     w2));
    _data[_dataWrd + 1] = clearLeftBits (_data[_dataWrd + 1], w2) | (value << (64 - w2));

    _dataPos += width;
    _dataWrd += 1;
    _dataBit  = 64 - w2;
  }

  //  updateBit() isn't needed; it's handled in the special cases.

  updateLen();

  return(width);
}



uint32
stuffedBits::setBinary(uint32 width, uint64 number, uint64 *values) {
  uint32 size = 0;

  if (width == 0)
    return(0);

  for (uint64 ii=0; ii<number; ii++)
    size += setBinary(width, values[ii]);

  return(size);
}





////////////////////////////////////////
//  ELIAS GAMMA CODED DATA
//
//  Unary coded length of binary data, then binary data of that length.
//  Works only on positive (non-zero) integers.
//
uint64
stuffedBits::getEliasGamma(void) {
  uint32  N = getUnary();
  uint64  V = getBinary(N);

  V |= (uint64)1 << N;

  return(V);
}



uint64 *
stuffedBits::getEliasGamma(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getEliasGamma();

  return(values);
}



uint32
stuffedBits::setEliasGamma(uint64 value) {
  uint32 size = 0;
  uint64 N    = countNumberOfBits64(value) - 1;

  assert(value > 0);

  size += setUnary(N);
  size += setBinary(N, value);

  return(size);
}



uint32
stuffedBits::setEliasGamma(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setEliasGamma(values[ii]);

  return(size);
}



////////////////////////////////////////
//  ELIAS DELTA CODED DATA
//
//  Similar to the gamma code, except the number of bits itself
//  is gamma coded.  An optimization can drop the high order bit (it's always 1)
//  from the binary coded data.  We don't do that.
//
uint64
stuffedBits::getEliasDelta(void) {
  uint32  N = getEliasGamma() - 1;
  uint64  V = getBinary(N);

  V |= (uint64)1 << N;

  return(V);
}



uint64 *
stuffedBits::getEliasDelta(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getEliasDelta();

  return(values);
}




uint32
stuffedBits::setEliasDelta(uint64 value) {
  uint32 size = 0;
  uint32 N    = countNumberOfBits64(value);

  assert(value > 0);

  size += setEliasGamma(N);
  size += setBinary(N-1, value);

  return(size);
}



uint32
stuffedBits::setEliasDelta(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setEliasDelta(values[ii]);

  return(size);
}





////////////////////////////////////////
//  FIBONACCI CODED DATA


uint64
stuffedBits::getZeckendorf(void) {
  uint64  value = 0;
  uint32  ff    = 1;

  //  The first bit in the official representation, representing the
  //  redundant value 1, is always zero, and we don't save it.  Thus, start
  //  decoding at ff=1.

  bool tbit = getBit();
  bool nbit = getBit();

  while (true) {
    if (tbit)
      value += __fibNumV[ff];

    if (tbit && nbit)
      break;

    ff++;

    tbit = nbit;
    nbit = getBit();
  }

  return(value);
}



uint64 *
stuffedBits::getZeckendorf(uint64 number, uint64 *values) {

  if (values == NULL)
    values = new uint64 [number];

  for (uint64 ii=0; ii<number; ii++)
    values[ii] = getZeckendorf();

  return(values);
}




uint32
stuffedBits::setZeckendorf(uint64 value) {
  uint32  ff = 0;

  uint64  word1 = 0;   uint32  wlen1 = 0;
  uint64  word2 = 0;   uint32  wlen2 = 0;

  //  Find the largest Fibonacci number smaller than our value.
  //  Probably should be binary searching for this.

  while ((ff < 93) && (__fibNumV[ff] <= value))
    ff++;

  //  For each smaller Fibonacci number:
  //    If the Fibonacci number is more than the value, it's not used in the
  //    encoding.  Push on a zero.
  //
  //    Otherwise, it is used in the encoding.  Push on a 1, and remove the
  //    fib number from our value.
  //
  while (ff-- > 0) {
    word2 <<= 1;                       //  Make space for the new bit.

    if (__fibNumV[ff] <= value) {     //  If used in the encoding,
      word2 |= 1;                     //  set the bit and remove
      value -= __fibNumV[ff];         //  it from the value.
    }

    if (++wlen2 > 60) {               //  If we're running outta
      word1 = word2;                  //  bits in the word, save it
      wlen1 = wlen2;                  //  to the first word to output.
      wlen2 = 0;                      //
    }
  }

  //  Reverse the words so we see the low bit first, then push on a
  //  terminating 1 so we end the string with a pair of 1 bits.
  //
  //  The lower bits, in word2, can have the (post-reverse) left-most bit,
  //  representing the redundant 1, stripped off.
  //
  //  An annoying special case oocurs when there are exactly 60 bits in the
  //  encoding: word2 is now empty!

  if (wlen1 == 0) {
    word2 = reverseBits64(word2);

    word2 >>= (64 - wlen2);

    word2 <<= 1;
    word2  |= 1;
    wlen2  += 1;

    wlen2  -= 1;  //  Strip off left-most bit.  Go Optimizer, Go!

    setBinary(wlen2, word2);
  }

  else if (wlen2 == 0) {
    word1 = reverseBits64(word1);

    word1 >>= (64 - wlen1);

    word1 <<= 1;
    word1  |= 1;
    wlen1  += 1;

    wlen1  -= 1;  //  Strip off left-most bit.  Go Optimizer, Go!

    setBinary(wlen1, word1);
  }

  else {
    word2 = reverseBits64(word2);
    word1 = reverseBits64(word1);

    word2 >>= (64 - wlen2);
    word1 >>= (64 - wlen1);

    word1 <<= 1;
    word1  |= 1;
    wlen1  += 1;

    wlen2  -= 1;

    setBinary(wlen2, word2);
    setBinary(wlen1, word1);
  }


  return(wlen1 + wlen2);
}



uint32
stuffedBits::setZeckendorf(uint64 number, uint64 *values) {
  uint32  size = 0;

  for (uint64 ii=0; ii<number; ii++)
    size += setZeckendorf(values[ii]);

  return(size);
}

