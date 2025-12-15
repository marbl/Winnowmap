
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

#ifndef MERYLUTIL_BITS_V1_STUFFEDBITS_H
#define MERYLUTIL_BITS_V1_STUFFEDBITS_H

#include "types.H"
#include "arrays.H"
#include "files.H"
#include "system.H"

namespace merylutil::inline bits::inline v1 {

class stuffedBits {
public:
  stuffedBits(uint64 nBits=16 * 1024 * 1024 * 8);
  stuffedBits(const char *inputName);
  stuffedBits(FILE *inFile);
  stuffedBits(readBuffer *B);
  ~stuffedBits();

  //  Debugging.

  char    *displayWord(uint64 w) {
    static char dd[129];
    return(toBin(_data[w], dd));
  };

  //  Files.

  void     dump(FILE *F, writeBuffer *B);
  void     dumpToBuffer(writeBuffer *B)   {        dump(nullptr, B);  }
  void     dumpToFile(FILE *F)            {        dump(F, nullptr);  }

  bool     load(FILE *F, readBuffer *B);
  bool     loadFromBuffer(readBuffer *B)  { return(load(nullptr, B)); }
  bool     loadFromFile(FILE *F)          { return(load(F, nullptr)); }

  //  Management of the read/write head.

  void     setPosition(uint64 position);
  uint64   getPosition(void);

  uint64   getLength(void);

  void     byteAlign(void);

  //  SINGLE BITS

  bool     getBit(void);           //  get a bit.
  bool     testBit(void);          //  get a bit, but don't move position.
  void     setBit(bool on=true);   //  set a bit.

  //  UNARY CODED DATA

  uint64   getUnary(void);
  uint64  *getUnary(uint64 number, uint64 *values);

  uint32   setUnary(uint64 value);
  uint32   setUnary(uint64 number, uint64 *values);

  //  BINARY CODED DATA

  uint64   getBinary(uint32 width);
  uint64  *getBinary(uint32 width, uint64 number, uint64 *values=NULL);

  uint32   setBinary(uint32 width, uint64 value);
  uint32   setBinary(uint32 width, uint64 number, uint64 *values);

  //  ELIAS GAMMA CODED DATA
  //
  //    Unary coded length of binary data, then binary data of that length.
  //    Works only on positive (non-zero) integers.

  uint64   getEliasGamma(void);
  uint64  *getEliasGamma(uint64 number, uint64 *values=NULL);

  uint32   setEliasGamma(uint64 value);
  uint32   setEliasGamma(uint64 number, uint64 *values);

  //  ELIAS DELTA CODED DATA
  //
  //    Similar to the gamma code, except the number of bits itself is gamma
  //    coded.  An optimization can drop the high order bit (it's always 1)
  //    from the binary coded data.  We don't do that.

  uint64   getEliasDelta(void);
  uint64  *getEliasDelta(uint64 number, uint64 *values=NULL);

  uint32   setEliasDelta(uint64 value);
  uint32   setEliasDelta(uint64 number, uint64 *values);

  //  ELIAS OMEGA CODED DATA - the omega code looks hard to implement - the
  //  encoding and decoding streams are backwards from each other.  The idea
  //  is:
  //
  //    push the binary representation of the value onto a stack.
  //    set value to one less than the number of bits emitted last.
  //    loop until the value is 1.
  //
  //  The stream is constructed by emitting the words on the stack, and
  //  terminating the stream with an extra '0'.
  //
#if 0
  uint64   getEliasOmega(void);
  uint64  *getEliasOmega(uint64 number, uint64 *values=NULL);

  uint32   setEliasOmega(uint64 value);
  uint32   setEliasOmega(uint64 number, uint64 *values);
#endif

  //  GOLOMB CODED DATA
  //
  //  Pick m.  For m == power_of_two, this is RICE CODED DATA.
  //
  //    q = floow(n/m).
  //    r = n-qm
  //    c = ceil(log_2 m)
  //
  //  Unary encode q, binary encode r.
  //
  //  The first 2^c-m values are encoded as c-1 bit values, starting with 00...00,
  //  The rest as c-bit numbers, ending with 11...11
  //


  //  FIBONACCI CODED DATA
  //
  //  A Fibonacci number is F(n) = F(n-1) + F(n-2), wher F(0) = F(1) = 1.
  //
  //  The Zeckendorf representation of a number encodes it such that no
  //  two consecurive Fibonacci numbers are used.  From the definition
  //  of a Fibonacci number, any pattern "100" can be replaced with "011".
  //  A number encoded after this transofmration is in the Fibonacci
  //  representation ("Zeckendorf representation" seems to be a real thing,
  //  I just made up "Fibonacci representation" - the two terms seem to
  //  be used synonymously in the real world).
  //
  //  Once encoded, it's added to the bit stream reversed.
  //
  //  For the Zeckendorf representation, a single 1-bit is added, terminating
  //  the number with the last '1' bit of data, followed immediately by
  //  another '1' bit.  (Because, by definition, there are no two adjacent
  //  set bits in the encoded number).
  //
  //  For the Fibonacci representation, we need to append two '0' bits.
  //  (Because, by definition, there are no two adjacent unset bits in the
  //  representation).  BUT, this representation saves at most one bit
  //  (replacing 100 at the start of the string by 011), but the savings
  //  is lost by the extra stop bit we need.
  //
  uint64   getZeckendorf(void);
  uint64  *getZeckendorf(uint64 number, uint64 *values=NULL);

  uint32   setZeckendorf(uint64 value);
  uint32   setZeckendorf(uint64 number, uint64 *values);

  //  Old meryl uses preDecrement() when using compressed bucket counting.
  //  Nothing else in canu uses these, and they're painful, so left unimplemented.
#if 0
  uint64   preIncrementBinary(uint64 width, uint64 position);
  uint64   postIncrementBinary(uint64 width, uint64 position);
  uint64   preDecrementBinary(uint64 width, uint64 position);
  uint64   postDecrementBinary(uint64 width, uint64 position);
#endif


private:
  void     updateLen(void);    //  Update block len after a write.
  void     updateBit(void);    //  Move to next block if we're at the end of the current.

  void     moveToNextBlock(uint64 wordLen);             //  move to next if 'wordLen' can't be read
  void     ensureSpaceInCurrentBlock(uint64 wordLen);   //  move to next if 'wordLen' can't be written

  uint64   roundMaxSizeUp(uint64 numbits);              //  Round up to next multiple of 64.
  void     allocateBlock(void);                         //  Allocate and init a new block, if needed.

  void     eraseBlocks(void);                           //  Resets allocated blocks to size zero.

  struct _dBlock {
    uint64  _bgn = 0;        //  Starting position, in the global file, of this block.
    uint64  _len = 0;        //  Length of the data in this block, in BITS.
    uint64  _max = 0;        //  Allocated size of this block, in WORDS.
    uint64 *_dat = nullptr;  //  Just piles of bits.  Nothing interesting here.
  };

  uint64    _maxBits = 0;             //  Allocated length of each block (in BITS).

  uint32    _blocksMax = 0;           //  Number of blocks we can allocate.
  _dBlock  *_blocks    = nullptr;     //  Blocks!

  uint64    _dataPos    = 0;          //  Position in this block, in BITS.
  uint64   *_data       = nullptr;    //  Pointer to the data in the currently active data block.

  uint32    _dataBlk    = 0;          //  Active data block.
  uint64    _dataWrd    = 0;          //  Active word in the active data block.
  uint64    _dataBit    = 64;         //  Active bit in the active word in the active data block (aka, number of bits left in this word)
};



//  For writing, update the length of the block to the maximum of where
//  we're at now and the existing length.
//
inline
void
stuffedBits::updateLen(void) {
  _blocks[_dataBlk]._len = std::max(_dataPos, _blocks[_dataBlk]._len);
}


//  For both reading and writing, move to the next word if we're at the end
//  of the current one.
//
inline
void
stuffedBits::updateBit(void) {
  if (_dataBit == 0) {
    _dataWrd += 1;
    _dataBit  = 64;
  }
}


//  Ensure that a read of the next sub-word of length readLength is present
//  entirely in the current block.  Move to the next block if not.
//
inline
void
stuffedBits::moveToNextBlock(uint64 wordLen) {
  assert(_dataBit != 0);
  assert(_dataBit <= 64);
  assert(_blocksMax > 0);
  assert(_dataBlk < _blocksMax);

  //  The word is in this block, we need to do nothing.

  if (_dataPos + wordLen <= _blocks[_dataBlk]._len)
    return;

  //  If we're not at the end of the current block, something is amiss;
  //  words do not span blocks.

  if (_dataPos != _blocks[_dataBlk]._len)
    fprintf(stderr, "stuffedBits()--  ERROR: _dataPos=" F_U64 " != _blocks[%u]._len=" F_U64 "\n",
            _dataPos, _dataBlk, _blocks[_dataBlk]._len);
  assert(_dataPos == _blocks[_dataBlk]._len);

  //  Move to the next block, failing if there are no more blocks with data

  _dataBlk += 1;

  if (_dataBlk >= _blocksMax)
    fprintf(stderr, "ERROR: _dataBlk = %u  _blocksMax = %u\n", _dataBlk, _blocksMax);
  assert(_dataBlk < _blocksMax);

  if (_blocks[_dataBlk]._len == 0)
    fprintf(stderr, "stuffedBits()-- ERROR: no data in block at _dataBlk=%u\n", _dataBlk);
  assert(_blocks[_dataBlk]._len != 0);

  if (_blocks[_dataBlk]._dat == nullptr)
    fprintf(stderr, "stuffedBits()-- ERROR: no data block at _dataBlk=%u\n", _dataBlk);
  assert(_blocks[_dataBlk]._dat != nullptr);

  //  Reset the various pointers to the start of the current block.

  _dataPos  = 0;
  _data     = _blocks[_dataBlk]._dat;

  _dataWrd  = 0;
  _dataBit  = 64;
}


//  Ensure that a write of a sub-word of length wordLen will occur entirely
//  in the current block.  Move to the next block if not.
//
//  Assumes the current block exists.
//
inline
void
stuffedBits::ensureSpaceInCurrentBlock(uint64 wordLen) {
  assert(_dataBit != 0);
  assert(_dataBit <= 64);
  assert(_blocksMax > 0);
  assert(_dataBlk < _blocksMax);

  //  The word will fit in this block, we need to do nothing.

  if (_dataPos + wordLen <= _blocks[_dataBlk]._max)
    return;

  //  Otherwise, there isn't enough space in the current block for a write
  //  of 'wordLen' bits.  Terminate the current block and move to the next.

  _blocks[_dataBlk]._len = _dataPos;

  _dataBlk++;

  allocateBlock();
}


//  The allocated blocks (_dat and _max) need to be a multiple of 64 so we
//  can use a simple shift to get the number of words to allocate.  This function
//  will round up to the next multiple, with a special case for zero.
//
inline
uint64
stuffedBits::roundMaxSizeUp(uint64 numbits) {

  if (numbits == 0)                           //  If zero, set it to a few pages minus a
    numbits = (8 * getPageSize() - 32) * 8;   //  few pointers for the allocator to use.

  if ((numbits & 0x3fllu) != 0)               //  If not a multiple of 64, bump
    numbits = (numbits | 0x3fllu) + 1;        //  up to the next multiple.

  assert(numbits % 64 == 0);                  //  Be paranoid.  Be very paranoid.

  return(numbits);
}


//  Allocate a new block and initialize it, if needed.
//
inline
void
stuffedBits::allocateBlock(void) {

  //  Allocate another 32 blocks if _dataBlk >= _blocksMax.
  increaseArray(_blocks, _dataBlk, _blocksMax, 32, _raAct::copyData | _raAct::clearNew);

  //  Initialize the position and length of the block.
  _blocks[_dataBlk]._bgn = (_dataBlk == 0) ? 0 : _blocks[_dataBlk - 1]._bgn + _blocks[_dataBlk - 1]._len;
  _blocks[_dataBlk]._len =                   0;

  //  Allocate space for the data, if needed.
  if (_blocks[_dataBlk]._dat == nullptr) {
    _blocks[_dataBlk]._max = _maxBits;
    _blocks[_dataBlk]._dat = new uint64 [ _blocks[_dataBlk]._max / 64 ];
  }

  assert(_blocks[_dataBlk]._max % 64 == 0);

  //  Clear it.
  memset(_blocks[_dataBlk]._dat, 0, sizeof(uint64) * _blocks[_dataBlk]._max / 64);

  //  Set the various pointers to the start of the current block.
  _dataPos  = 0;
  _data     = _blocks[_dataBlk]._dat;

  assert(_data != nullptr);

  _dataWrd  = 0;
  _dataBit  = 64;
}


//  Resets any allocated blocks to have size zero.
//
inline
void
stuffedBits::eraseBlocks(void) {

  if (_blocksMax == 0)
    return;

  for (uint32 ii=0; (ii < _blocksMax); ii++)
    _blocks[ii]._len = 0;

  //  Set the various pointers to the start of the current block.
  _dataPos  = 0;
  _data     = _blocks[_dataBlk]._dat;

  assert(_data != nullptr);

  _dataWrd  = 0;
  _dataBit  = 64;
}

}  //  namespace merylutil::bits::v1

#endif  //  MERYLUTIL_BITS_V1_STUFFEDBITS_H
