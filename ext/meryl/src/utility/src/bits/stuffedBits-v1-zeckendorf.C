
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
#include "stuffedBits-v1.H"

//  stuffedBits operations on Zeckendorf (Fibonacci) coded data.

namespace merylutil::inline bits::inline v1 {

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

}  //  namespace merylutil::bits::v1
