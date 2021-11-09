
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
#include "strings.H"
#include "mt19937ar.H"
#include "bitsTest-MaskData.H"

#include <algorithm>

char          b1[65];
char          b2[65];
char          b3[65];

void
testMasks(bool verbose) {
  uint128  m128;
  uint64   m64;
  uint32   m32;
  uint16   m16;
  uint8    m8;

  for (uint32 ii=0; ii<=128; ii++) {
    if ((strcmp(maskstr[ii][0], toHex(buildLowBitMask<uint128>(ii))) != 0) || (strcmp(maskstr[ii][1], toHex(buildHighBitMask<uint128>(ii))) != 0) ||
        (strcmp(maskstr[ii][2], toHex(buildLowBitMask<uint64> (ii))) != 0) || (strcmp(maskstr[ii][3], toHex(buildHighBitMask<uint64> (ii))) != 0) ||
        (strcmp(maskstr[ii][4], toHex(buildLowBitMask<uint32> (ii))) != 0) || (strcmp(maskstr[ii][5], toHex(buildHighBitMask<uint32> (ii))) != 0) ||
        (strcmp(maskstr[ii][6], toHex(buildLowBitMask<uint16> (ii))) != 0) || (strcmp(maskstr[ii][7], toHex(buildHighBitMask<uint16> (ii))) != 0) ||
        (strcmp(maskstr[ii][8], toHex(buildLowBitMask<uint8>  (ii))) != 0) || (strcmp(maskstr[ii][9], toHex(buildHighBitMask<uint8>  (ii))) != 0)) {
      fprintf(stderr, "Failed at ii=%u\n", ii);
      fprintf(stderr, "Expected: %s %s  %s %s  %s %s  %s %s  %s %s\n",
              maskstr[ii][0], maskstr[ii][1],
              maskstr[ii][2], maskstr[ii][3],
              maskstr[ii][4], maskstr[ii][5],
              maskstr[ii][6], maskstr[ii][7],
              maskstr[ii][8], maskstr[ii][9]);
      fprintf(stderr, "Computed: %s %s  %s %s  %s %s  %s %s  %s %s\n",
              toHex(buildLowBitMask<uint128>(ii)), toHex(buildHighBitMask<uint128>(ii)),
              toHex(buildLowBitMask<uint64> (ii)), toHex(buildHighBitMask<uint64> (ii)),
              toHex(buildLowBitMask<uint32> (ii)), toHex(buildHighBitMask<uint32> (ii)),
              toHex(buildLowBitMask<uint16> (ii)), toHex(buildHighBitMask<uint16> (ii)),
              toHex(buildLowBitMask<uint8>  (ii)), toHex(buildHighBitMask<uint8>  (ii)));
    }
  }

  if (verbose) {
    fprintf(stdout, "\n");
    fprintf(stdout, "/******************************************************************************\n");
    fprintf(stdout, " *\n");
    fprintf(stdout, " *  This file is part of meryl-utility, a collection of miscellaneous code\n");
    fprintf(stdout, " *  used by Meryl, Canu and others.\n");
    fprintf(stdout, " *\n");
    fprintf(stdout, " *  This software is based on:\n");
    fprintf(stdout, " *    'Canu' v2.0              (https://github.com/marbl/canu)\n");
    fprintf(stdout, " *  which is based on:\n");
    fprintf(stdout, " *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)\n");
    fprintf(stdout, " *    the 'kmer package' r1994 (http://kmer.sourceforge.net)\n");
    fprintf(stdout, " *\n");
    fprintf(stdout, " *  Except as indicated otherwise, this is a 'United States Government Work',\n");
    fprintf(stdout, " *  and is released in the public domain.\n");
    fprintf(stdout, " *\n");
    fprintf(stdout, " *  File 'README.licenses' in the root directory of this distribution\n");
    fprintf(stdout, " *  contains full conditions and disclaimers.\n");
    fprintf(stdout, " */\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "//  This is the output of bitsTest -verbose.  It is (supposedly) the correct\n");
    fprintf(stdout, "//  output of the mask building functions.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "char const *maskstr[129][10] = {\n");
    for (uint32 ii=0; ii<=128; ii++)
      fprintf(stdout, "  { \"%s\", \"%s\",  \"%s\", \"%s\",  \"%s\", \"%s\",  \"%s\", \"%s\",  \"%s\", \"%s\" },\n",
              toHex(buildLowBitMask<uint128>(ii)), toHex(buildHighBitMask<uint128>(ii)),
              toHex(buildLowBitMask<uint64> (ii)), toHex(buildHighBitMask<uint64> (ii)),
              toHex(buildLowBitMask<uint32> (ii)), toHex(buildHighBitMask<uint32> (ii)),
              toHex(buildLowBitMask<uint16> (ii)), toHex(buildHighBitMask<uint16> (ii)),
              toHex(buildLowBitMask<uint8>  (ii)), toHex(buildHighBitMask<uint8>  (ii)));
    fprintf(stdout, "};\n");
  }
}



void
testLogBaseTwo(bool verbose) {
  assert(0 == countNumberOfBits64(0));
  assert(1 == countNumberOfBits64(1));
  assert(2 == countNumberOfBits64(2));
  assert(2 == countNumberOfBits64(3));
  assert(3 == countNumberOfBits64(4));

  for (uint64 ii=2, val=2; val > 0; ii++, val <<= 1) {
    assert(ii-1 == countNumberOfBits64(val-1));
    assert(ii   == countNumberOfBits64(val));
    assert(ii   == countNumberOfBits64(val+1));
  }

  assert(64 == countNumberOfBits64(uint64max - 1));
  assert(64 == countNumberOfBits64(uint64max));
  assert(0  == countNumberOfBits64(uint64max + 1));
}



void
testSaveClear(bool verbose) {
  uint64  bb, rr, rl;

  bb = uint64max;
  rr = uint64max;
  rl = uint64max;
  for (uint32 ii=0; ii<65; ii++) {
    if (verbose)
      fprintf(stderr, "%2u clearLeft %s clearRight %s\n", ii,
              displayWord(clearLeftBits (bb, ii), b1),
              displayWord(clearRightBits(bb, ii), b2));

    assert(clearLeftBits (bb, ii) == rr);   rr >>= 1;
    assert(clearRightBits(bb, ii) == rl);   rl <<= 1;

    assert(clearLeftBits (bb, ii) == buildLowBitMask<uint64>(64 - ii));
    assert(clearRightBits(bb, ii) == buildHighBitMask<uint64>(64 - ii));
  }

  bb = uint64max;
  rr = uint64max;
  rl = uint64max;
  for (uint32 ii=0; ii<65; ii++) {
    if (verbose)
      fprintf(stderr, "%2u  saveLeft %s  saveRight %s\n", ii,
              displayWord(saveLeftBits (bb, ii), b1),
              displayWord(saveRightBits(bb, ii), b2));

    assert(saveLeftBits (bb, ii) == ~rr);   rr >>= 1;
    assert(saveRightBits(bb, ii) == ~rl);   rl <<= 1;

    assert(saveLeftBits (bb, ii) == buildHighBitMask<uint64>(ii));
    assert(saveRightBits(bb, ii) == buildLowBitMask<uint64>(ii));
  }

  bb = uint64max;
  rr = uint64max << 10;
  rl = uint64max;
  for (uint32 ii=0; ii<65; ii++) {
    if (verbose)
      fprintf(stderr, "%2u  saveMid  %s clearMid   %s\n", ii,
              displayWord(saveMiddleBits (bb, ii, 10), b1),
              displayWord(clearMiddleBits(bb, ii, 10), b2));

    assert(saveMiddleBits (bb, ii, 10) ==  rr);
    assert(clearMiddleBits(bb, ii, 10) == ~rr);

    rr >>= 1;
    rr  &= 0xfffffffffffffc00llu;

    rl = buildLowBitMask<uint64>(64 - ii) & 0xfffffffffffffc00llu;

    assert(saveMiddleBits (bb, ii, 10) ==  rl);
    assert(clearMiddleBits(bb, ii, 10) == ~rl);
  }

  bb = uint64max;
  rr = uint64max >> 10;
  rl = uint64max;
  for (uint32 ii=0; ii<65; ii++) {
    if (verbose)
      fprintf(stderr, "%2u  saveMid  %s clearMid   %s\n", ii,
              displayWord(saveMiddleBits (bb, 10, ii), b1),
              displayWord(clearMiddleBits(bb, 10, ii), b2));

    assert(saveMiddleBits (bb, 10, ii) ==  rr);
    assert(clearMiddleBits(bb, 10, ii) == ~rr);

    rr <<= 1;
    rr  &= 0x003fffffffffffffllu;

    rl = buildHighBitMask<uint64>(64 - ii) & 0x003fffffffffffffllu;

    assert(saveMiddleBits (bb, 10, ii) ==  rl);
    assert(clearMiddleBits(bb, 10, ii) == ~rl);
  }
}



//  Just a useful report of the fibonacci numbers.
void
testFibonacciNumbers(bool verbose) {
  uint64  _fibDataMax = 93;
  uint64 *_fibData    = new uint64 [_fibDataMax + 1];

  _fibData[0] = 1;
  _fibData[1] = 1;

  for (uint32 ii=2; ii<_fibDataMax; ii++) {
    _fibData[ii] = _fibData[ii-1] + _fibData[ii-2];

    if (verbose)
      fprintf(stderr, "%4u -- %22lu = %22lu + %22lu -- %22lu\n",
              ii, _fibData[ii], _fibData[ii-1], _fibData[ii-2], UINT64_MAX - _fibData[ii]);

    assert(_fibData[ii] == fibonacciNumber(ii));
    assert(_fibData[ii] > _fibData[ii-1]);
  }
}



void
testExpandCompress(bool verbose, uint32 limit) {
  uint64   o = 0x0000000000000003llu;
  uint64   e = 0;
  uint64   c = 0;

  for (uint32 xx=0; xx<limit; xx++) {
    e = expandTo3(o);
    c = compressTo2(e);

    if (c != o)
      fprintf(stdout, "orig 0x%s expanded %022lo compressed 0x%s FAIL\n", toHex(o), e, toHex(c));
    assert(c == o);

    if (verbose)
      fprintf(stdout, "orig 0x%s expanded %021lo compressed 0x%s\n", toHex(o), e, toHex(c));

    o <<= 2;
  }
}



struct tbaDat {
  uint64  key;
  uint64  pos;
};


void
testBitArray(bool verbose, uint64 nbits, uint64 length) {
  mtRandom          mt;
  tbaDat           *list = new tbaDat [nbits];
  bitArray         *ba   = new bitArray(length);

  if (verbose)
    fprintf(stderr, "Testing bitArray of length %lu with %lu bits set.\n", length, nbits);

  //  Generate a list of random positions.

  for (uint64 xx=0; xx<nbits; xx++) {
    list[xx].key = mt.mtRandom64();
    list[xx].pos = mt.mtRandom64() % (length-2) + 1;
  }

  //  Sort by position, remove duplicates by setting the key to an invalid
  //  value; pos must remain valid so we can find multiple duplicates.

  auto byPos = [](tbaDat const &a, tbaDat const &b) { return(a.pos < b.pos); };
  auto byKey = [](tbaDat const &a, tbaDat const &b) { return(a.key < b.key); };

  std::sort(list, list+nbits, byPos);

  for (uint64 xx=1; xx<nbits; xx++)
    if (list[xx].pos == list[xx-1].pos)
      list[xx].key = uint64max;

  //  Sort again by key to randomize the positions.

  std::sort(list, list+nbits, byKey);

  //  Set all those bits, checking that the bits before/after do not change.

  uint64  numset   = 0;
  uint64  numfound = 0;

  for (uint64 xx=0; xx<nbits; xx++) {
    if (list[xx].key == uint64max)
      continue;

    numset++;

    bool  a = ba->getBit(list[xx].pos-1);
    bool  c = ba->getBit(list[xx].pos+1);

    assert(0 == ba->getBit(list[xx].pos));
    ba->setBit(list[xx].pos, 1);
    assert(1 == ba->getBit(list[xx].pos));

    {
      uint64 nf = 0;

      for (uint64 ll=0; ll<length; ll++)
        if (1 == ba->getBit(ll))
          nf++;

      assert(nf == xx+1);
    }

    assert(a == ba->getBit(list[xx].pos-1));
    assert(c == ba->getBit(list[xx].pos+1));

    ba->flipBit(list[xx].pos-1);
    assert(a != ba->getBit(list[xx].pos-1));
    assert(1 == ba->getBit(list[xx].pos));
    assert(c == ba->getBit(list[xx].pos+1));

    ba->flipBit(list[xx].pos+1);
    assert(a != ba->getBit(list[xx].pos-1));
    assert(1 == ba->getBit(list[xx].pos));
    assert(c != ba->getBit(list[xx].pos+1));

    ba->flipBit(list[xx].pos-1);
    assert(a == ba->getBit(list[xx].pos-1));
    assert(1 == ba->getBit(list[xx].pos));
    assert(c != ba->getBit(list[xx].pos+1));

    ba->flipBit(list[xx].pos+1);
    assert(a == ba->getBit(list[xx].pos-1));
    assert(1 == ba->getBit(list[xx].pos));
    assert(c == ba->getBit(list[xx].pos+1));

    {
      uint64 nf = 0;

      for (uint64 ll=0; ll<length; ll++)
        if (1 == ba->getBit(ll))
          nf++;

      assert(nf == xx+1);
    }
  }

  //  Count the number of bits that are set.

  for (uint64 xx=0; xx<length; xx++) {
    if (1 == ba->getBit(xx))
      numfound++;
  }

  assert(numset == numfound);

  //  Toggle all those bits.

  for (uint64 xx=0; xx<nbits; xx++) {
    if (list[xx].key == uint64max)
      continue;

    assert(1 == ba->getBit(list[xx].pos));
    assert(1 == ba->flipBit(list[xx].pos));
    assert(0 == ba->getBit(list[xx].pos));
  }

  //  Check all bits are zero...and set them to one.

  for (uint64 xx=0; xx<length; xx++) {
    assert(0 == ba->getBit(xx));
    ba->setBit(xx, 1);
  }

  //  Clear the array.

  ba->clear();

  //  And check they're all zero again.

  for (uint64 xx=0; xx<length; xx++)
    assert(0 == ba->getBit(xx));

  //  Done.

  delete [] list;
  delete    ba;
}



void
testWordArray(bool verbose, uint64 length, uint32 wordSize, uint32 arraySize) {
  wordArray  *wa      = new wordArray(wordSize, arraySize, false);
  uint64      mask    = buildLowBitMask<uint64>(wordSize);
  uint64      testLen = 111;
  bool        dodump  = false;

  if (verbose)
    fprintf(stderr, "Testing wordArray with %lu words of %u bits, blocks of %u bits.\n", length, wordSize, arraySize);

  //  Set words to something we can check later.

  for (uint32 ii=0; ii<length; ii++)
    wa->set(ii, ii);

  //  And check.

  for (uint32 ii=0; ii<length; ii++)
    assert(wa->get(ii) == (ii & mask));

  delete wa;
}



void
testWordArraySpeed(bool verbose, uint64 length, uint32 wordSize) {
  wordArray  *wa = new wordArray(wordSize * 2, 1024 * 1024 * 1024, false);
  uint64     *a1 = new uint64 [length];
  uint64     *a2 = new uint64 [length];
  uint64     *a3 = new uint64 [length];

  uint64     dat;
  uint64    mask = buildLowBitMask<uint64>(wordSize);

  fprintf(stderr, "Testing wordArray speed vs 2x arrays of length %lu and width %u.\n",
          length, wordSize);

  double  space = (length * wordSize / 64 + length + length) * 8 / 1024 / 1024 / 1024.0;
  double  start = 0;

  mtRandom  mt;

  fprintf(stderr, "Total space required: %.3f GB.\n", space);

  //  Load data into the wordArray.

  fprintf(stderr, "Loading wordArray.\n");
  mt.mtSetSeed(625);
  start = getTime();

  for (uint64 ii=0; ii<length; ii++) {
    dat   = (ii & mask);
    dat <<=  wordSize;
    dat  |= (ii & mask);

    wa->set(ii, dat);
  }
  fprintf(stderr, "  %.3f seconds\n", getTime() - start);

  //  Load data into the arrays.

  fprintf(stderr, "Loading arrays.\n");
  mt.mtSetSeed(625);
  start = getTime();

  for (uint64 ii=0; ii<length; ii++) {
    dat  = (ii & mask);

    a1[ii] = dat;
    a2[ii] = dat;
    a3[ii] = dat;
  }
  fprintf(stderr, "  %.3f seconds\n", getTime() - start);

  //  Retrieve data from the wordArray.

  fprintf(stderr, "Reading wordArray.\n");
  mt.mtSetSeed(625);
  start = getTime();

  for (uint64 ii=0; ii<length; ii++) {
    dat   = (ii & mask);
    dat <<=  wordSize;
    dat  |= (ii & mask);

    assert(dat == wa->get(ii));
  }
  fprintf(stderr, "  %.3f seconds\n", getTime() - start);

  fprintf(stderr, "Reading wordArray randomly.\n");
  mt.mtSetSeed(625);
  start = getTime();

  for (uint64 jj=0; jj<length; jj++) {
    uint64 ii = mt.mtRandom64() % length;

    dat   = (ii & mask);
    dat <<=  wordSize;
    dat  |= (ii & mask);

    assert(dat == wa->get(ii));
  }
  fprintf(stderr, "  %.3f seconds\n", getTime() - start);

  //  Retrieve data from the arrays.

  fprintf(stderr, "Reading arrays.\n");
  mt.mtSetSeed(625);
  start = getTime();

  for (uint64 ii=0; ii<length; ii++) {
    dat  = (ii & mask);

    assert(dat == a1[ii]);
    assert(dat == a2[ii]);
    assert(dat == a3[ii]);
  }
  fprintf(stderr, "  %.3f seconds\n", getTime() - start);

  fprintf(stderr, "Reading arrays randomly.\n");
  mt.mtSetSeed(625);
  start = getTime();

  for (uint64 jj=0; jj<length; jj++) {
    uint64 ii = mt.mtRandom64() % length;

    dat  = (ii & mask);

    assert(dat == a1[ii]);
    assert(dat == a2[ii]);
    assert(dat == a3[ii]);
  }
  fprintf(stderr, "  %.3f seconds\n", getTime() - start);

  fprintf(stderr, "Passed!\n");

  delete    wa;
  delete [] a1;
  delete [] a2;
  delete [] a3;
}



void
testUnary(bool verbose, uint64 length, uint32 maxWidth) {
  uint64      maxN   = length;
  uint64      Nbits  = 0;
  uint32     *random = new uint32 [maxN];
  mtRandom    mt;

  if (verbose)
    fprintf(stderr, "Testing stuffedBits unary encoding with %lu numbers between 1 and %u.\n", length, maxWidth);

  for (uint64 ii=0; ii<maxN; ii++) {
    random[ii]  = mt.mtRandom32() % maxWidth;
    Nbits      += random[ii] + 1;
  }

  if (verbose)
    fprintf(stderr, "Setting  %lu numbers with total length %lu bits.\n", maxN, Nbits);

  stuffedBits *bits = new stuffedBits;

  for (uint64 position=0, ii=0; ii<maxN; ii++) {
    bits->setUnary(random[ii]);

    position += random[ii] + 1;

    assert(bits->getPosition() == position);
  }

  if (verbose) {
    fprintf(stderr, "Testing  %lu numbers with total length %lu bits encoded into %lu bits.\n", maxN, Nbits, bits->getLength());
    fprintf(stderr, "\n");
    fprintf(stderr, "                     value  b  data word\n");
    fprintf(stderr, "--  ---------------------- --  ----------------------------------------------------------------\n");
    for (uint32 ii=0; ii<10; ii++)
      fprintf(stderr, "%02u  %22u     %s\n", ii, random[ii], bits->displayWord(ii));
    fprintf(stderr, "--  ---------------------- --  ----------------------------------------------------------------\n");
    fprintf(stderr, "\n");
  }

  bits->setPosition(0);

  for (uint64 position=0, ii=0; ii<maxN; ii++) {
    uint64  b = bits->getUnary();

    position += random[ii] + 1;

    assert(random[ii] == b);
    assert(bits->getPosition() == position);
  }

  delete    bits;
  delete [] random;
}





void
testBinary(bool verbose, uint64 length, uint32 maxWidth) {
  uint64      maxN     = length;
  uint64      Nbits    = 0;
  uint32     *width    = new uint32 [maxN];
  uint64     *random   = new uint64 [maxN];
  mtRandom    mt;

  if (verbose)
    fprintf(stderr, "Testing stuffedBits binary encoding with %lu numbers of %u bits each.\n", length, maxWidth);

  for (uint64 ii=0; ii<maxN; ii++) {
    width[ii]   = mt.mtRandom32() % maxWidth;
    random[ii]  = mt.mtRandom64() & (((uint64)1 << width[ii]) - 1);
    Nbits      +=  width[ii];
  }

  if (verbose)
    fprintf(stderr, "Setting  %lu numbers with total length %lu bits.\n", maxN, Nbits);

  stuffedBits *bits = new stuffedBits;

  for (uint64 position=0, ii=0; ii<maxN; ii++) {
    bits->setBinary(width[ii], random[ii]);

    position += width[ii];

    assert(bits->getPosition() == position);
  }

  if (verbose) {
    fprintf(stderr, "Testing  %lu numbers with total length %lu bits encoded into %lu bits.\n", maxN, Nbits, bits->getLength());
    fprintf(stderr, "\n");
    fprintf(stderr, "                     value  b  data word\n");
    fprintf(stderr, "--  ---------------------- --  ----------------------------------------------------------------\n");
    for (uint32 ii=0; ii<10; ii++)
      fprintf(stderr, "%02u  %22lu %2u  %s\n", ii, random[ii], width[ii], bits->displayWord(ii));
    fprintf(stderr, "--  ---------------------- --  ----------------------------------------------------------------\n");
    fprintf(stderr, "\n");
  }

  bits->setPosition(0);

  for (uint64 position=0, ii=0; ii<maxN; ii++) {
    uint64  b = bits->getBinary(width[ii]);

    position += width[ii];

    assert(random[ii] == b);
    assert(bits->getPosition() == position);
  }

  {
    char          N[FILENAME_MAX+1];

    snprintf(N, FILENAME_MAX, "bitsTest-binary-%02u.sb", maxWidth);

    if (verbose)
      fprintf(stderr, "Writing  %lu numbers with total length %lu bits encoded into %lu bits.\n", maxN, Nbits, bits->getLength());    

    writeBuffer *Bw = new writeBuffer(N, "w");
    bits->dumpToBuffer(Bw);
    delete Bw;

    if (verbose)
      fprintf(stderr, "Forgetting stuffedBits.\n");

    delete bits;
    bits = nullptr;

    if (verbose)
      fprintf(stderr, "Reading  %lu numbers with total length %lu bits.\n", maxN, Nbits);

    readBuffer *Br = new readBuffer(N);
    bits = new stuffedBits;
    assert(bits->loadFromBuffer(Br) == true);
    delete Br;
  }

  if (verbose)
    fprintf(stderr, "Testing  %lu numbers with total length %lu bits encoded into %lu bits.\n", maxN, Nbits, bits->getLength());

  for (uint64 position=0, ii=0; ii<maxN; ii++) {
    uint64  b = bits->getBinary(width[ii]);

    position += width[ii];

    assert(random[ii] == b);
    assert(bits->getPosition() == position);
  }

  delete    bits;
  delete [] random;
  delete [] width;
}









void
testPrefixFree(bool verbose, uint64 length, uint32 type) {
  uint64      maxN   = length;
  uint64      Nbits  = 0;
  uint32     *width  = new uint32 [maxN];
  uint64     *random = new uint64 [maxN];
  uint64      histo[65];
  uint64      masks[65];
  mtRandom    mt;

  if ((verbose) && (type == 0))
    fprintf(stderr, "Testing stuffedBits EliasGamma encoding with %lu numbers.\n", length);
  if ((verbose) && (type == 1))
    fprintf(stderr, "Testing stuffedBits EliasDelta encoding with %lu numbers.\n", length);
  if ((verbose) && (type == 2))
    fprintf(stderr, "Testing stuffedBits Zeckendorf encoding with %lu numbers.\n", length);

  for (uint32 ii=0; ii<65; ii++) {
    histo[ii] = 0;
    masks[ii] = buildLowBitMask<uint64>(ii);
  }

  for (uint64 ii=0; ii<maxN; ii++) {
    width[ii]  = mt.mtRandom32() % 64 + 1;
    random[ii] = mt.mtRandom64() & masks[width[ii]];

    Nbits += width[ii];
    histo[width[ii]]++;

    if (random[ii] == 0)
      ii--;
  }

  if (verbose)
    fprintf(stderr, "Setting  %lu numbers with total length %lu bits.\n", maxN, Nbits);

  stuffedBits *bits = new stuffedBits;

  for (uint64 ii=0; ii<maxN; ii++) {
    switch (type) {
      case 0:
        bits->setEliasGamma(random[ii]);
        break;
      case 1:
        bits->setEliasDelta(random[ii]);
        break;
      case 2:
        bits->setZeckendorf(random[ii]);
        break;
    }
  }

  if (verbose) {
    fprintf(stderr, "Testing  %lu numbers with total length %lu bits encoded into %lu bits.\n", maxN, Nbits, bits->getLength());
    fprintf(stderr, "\n");
    fprintf(stderr, "                     value  b  data word\n");
    fprintf(stderr, "--  ---------------------- --  ----------------------------------------------------------------\n");
    for (uint32 ii=0; ii<10; ii++)
      fprintf(stderr, "%02u  %22lu %2u  %s\n", ii, random[ii], width[ii], bits->displayWord(ii));
    fprintf(stderr, "--  ---------------------- --  ----------------------------------------------------------------\n");
    fprintf(stderr, "\n");
  }

  bits->setPosition(0);

  for (uint64 ii=0; ii<maxN; ii++) {
    uint64  b = 0;

    switch (type) {
      case 0:
        b = bits->getEliasGamma();
        break;
      case 1:
        b = bits->getEliasDelta();
        break;
      case 2:
        b = bits->getZeckendorf();
        break;
    }

    if (b != random[ii])
      fprintf(stderr, "Failed at ii %lu expect random=%lu got b=%lu\n",
              ii, random[ii], b);
    assert(random[ii] == b);
  }

  delete    bits;
  delete [] random;
  delete [] width;
}



void
testIO(bool verbose, uint64 length) {
}




int
main(int argc, char **argv) {
  bool   verbose = false;
  uint64 length  = 10 * 1000000;
  uint32 width   = 21;
  int32  arg     = 1;
  int32  err     = 0;

  bool  tBitArray       = false;
  bool  tWordArray      = false;
  bool  tWordArraySpeed = false;
  bool  tUnary          = false;
  bool  tBinary         = false;
  bool  tEliasGamma     = false;
  bool  tEliasDelta     = false;
  bool  tZeckendorf     = false;

  omp_set_num_threads(1);

  testMasks(false);
  testLogBaseTwo(false);
  testSaveClear(false);
  testExpandCompress(false, 21);

  while (arg < argc) {
    if      (strcmp(argv[arg], "-verbose") == 0) {
      verbose = true;
    }
    else if (strcmp(argv[arg], "-length") == 0) {
      length = strtouint64(argv[++arg]) * 1000000;
    }
    else if (strcmp(argv[arg], "-width") == 0) {
      width = strtouint32(argv[++arg]);
    }
    else if (strcmp(argv[arg], "-threads") == 0) {
      omp_set_num_threads(strtouint32(argv[++arg]));
    }

    else if (strcmp(argv[arg], "-masks") == 0) {
      testMasks(verbose);
    }
    else if (strcmp(argv[arg], "-logbasetwo") == 0) {
      testLogBaseTwo(verbose);
    }
    else if (strcmp(argv[arg], "-clear") == 0) {
      testSaveClear(verbose);
    }
    else if (strcmp(argv[arg], "-expand") == 0) {
      testExpandCompress(verbose, 21);
    }
    else if (strcmp(argv[arg], "-expandfail") == 0) {
      testExpandCompress(verbose, 22);
    }
    else if (strcmp(argv[arg], "-fibonacci") == 0) {
      testFibonacciNumbers(verbose);
    }

    else if (strcmp(argv[arg], "-all") == 0) {
      tBitArray   = true;
      tWordArray  = true;
      tUnary      = true;
      tBinary     = true;
      tEliasGamma = true;
      tEliasDelta = true;
      tZeckendorf = true;
    }

    else if (strcmp(argv[arg], "-bitarray") == 0) {
      tBitArray = true;
    }

    else if (strcmp(argv[arg], "-wordarray") == 0) {
      tWordArray = true;
    }
    else if (strcmp(argv[arg], "-wordarrayspeed") == 0) {
      tWordArraySpeed = true;
    }

    else if (strcmp(argv[arg], "-unary") == 0) {
      tUnary = true;
    }
    else if (strcmp(argv[arg], "-binary") == 0) {
      tBinary = true;
    }
    else if (strcmp(argv[arg], "-eliasgamma") == 0) {
      tEliasGamma = true;
    }
    else if (strcmp(argv[arg], "-eliasdelta") == 0) {
      tEliasDelta = true;
    }
    else if (strcmp(argv[arg], "-zeckendorf") == 0) {
      tZeckendorf = true;
    }

    else {
      err++;
    }

    arg++;
  }

  if (argc == 1)
    err++;

  if (err) {
    fprintf(stderr, "OPTIONS\n");
    fprintf(stderr, "  -verbose           Report what the tests are doing.\n");
    fprintf(stderr, "  -length            Change the length of the test, in millions.  Default 10.\n");
    fprintf(stderr, "  -threads           Use multiple threads, for -unary and -binary.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "BASIC TESTS\n");
    fprintf(stderr, "  These run automatically, every time (except -expandfail).  Running\n");
    fprintf(stderr, "  explicitly with -verbose gives more details.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -masks             buildLowBitMask() and buildHighBitMask()\n");
    fprintf(stderr, "  -logbasetwo        countNumberOfBits64()\n");
    fprintf(stderr, "  -clear             [save|clear][Left|Middle|Right]Bits()\n");
    fprintf(stderr, "  -expand            expandTo3() and compressTo2()\n");
    fprintf(stderr, "  -expandfail        expandTo3() and compressTo2(), success if assert() fails!\n");
    fprintf(stderr, "  -fibonacci         fibonacciNumber()\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "ENCODING TESTS\n");
    fprintf(stderr, "  Test bitArray, wordArray and stuffedBits.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -all               Run all of the below, except -wordarrayspeed.\n");
    fprintf(stderr, "                       (about five minutes with no threads)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -bitarray          bitArray::setBit(), bitArray::getBit() and bitArray::flipBit()\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -wordarray         wordArray::set() and wordArray::get()\n");
    fprintf(stderr, "  -wordarrayspeed    wordArray speed against plain arrays\n");
    fprintf(stderr, "    -bits N          set size of word in speed test\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -unary             stuffedBits::setUnary() for values up to 8193\n");
    fprintf(stderr, "  -binary            stuffedBits::setBinary() for all widths up to 64\n");
    fprintf(stderr, "                     stuffedBits::dumpToFile() and stuffedBits::loadFromFile()\n");
    fprintf(stderr, "                       (by far the slowest, benefits from -threads)\n");
    fprintf(stderr, "                       (also tests input/output)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -eliasgamma        stuffedBits::setEliasGamma()\n");
    fprintf(stderr, "  -eliasdelta        stuffedBits::setEliasDelta()\n");
    fprintf(stderr, "  -zeckendorf        stuffedBits::setZeckendorf()\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "The basic tests are always run, silently, regardless of options.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Program crashes if any test fails (except for -expandfail, which is supposed to fail).\n");
    return(0);
  }


  if (tBitArray) {
    testBitArray(verbose,     0, 10000);
    testBitArray(verbose,     1, 10000);
    testBitArray(verbose,   100, 10000);
    testBitArray(verbose,  1000, 10000);
    testBitArray(verbose, 10000, 10000);
    testBitArray(verbose, 40000, 10000);
  }

  if (tWordArray) {
    testWordArray(verbose,    1000,   1,       256);
    testWordArray(verbose,   10000,   2,       384);
    testWordArray(verbose,   10000,   3,       512);
    testWordArray(verbose,   10000,   5,       768);
    testWordArray(verbose,  100000,   7, 1 * 32768);   //  4 KB
    testWordArray(verbose,  100000,   9, 4 * 32768);
    testWordArray(verbose,  100000,  13, 4 * 32768);
    testWordArray(verbose,  200000,  17, 4 * 32768);
    testWordArray(verbose,  200000,  19, 4 * 32768);
    testWordArray(verbose,  200000,  23, 4 * 32768);
    testWordArray(verbose,  300000,  31, 4 * 32768);
    testWordArray(verbose,  300000,  37, 4 * 32768);
    testWordArray(verbose,  300000,  41, 4 * 32768);
    testWordArray(verbose,  400000,  43, 4 * 32768);
    testWordArray(verbose,  400000,  47, 4 * 32768);
    testWordArray(verbose,  400000,  53, 4 * 32768);
    testWordArray(verbose,  500000,  59, 4 * 32768);
    testWordArray(verbose,  500000,  61, 4 * 32768);
    testWordArray(verbose,  500000,  67, 4 * 32768);
    testWordArray(verbose,  600000,  71, 4 * 32768);
    testWordArray(verbose,  600000,  73, 4 * 32768);
    testWordArray(verbose,  600000,  79, 4 * 32768);
    testWordArray(verbose,  700000,  83, 4 * 32768);
    testWordArray(verbose,  700000,  89, 6 * 32768);
    testWordArray(verbose,  700000,  97, 6 * 32768);
    testWordArray(verbose,  700000, 101, 6 * 32768);
    testWordArray(verbose,  800000, 103, 6 * 32768);
    testWordArray(verbose,  800000, 107, 8 * 32768);
    testWordArray(verbose,  800000, 109, 8 * 32768);
    testWordArray(verbose,  800000, 113, 8 * 32768);
    testWordArray(verbose,  900000, 127, 8 * 32768);

    testWordArray(verbose,  500000,  32, 4 * 32768);
    testWordArray(verbose,  500000,  64, 6 * 32768);
    testWordArray(verbose,  500000, 128, 8 * 32768);
  }

  if (tWordArraySpeed) {
    testWordArraySpeed(verbose, length, width);
  }

  if (tUnary) {
    uint32  sizes[14] = { 1, 2, 3, 4, 62, 63, 64, 65, 126, 127, 128, 129, 257, 8193 };

    //#pragma omp parallel for
    for (uint32 ss=0; ss < 14; ss++)
      testUnary(verbose, length, sizes[ss]);
  }

  if (tBinary) {
    //#pragma omp parallel for
    for (uint32 xx=1; xx<=64; xx++)
      testBinary(verbose, length, xx);
  }

  if (tEliasGamma) {
    testPrefixFree(verbose, length, 0);
  }

  if (tEliasDelta) {
    testPrefixFree(verbose, length, 1);
  }

  if (tZeckendorf) {
    testPrefixFree(verbose, length, 2);
  }

  return(0);
}
