
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

#include "types.H"
#include "strings.H"

//
//  Test strings for test_strto()
//

char const *minu128 =  "0";
char const *maxu128 =  "340282366920938463463374607431768211455";
char const *min128  = "-170141183460469231731687303715884105728";
char const *max128  = "+170141183460469231731687303715884105727";

char const *minu64 =  "0";
char const *maxu64 =  "18446744073709551615";
char const *min64  = "-9223372036854775808";
char const *max64  = "+9223372036854775807";

char const *minu32 =  "0";
char const *maxu32 =  "4294967295";
char const *min32  = "-2147483648";
char const *max32  =  "2147483647";

char const *minu16 =  "0";
char const *maxu16 =  "65535";
char const *min16  = "-32768";
char const *max16  =  "32767";

char const *minu8 =  "0";
char const *maxu8 =  "255";
char const *min8  = "-128";
char const *max8  =  "127";

//
//  Test ranges for test_decodeRange()
//

using merylutil::ansiCode;

template<typename R>
void
test_asciiXXXtoInt(R (*XXXtoInt)(char d), bool (*isXXXdigit)(char d), char const *name, char const *type=nullptr) {
  fprintf(stdout, "|--Full %s table.\n", name);
  if (type)
    fprintf(stdout, "|  (undefined for non-%s inputs)\n", type);
  fprintf(stdout, "|\n");
  fprintf(stdout, "|     ");
  for (uint32 jj=0; jj<16; jj++) {
    fprintf(stdout, "    .%x", jj);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "|     ");
  for (uint32 jj=0; jj<16; jj++) {
    fprintf(stdout, "    --");
  }
  fprintf(stdout, "\n");

  const char *inv = makeAnsiEscapeSequence( (const ansiCode[]) { ansiCode::Green, ansiCode::Bold, ansiCode::END } );
  const char *nml = makeAnsiEscapeSequence( (const ansiCode[]) { ansiCode::Normal,                ansiCode::END } );

  for (uint32 ii=0; ii<16; ii++) {
    fprintf(stdout, "|  %x. |", ii);
    for (uint32 jj=0; jj<16; jj++) {
      uint8 ch = jj << 4 | ii;
      if ((*isXXXdigit)(ch))
        fprintf(stdout, " %s%c=%02xh%s", inv, integerToLetter(ch), (*XXXtoInt)(ch), nml);
      else
        fprintf(stdout, " %c=%02xh", integerToLetter(ch), (*XXXtoInt)(ch));
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "|--\n");
  fprintf(stdout, "\n");
}


void
test_asciiBinToInt(bool display=false) {
  if (display)
    test_asciiXXXtoInt(asciiBinToInteger, isBinDigit, "asciiBinToInteger()", "binary");

  assert(asciiBinToInteger('0') == 0x00);
  assert(asciiBinToInteger('1') == 0x01);
}

void
test_asciiOctToInt(bool display=false) {
  if (display)
    test_asciiXXXtoInt(asciiOctToInteger, isOctDigit, "asciiOctToInteger()", "octal");

  assert(asciiOctToInteger('0') == 0x00);
  assert(asciiOctToInteger('1') == 0x01);
  assert(asciiOctToInteger('2') == 0x02);
  assert(asciiOctToInteger('3') == 0x03);
  assert(asciiOctToInteger('4') == 0x04);
  assert(asciiOctToInteger('5') == 0x05);
  assert(asciiOctToInteger('6') == 0x06);
  assert(asciiOctToInteger('7') == 0x07);
}

void
test_asciiDecToInt(bool display=false) {
  if (display)
    test_asciiXXXtoInt(asciiDecToInteger, isDecDigit, "asciiDecToInteger()", "decimal");

  assert(asciiDecToInteger('0') == 0x00);
  assert(asciiDecToInteger('1') == 0x01);
  assert(asciiDecToInteger('2') == 0x02);
  assert(asciiDecToInteger('3') == 0x03);
  assert(asciiDecToInteger('4') == 0x04);
  assert(asciiDecToInteger('5') == 0x05);
  assert(asciiDecToInteger('6') == 0x06);
  assert(asciiDecToInteger('7') == 0x07);
  assert(asciiDecToInteger('8') == 0x08);
  assert(asciiDecToInteger('9') == 0x09);
}

void
test_asciiHexToInt(bool display=false) {

  if (display)
    test_asciiXXXtoInt(asciiHexToInteger, isHexDigit, "asciiHexToInteger()", "hexadecimal");

  assert(asciiHexToInteger('0') == 0x00);
  assert(asciiHexToInteger('1') == 0x01);
  assert(asciiHexToInteger('2') == 0x02);
  assert(asciiHexToInteger('3') == 0x03);
  assert(asciiHexToInteger('4') == 0x04);
  assert(asciiHexToInteger('5') == 0x05);
  assert(asciiHexToInteger('6') == 0x06);
  assert(asciiHexToInteger('7') == 0x07);
  assert(asciiHexToInteger('8') == 0x08);
  assert(asciiHexToInteger('9') == 0x09);
  assert(asciiHexToInteger('A') == 0x0A);
  assert(asciiHexToInteger('B') == 0x0B);
  assert(asciiHexToInteger('C') == 0x0C);
  assert(asciiHexToInteger('D') == 0x0D);
  assert(asciiHexToInteger('E') == 0x0E);
  assert(asciiHexToInteger('F') == 0x0F);
  assert(asciiHexToInteger('a') == 0x0a);
  assert(asciiHexToInteger('b') == 0x0b);
  assert(asciiHexToInteger('c') == 0x0c);
  assert(asciiHexToInteger('d') == 0x0d);
  assert(asciiHexToInteger('e') == 0x0e);
  assert(asciiHexToInteger('f') == 0x0f);
}


bool
test_strto(bool verbose) {

  if (verbose)
    fprintf(stdout, "Testing conversion of string to unsigned integers.\n");

  assert(strtouint128(minu128) == uint128min);
  assert(strtouint128(maxu128) == uint128max);

  assert(strtouint64(minu64) == uint64min);
  assert(strtouint64(maxu64) == uint64max);

  assert(strtouint32(minu32) == uint32min);
  assert(strtouint32(maxu32) == uint32max);

  assert(strtouint16(minu16) == uint16min);
  assert(strtouint16(maxu16) == uint16max);

  assert(strtouint8(minu8) == uint8min);
  assert(strtouint8(maxu8) == uint8max);

  if (verbose)
    fprintf(stdout, "Testing conversion of string to signed integers.\n");

  assert(strtoint128(min128) == int128min);
  assert(strtoint128(max128) == int128max);

  assert(strtoint64(min64) == int64min);
  assert(strtoint64(max64) == int64max);

  assert(strtoint32(min32) == int32min);
  assert(strtoint32(max32) == int32max);

  assert(strtoint16(min16) == int16min);
  assert(strtoint16(max16) == int16max);

  assert(strtoint8(min8) == int8min);
  assert(strtoint8(max8) == int8max);

  if (verbose)
    fprintf(stdout, "Testing conversion of string to signed integers.\n");

  char  outstr[32], *outp;

  for (uint32 tt=0; tt<32; tt++)   //  Initialize string to junk.
    outstr[tt] = 100;              //

  for (uint64 tt=0; tt<999999999999999; ) {
    outp = toDec(tt, outstr);

    if (verbose)
      fprintf(stdout, "  %2u digits - %s\n", (int32)(outp - outstr), outstr);

    assert(outp[0] == 0);     //  The output pointer should point to a NUL
    assert(outp[1] == 100);   //  byte followed by uninitialized junk.
    assert(outp[2] == 100);   //

    assert(strtouint64(outstr) == tt);   //  The outstr should decode back to the input integer.

    if      (tt < 10)              assert(outstr+1   == outp);   //  Check outp is at the expected
    else if (tt < 100)             assert(outstr+2   == outp);   //  location for the number of
    else if (tt < 1000)            assert(outstr+3   == outp);   //  digits in the input.
    else if (tt < 10000)           assert(outstr+4   == outp);
    else if (tt < 100000)          assert(outstr+5   == outp);
    else if (tt < 1000000)         assert(outstr+6   == outp);
    else if (tt < 10000000)        assert(outstr+7   == outp);
    else if (tt < 100000000)       assert(outstr+8   == outp);
    else if (tt < 1000000000)      assert(outstr+9   == outp);
    else if (tt < 10000000000)     assert(outstr+10  == outp);
    else if (tt < 100000000000)    assert(outstr+11  == outp);
    else if (tt < 1000000000000)   assert(outstr+12  == outp);
    else if (tt < 10000000000000)  assert(outstr+13  == outp);
    else if (tt < 100000000000000) assert(outstr+14  == outp);

    if      (tt < 10)            tt += 1;
    else if (tt < 100)           tt += 1;
    else if (tt < 1000)          tt += 3;
    else if (tt < 10000)         tt += 7;
    else if (tt < 100000)        tt += 73;
    else if (tt < 1000000)       tt += 337;
    else if (tt < 10000000)      tt += 773;
    else if (tt < 100000000)     tt += 3733;
    else if (tt < 1000000000)    tt += 33377;
    else if (tt < 10000000000)   tt += 333337;
    else if (tt < 100000000000)  tt += 33373777;
    else if (tt < 1000000000000) tt += 333377777;
    else                         tt += 3377737777;
  }

  if (verbose)
    fprintf(stdout, "\nTests passed.\n");

  return(true);
}



void
test_isASCII(void) {
  assert(isNUL  ('\0') == true);
  assert(isBEL  ('\a') == true);
  assert(isBS   ('\b') == true);
  assert(isTab  ('\t') == true);
  assert(isLF   ('\n') == true);
  assert(isVT   ('\v') == true);
  assert(isFF   ('\f') == true);
  assert(isCR   ('\r') == true);
  assert(isSpace( ' ') == true);

  assert(isEndOfLine('\n') == true);
  assert(isEndOfLine('\r') == true);

  assert(isWhiteSpace( ' ') == true);
  assert(isWhiteSpace('\t') == true);
  assert(isWhiteSpace('\n') == true);
  assert(isWhiteSpace('\r') == true);
  assert(isWhiteSpace('\f') == true);
  assert(isWhiteSpace('\v') == true);
}



//template<typename integerType>
typedef uint32 integerType;
void
test_decodeRange(const char *range=nullptr) {
  char          str[1024];
  integerType   minv;
  integerType   maxv;

  std::vector<integerType>  mins;
  std::vector<integerType>  maxs;

  if (range == nullptr) {
    decodeRange("0-0",   minv, maxv);  assert(minv == 0);   assert(maxv == 0);
    decodeRange("0-1",   minv, maxv);  assert(minv == 0);   assert(maxv == 1);
    decodeRange("1-0",   minv, maxv);  assert(minv == 1);   assert(maxv == 0);
    decodeRange("1-1",   minv, maxv);  assert(minv == 1);   assert(maxv == 1);

    decodeRange("00-01", minv, maxv);  assert(minv == 0);   assert(maxv == 1);
    decodeRange("01-00", minv, maxv);  assert(minv == 1);   assert(maxv == 0);
    decodeRange("1-100", minv, maxv);  assert(minv == 1);   assert(maxv == 100);
    decodeRange("100-1", minv, maxv);  assert(minv == 100); assert(maxv == 1);

    decodeRange("0-0,1-1,22-10", mins, maxs);
    assert(mins[0] == 0);    assert(maxs[0] == 0);
    assert(mins[1] == 1);    assert(maxs[1] == 1);
    assert(mins[2] == 22);   assert(maxs[2] == 10);
  }
  else {
    decodeRange(range, mins, maxs);

    for (uint64 ii=0; ii<mins.size(); ii++)
      fprintf(stdout, "%u %u\n", mins[ii], maxs[ii]);
  }
}



int
main(int argc, char **argv) {
  int32 arg=1;
  int32 err=0;

  omp_set_num_threads(1);

  test_asciiBinToInt(false);
  test_asciiOctToInt(false);
  test_asciiDecToInt(false);
  test_asciiHexToInt(false);
  test_strto(false);
  test_isASCII();

  while (arg < argc) {
    if      (strcmp(argv[arg], "-c")  == 0)   test_strto(true);
    else if (strcmp(argv[arg], "-r")  == 0)   test_decodeRange();
    else if (strcmp(argv[arg], "-R")  == 0)   test_decodeRange(argv[++arg]);

    else if (strcmp(argv[arg], "-ab") == 0)   test_asciiBinToInt(true);
    else if (strcmp(argv[arg], "-ao") == 0)   test_asciiOctToInt(true);
    else if (strcmp(argv[arg], "-ad") == 0)   test_asciiDecToInt(true);
    else if (strcmp(argv[arg], "-ah") == 0)   test_asciiHexToInt(true);

    else if (strcmp(argv[arg], "-ieol")     == 0)   test_asciiXXXtoInt(isEndOfLine,  isEndOfLine,  "isEndOfLine()");
    else if (strcmp(argv[arg], "-iwhite")   == 0)   test_asciiXXXtoInt(isWhiteSpace, isWhiteSpace, "isWhiteSpace()");
    else if (strcmp(argv[arg], "-ivisible") == 0)   test_asciiXXXtoInt(isVisible,    isVisible,    "isVisible()");
    else if (strcmp(argv[arg], "-iletter")  == 0)   test_asciiXXXtoInt(isLetter,     isLetter,     "isLetter()");

    else {
      fprintf(stderr, "usage: %s [test]\n", argv[0]);
      fprintf(stderr, "  -c         run test converting strings to integers.\n");
      fprintf(stderr, "  -r         run test decoding ranges.\n");
      fprintf(stderr, "  -R R       decode range R and print values.\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -ab        show table of conversion of asciiBinToInteger().\n");
      fprintf(stderr, "  -ao        show table of conversion of asciiOctToInteger().\n");
      fprintf(stderr, "  -ad        show table of conversion of asciiDecToInteger().\n");
      fprintf(stderr, "  -ah        show table of conversion of asciiHexToInteger().\n");
      fprintf(stderr, "\n");
      fprintf(stderr, "  -ieol      show table of isEndOfLine().\n");
      fprintf(stderr, "  -iwhite    show table of isEndOfLine().\n");
      fprintf(stderr, "  -ivisible  show table of isEndOfLine().\n");
      fprintf(stderr, "  -iletter   show table of isEndOfLine().\n");
      fprintf(stderr, "\n");
      return 1;
    }

    arg++;
  }

  return 0;
}
