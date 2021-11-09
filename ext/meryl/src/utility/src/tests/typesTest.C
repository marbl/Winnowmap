
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

char const *minu128 = "0";
char const *maxu128 = "340282366920938463463374607431768211455";

char const *minu64 = "0";
char const *maxu64 = "18446744073709551615";

char const *minu32 = "0";
char const *maxu32 = "4294967295";

char const *minu16 = "0";
char const *maxu16 = "65535";

char const *minu8 = "0";
char const *maxu8 = "255";


char const *min128 = "-170141183460469231731687303715884105728";
char const *max128 = "+170141183460469231731687303715884105727";

char const *min64 = "-9223372036854775808";
char const *max64 = "+9223372036854775807";

char const *min32 = "-2147483648";
char const *max32 =  "2147483647";

char const *min16 = "-32768";
char const *max16 =  "32767";

char const *min8 = "-128";
char const *max8 =  "127";




bool
test_strto(void) {

  fprintf(stderr, "Testing conversion of string to unsigned integers.\n");

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

  fprintf(stderr, "Testing conversion of string to signed integers.\n");

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

  fprintf(stderr, "Testing conversion of string to signed integers.\n");

  char  outstr[13], *outp;

  for (uint32 tt=0; tt<13; tt++)
    outstr[tt] = 100;

  for (uint32 tt=0; 1; tt++) {
    outp = toDec(tt, outstr);

    if ((tt < 25) || (tt > uint32max-25) || ((tt % 1844751) == 0))
      fprintf(stderr, "tt %10u len %d -- out '%s'\n", tt, (int32)(outp - outstr), outstr);

    if      (tt < 10)            assert(outstr+1  == outp);
    else if (tt < 100)           assert(outstr+2  == outp);
    else if (tt < 1000)          assert(outstr+3  == outp);
    else if (tt < 10000)         assert(outstr+4  == outp);
    else if (tt < 100000)        assert(outstr+5  == outp);
    else if (tt < 1000000)       assert(outstr+6  == outp);
    else if (tt < 10000000)      assert(outstr+7  == outp);
    else if (tt < 100000000)     assert(outstr+8  == outp);
    else if (tt < 1000000000)    assert(outstr+9  == outp);

    assert(outp[0] == 0);
    assert(outp[1] == 100);
    assert(outp[2] == 100);

    assert(strtouint32(outstr) == tt);

    if (tt == uint32max)   //  Otherwise we look infinitely.  (Actually, we fail the
      break;               //  '==100' assert above when we wrap around.)
  }

  fprintf(stderr, "Tests passed.\n");

  return(true);
}







int
main(int argc, char **argv) {
  int32 arg=1;
  int32 err=0;

  omp_set_num_threads(1);

  while (arg < argc) {
    if      (strcmp(argv[arg], "-h") == 0) {
      err++;
    }

    else if (strcmp(argv[arg], "-something") == 0) {
      //testSomething();
    }

    else {
      err++;
    }

    arg++;
  }

  if (err)
    fprintf(stderr, "ERROR: didn't parse command line.\n"), exit(1);


  test_strto();

  exit(0);
}
