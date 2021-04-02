
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
