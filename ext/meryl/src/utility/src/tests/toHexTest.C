
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
#include "mt19937ar.H"



uint128
decodeHex(char const *str) {
  uint128   ret = 0;

  for (uint32 ii=0; str[ii]; ii++) {
    ret <<= 4;

    if (str[ii] <= '9')
      ret |= (str[ii] - '0');
    else
      ret |= (str[ii] - 'a' + 10);
  }

  return(ret);
}


uint128
decodeDec(char const *str) {
  uint128   ret = 0;

  for (uint32 ii=0; str[ii]; ii++) {
    ret *= 10;
    ret += (str[ii] - '0');
  }

  return(ret);
}


uint128
decodeOct(char const *str) {
  uint128   ret = 0;

  for (uint32 ii=0; str[ii]; ii++) {
    ret <<= 3;
    ret |= (str[ii] - '0');
  }

  return(ret);
}


uint128
decodeBin(char const *str) {
  uint128   ret = 0;

  for (uint32 ii=0; str[ii]; ii++) {
    ret <<= 1;
    ret |= (str[ii] - '0');
  }

  return(ret);
}




int
main(int argc, char **argv) {
  uint128  u1 = 0;   char const *s1 = nullptr;
  uint64   u2 = 0;   char const *s2 = nullptr;
  uint32   u3 = 0;   char const *s3 = nullptr;
  uint16   u4 = 0;   char const *s4 = nullptr;

  uint32   limit = (argc > 1) ? strtouint32(argv[1]) : 1;
  uint32   err   = (argc > 1) ? 0                    : 1;

  mtRandom mt;

  for (uint32 ii=0; ii<limit; ii++) {
    uint32 sh1 = (mt.mtRandom32() % 128);
    uint32 sh2 = (mt.mtRandom32() %  64);
    uint32 sh3 = (mt.mtRandom32() %  32);
    uint32 sh4 = (mt.mtRandom32() %  16);

    u1 = mt.mtRandom64();  u1 <<= 64;   u1 |= mt.mtRandom64();
    u2 = mt.mtRandom64();
    u3 = mt.mtRandom32();
    u4 = mt.mtRandom32();

    if (mt.mtRandom32() & 0x01) {
      uint128 U1 = u1 >> sh1;   //  Isn't for efficiency, but to get around
      uint64  U2 = u2 >> sh2;   //  'u4 >> sh4' being promoted to an 'int' and
      uint32  U3 = u3 >> sh3;   //  then not finding a version of toHex() that
      uint16  U4 = u4 >> sh4;   //  works.  Only happens on u4.

      if (err) {
        fprintf(stderr, "128: %32s >>%2u %s\n", toHex(u1), sh1, toHex(U1));
        fprintf(stderr, " 64: %32s >>%2u %s\n", toHex(u2), sh2, toHex(U2));
        fprintf(stderr, " 32: %32s >>%2u %s\n", toHex(u3), sh3, toHex(U3));
        fprintf(stderr, " 16: %32s >>%2u %s\n", toHex(u4), sh4, toHex(U4));
        fprintf(stderr, "\n");
      }

      u1 = U1;
      u2 = U2;
      u3 = U3;
      u4 = U4;
    }
    else {
      uint128 U1 = u1 << sh1;
      uint64  U2 = u2 << sh2;
      uint32  U3 = u3 << sh3;
      uint16  U4 = u4 << sh4;

      if (err) {
        fprintf(stderr, "128: %32s <<%2u %s\n", toHex(u1), sh1, toHex(U1));
        fprintf(stderr, " 64: %32s <<%2u %s\n", toHex(u2), sh2, toHex(U2));
        fprintf(stderr, " 32: %32s <<%2u %s\n", toHex(u3), sh3, toHex(U3));
        fprintf(stderr, " 16: %32s <<%2u %s\n", toHex(u4), sh4, toHex(U4));
        fprintf(stderr, "\n");
      }

      u1 = U1;
      u2 = U2;
      u3 = U3;
      u4 = U4;
    }

    s1 = toHex(u1);
    s2 = toDec(u1);
    s3 = toOct(u1);
    s4 = toBin(u1);

    if (u1 != decodeHex(s1))   err += fprintf(stderr, "ERROR in HEX\n");
    if (u1 != decodeDec(s2))   err += fprintf(stderr, "ERROR in DEC\n");
    if (u1 != decodeOct(s3))   err += fprintf(stderr, "ERROR in OCT\n");
    if (u1 != decodeBin(s4))   err += fprintf(stderr, "ERROR in BIN\n");

    if (err) {
      fprintf(stderr, "128\n");
      fprintf(stderr, "  HEX %s\n", s1);
      fprintf(stderr, "  DEC %s\n", s2);
      fprintf(stderr, "  OCT %s\n", s3);
      fprintf(stderr, "  BIN %s\n", s4);
    }

    s1 = toHex(u2);
    s2 = toDec(u2);
    s3 = toOct(u2);
    s4 = toBin(u2);

    if (u2 != strtoull(s1, nullptr, 16))   fprintf(stderr, "ERROR in HEX\n");
    if (u2 != strtoull(s2, nullptr, 10))   fprintf(stderr, "ERROR in DEC\n");
    if (u2 != strtoull(s3, nullptr,  8))   fprintf(stderr, "ERROR in OCT\n");
    if (u2 != strtoull(s4, nullptr,  2))   fprintf(stderr, "ERROR in BIN\n");

    if (err) {
      fprintf(stderr, "64\n");
      fprintf(stderr, "  HEX %s\n", s1);
      fprintf(stderr, "  DEC %s\n", s2);
      fprintf(stderr, "  OCT %s\n", s3);
      fprintf(stderr, "  BIN %s\n", s4);
    }

    s1 = toHex(u3);
    s2 = toDec(u3);
    s3 = toOct(u3);
    s4 = toBin(u3);

    if (u3 != strtoull(s1, nullptr, 16))   fprintf(stderr, "ERROR in HEX\n");
    if (u3 != strtoull(s2, nullptr, 10))   fprintf(stderr, "ERROR in DEC\n");
    if (u3 != strtoull(s3, nullptr,  8))   fprintf(stderr, "ERROR in OCT\n");
    if (u3 != strtoull(s4, nullptr,  2))   fprintf(stderr, "ERROR in BIN\n");

    if (err) {
      fprintf(stderr, "32\n");
      fprintf(stderr, "  HEX %s\n", s1);
      fprintf(stderr, "  DEC %s\n", s2);
      fprintf(stderr, "  OCT %s\n", s3);
      fprintf(stderr, "  BIN %s\n", s4);
    }

    s1 = toHex(u4);
    s2 = toDec(u4);
    s3 = toOct(u4);
    s4 = toBin(u4);

    if (u4 != strtoull(s1, nullptr, 16))   fprintf(stderr, "ERROR in HEX\n");
    if (u4 != strtoull(s2, nullptr, 10))   fprintf(stderr, "ERROR in DEC\n");
    if (u4 != strtoull(s3, nullptr,  8))   fprintf(stderr, "ERROR in OCT\n");
    if (u4 != strtoull(s4, nullptr,  2))   fprintf(stderr, "ERROR in BIN\n");

    if (err) {
      fprintf(stderr, "16\n");
      fprintf(stderr, "  HEX %s\n", s1);
      fprintf(stderr, "  DEC %s\n", s2);
      fprintf(stderr, "  OCT %s\n", s3);
      fprintf(stderr, "  BIN %s\n", s4);
    }
  }

  return(0);
}
