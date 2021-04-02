
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

int
main(int argc, char **argv) {

  if (argc != 2) {
    fprintf(stderr, "usage: %s text-file\n", argv[0]);
    fprintf(stderr, "  Converts the input text-file into 64-bit and 32-bit integer\n");
    fprintf(stderr, "  constants, for use as magic numbers in data files.  If you then\n");
    fprintf(stderr, "  write this integer constant to a file, it'll appear as readable\n");
    fprintf(stderr, "  text in the file.  The input file is limited to 4 KB.\n");
    return(1);
  }

  uint32  ccLen = 0;
  uint32  ccMax = 4096;
  char   *cc    = new char [ccMax];
  FILE   *F;

  memset(cc, 0, ccMax);

  F = fopen(argv[1], "r");
  ccLen = fread(cc, sizeof(char), 4096, F);
  fclose(F);

  F = fopen(argv[1], "r");
  for (uint32 ii=0, nn=0; ii<ccLen; ii += 8, nn++) {
    uint64 u64;

    fread(&u64, sizeof(uint64), 1, F);       //  You can get away with only the char array,
    assert(u64 == *((uint64 *)(cc + ii)));   //  but I'm not sure what will happen on big-endian.

    fprintf(stdout, "uint64 u64_%02u = 0x%016lxllu;  //  %c%c%c%c%c%c%c%c\n",
            nn, u64, 
            integerToLetter(cc[ii+0]),
            integerToLetter(cc[ii+1]),
            integerToLetter(cc[ii+2]),
            integerToLetter(cc[ii+3]),
            integerToLetter(cc[ii+4]),
            integerToLetter(cc[ii+5]),
            integerToLetter(cc[ii+6]),
            integerToLetter(cc[ii+7]));
  }
  fclose(F);

  F = fopen(argv[1], "r");
  for (uint32 ii=0, nn=0; ii<ccLen; ii += 4, nn++) {
    uint32 u32;

    fread(&u32, sizeof(uint32), 1, F);
    assert(u32 == *((uint32 *)(cc + ii)));

    fprintf(stdout, "uint32 u32_%02u = 0x%08xlu;  //  %c%c%c%c\n",
            nn, u32,
            integerToLetter(cc[ii+0]),
            integerToLetter(cc[ii+1]),
            integerToLetter(cc[ii+2]),
            integerToLetter(cc[ii+3]));
  }
  fclose(F);

  delete [] cc;

  return(0);
}
