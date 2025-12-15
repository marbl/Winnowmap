
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

#include "sequence.H"

using namespace merylutil::sequence;

bool
testLoadSeq(void) {
  FILE *O;

  O = merylutil::openOutputFile("sequenceTest.data.fasta");
  fprintf(O, ">name\n");
  fprintf(O, "ACGT\n");
  fprintf(O, ">  name   \n");
  fprintf(O, "ACGT\n");
  fprintf(O, ">  name   flags\n");
  fprintf(O, "ACGT\n");
  fprintf(O, ">  name   f l a g s    \n");
  fprintf(O, "ACGT\n");
  merylutil::closeFile(O);

  dnaSeqFile *F = openSequenceFile("sequenceTest.data.fasta");
  dnaSeq      S;

  F->loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "")          == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  F->loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "")          == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  F->loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "flags")     == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  F->loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "f l a g s") == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  merylutil::unlink("sequenceTest.data.fasta");

  delete F;

  fprintf(stderr, "Success!\n");

  return true;
}



bool
encodeTest(void) {

  fprintf(stderr, "%c=%02x %c=%02x %c=%02x %c=%02x %c=%02x\n",
          decode2bitBase(0x00), encode2bitBase(decode2bitBase(0x00)),
          decode2bitBase(0x01), encode2bitBase(decode2bitBase(0x01)),
          decode2bitBase(0x02), encode2bitBase(decode2bitBase(0x02)),
          decode2bitBase(0x03), encode2bitBase(decode2bitBase(0x03)),
          decode2bitBase(0x04), encode2bitBase(decode2bitBase(0x04)));

  assert(encode2bitBase(decode2bitBase(0x00)) == 0x00);
  assert(encode2bitBase(decode2bitBase(0x01)) == 0x01);
  assert(encode2bitBase(decode2bitBase(0x02)) == 0x02);
  assert(encode2bitBase(decode2bitBase(0x03)) == 0x03);
  assert(encode2bitBase(decode2bitBase(0x04)) == 0x04);

  assert(encode2bitBase('a') == 0x00);
  assert(encode2bitBase('c') == 0x01);
  assert(encode2bitBase('t') == 0x02);
  assert(encode2bitBase('g') == 0x03);
  assert(encode2bitBase('n') == 0x04);

  assert(encode2bitBase('A') == 0x00);
  assert(encode2bitBase('C') == 0x01);
  assert(encode2bitBase('T') == 0x02);
  assert(encode2bitBase('G') == 0x03);
  assert(encode2bitBase('N') == 0x04);

  assert(decode2bitBase(0x00) == 'A');
  assert(decode2bitBase(0x01) == 'C');
  assert(decode2bitBase(0x02) == 'T');
  assert(decode2bitBase(0x03) == 'G');
  assert(decode2bitBase(0x04) == 'N');


  fprintf(stdout, "\n");
  fprintf(stdout, "     ASCII                    encode2bitBase()         decode2bitBase()\n");
  fprintf(stdout, "   \\_0123456789abcdef       \\_0123456789abcdef       \\_0123456789abcdef\n");
  for (uint32 ii=0; ii<16; ii++) {
    fprintf(stdout, "0x%01x| ", ii);
    for (uint32 jj=0; jj<16; jj++)
      fprintf(stdout, "%c", isVisible(ii*16 + jj) ? (ii*16 + jj) : '.');

    fprintf(stdout, "    0x%01x| ", ii);
    for (uint32 jj=0; jj<16; jj++)
      fprintf(stdout, "%c", encode2bitBase(ii*16 + jj) + '0');

    fprintf(stdout, "    0x%01x| ", ii);
    for (uint32 jj=0; jj<16; jj++)
      fprintf(stdout, "%c", decode2bitBase(ii*16 + jj));
    fprintf(stdout, "\n");
  }


  return true;
}



int
main(int argc, char **argv) {

  int arg=1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-l") == 0) {
      testLoadSeq();
    }
    else if (strcmp(argv[arg], "-e") == 0) {
      encodeTest();
    }
    else {
      fprintf(stderr, "usage: %s -l | -t\n", argv[0]);
      return 1;
    }

    arg++;
  }



  return 0;
}
