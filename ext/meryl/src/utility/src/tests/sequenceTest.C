
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

int
main(int argc, char **argv) {
  FILE *O;

  O = AS_UTL_openOutputFile("sequenceTest.data.fasta");
  fprintf(O, ">name\n");
  fprintf(O, "ACGT\n");
  fprintf(O, ">  name   \n");
  fprintf(O, "ACGT\n");
  fprintf(O, ">  name   flags\n");
  fprintf(O, "ACGT\n");
  fprintf(O, ">  name   f l a g s    \n");
  fprintf(O, "ACGT\n");
  AS_UTL_closeFile(O);

  dnaSeqFile  F("sequenceTest.data.fasta");
  dnaSeq      S;

  F.loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "")          == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  F.loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "")          == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  F.loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "flags")     == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  F.loadSequence(S);
  assert(strcmp(S.ident(), "name")      == 0);
  assert(strcmp(S.flags(), "f l a g s") == 0);
  assert(strcmp(S.bases(), "ACGT")      == 0);

  AS_UTL_unlink("sequenceTest.data.fasta");

  fprintf(stderr, "Success!\n");

  return(0);
}

