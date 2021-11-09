
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
#include "files.H"

uint32  lengths1[] = {
  12, 100,
  12, 100,
  12, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
  12, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 1,
  12, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  12, 10,
  0,
  15, 100, 1, 100,
  15, 100, 1, 100,
  15, 100, 1, 100,
  15, 100, 1, 100,
  0,
  15, 100, 1, 100,
  15, 100, 1, 100,
  15, 100, 1, 100,
  15, 100, 1, 100,
};

uint32  lengths2[] = {
  2, 10,
  2, 10, 1, 10,
  2, 10,
  2, 10,
  2, 10, 1, 10,
  2, 10, 1, 10,
  2, 10,
  2, 10, 1, 10
};


bool
checkLineLength(char const *filename, uint32 *linelength) {
  uint32   Lnum = 0;
  uint32   Llen = 0;
  uint32   Lmax = 0;
  char    *L    = nullptr;
  bool     pass = true;

  FILE *F = AS_UTL_openInputFile(filename);

  if (linelength == nullptr)
    fprintf(stdout, "uint32 *lengths = {\n");

  while (AS_UTL_readLine(L, Llen, Lmax, F) == true) {
    if (linelength == nullptr)
      fprintf(stdout, "%u,\n", Llen);
    else
      pass &= (Llen == linelength[Lnum]);

    Lnum++;
  }

  AS_UTL_closeFile(F);

  delete [] L;

  return(pass);
}


int
main(int argc, char **argv) {
  char   *seq = new char  [256];
  uint8  *qvs = new uint8 [256];
  char   *qlt = new char  [256];

  FILE   *F   = nullptr;
  bool    t1  = false;
  bool    t2  = false;

  //  Wikipedia says encoded QVs run from '!' (dec 33) to '~' (dec 126) for a
  //  range of 94.
  //
  //  Special obnoxious values are
  //    '@' dec 64 qvs 31
  //    '+' dec 43 qvc 10
  //    '>' dec 62 qvc 29

  for (uint32 ii=0; ii<256; ii++) {
    seq[ii] = 'A';
    qvs[ii] = ii % 94;
    qlt[ii] = '!' + qvs[ii];
  }

  {
    FILE *F = AS_UTL_openOutputFile("fasta-fastq.1.test");

    outputFASTA(F, seq, 100,   0, "name%s%d.", "_one_", 1);
    outputFASTA(F, seq, 100, 100, "name%s%d.", "_two_", 2);
    outputFASTA(F, seq, 100,  10, "name%s%d.", "_thr_", 3);
    outputFASTA(F, seq, 100,   9, "name%s%d.", "_for_", 4);
    outputFASTA(F, seq,  10,   1, "name%s%d.", "_fiv_", 5);
    outputFASTA(F, seq,  10,   0, "name%s%d.", "_six_", 6);

    fprintf(F, "\n");

    outputFASTQ(F, seq,    qvs,    100, "name%s%d.", "_oneqvs_", 1);
    outputFASTQ(F, seq+31, qvs+31, 100, "name%s%d.", "_twoqvs_", 2);
    outputFASTQ(F, seq+10, qvs+10, 100, "name%s%d.", "_thrqvs_", 3);
    outputFASTQ(F, seq+29, qvs+29, 100, "name%s%d.", "_forqvs_", 4);

    fprintf(F, "\n");

    outputFASTQ(F, seq,    qlt,    100, "name%s%d.", "_oneqlt_", 1);
    outputFASTQ(F, seq+31, qlt+31, 100, "name%s%d.", "_twoqlt_", 2);
    outputFASTQ(F, seq+10, qlt+10, 100, "name%s%d.", "_thrqlt_", 3);
    outputFASTQ(F, seq+29, qlt+29, 100, "name%s%d.", "_forqlt_", 4);

    AS_UTL_closeFile(F, "fasta-fastq.1.test");

    t1 = checkLineLength("fasta-fastq.1.test", lengths1);
  }

  {
    FILE *F = AS_UTL_openOutputFile("fasta-fastq.2.test");
    outputSequence(F, "1", seq, qvs, 10, false, false, false, 0);
    outputSequence(F, "2", seq, qvs, 10, false, false,  true, 0);
    outputSequence(F, "3", seq, qvs, 10, false,  true, false, 0);
    outputSequence(F, "4", seq, qvs, 10, false,  true,  true, 0);
    outputSequence(F, "5", seq, qvs, 10,  true, false, false, 0);
    outputSequence(F, "6", seq, qvs, 10,  true, false,  true, 0);
    outputSequence(F, "7", seq, qvs, 10,  true,  true, false, 0);
    outputSequence(F, "8", seq, qvs, 10,  true,  true,  true, 0);
    AS_UTL_closeFile(F, "fasta-fastq.2.test");

    t2 = checkLineLength("fasta-fastq.2.test", lengths2);
  }

  if (t1 == false)
    fprintf(stderr, "TEST 1 FAILED.\n");

  if (t2 == false)
    fprintf(stderr, "TEST 2 FAILED.\n");

  if ((t1 == true) &&
      (t2 == true)) {
    fprintf(stderr, "Success!\n");
    AS_UTL_unlink("fasta-fastq.1.test");
    AS_UTL_unlink("fasta-fastq.2.test");
    return(0);
  }

  return(1);
}
