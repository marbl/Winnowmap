
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

#include "system.H"
#include "types.H"
#include "files.H"
#include "sequence.H"

#include "align-ksw2-driver.H"



#if 0
void
print_aln(const char *qname, const char *tname, ksw_extz_t *ez) {

  fprintf(stdout, "%s %5d-%5d\n", qname, 0, ez->max_q);
  fprintf(stdout, "%s %5d-%5d\n", tname, 0, ez->max_t);
  fprintf(stdout, "score %d max %d\n", ez->score, ez->max);

	if (ez->n_cigar > 0) {
		for (int32 i = 0; i < ez->n_cigar; ++i)
			fprintf(stdout, "%d%c", ez->cigar[i]>>4, "MID"[ez->cigar[i]&0xf]);
    fprintf(stdout, "\n");
	}
}
#endif




int
main(int argc, char **argv) {
  dnaSeqFile  *fileA, *fileB;
  dnaSeq       dseqA,  dseqB;
  char const  *seqA = nullptr;
  char const  *seqB = nullptr;

  //fprintf(stderr, "A -> %2u -> %c\n", encode2bitBase('A'), decode2bitBase(0));
  assert(encode2bitBase('A') == 0);
  assert(decode2bitBase(0) == 'A');

  //fprintf(stderr, "C -> %2u -> %c\n", encode2bitBase('C'), decode2bitBase(1));
  assert(encode2bitBase('C') == 1);
  assert(decode2bitBase(1) == 'C');

  //fprintf(stderr, "G -> %2u -> %c\n", encode2bitBase('G'), decode2bitBase(3));
  assert(encode2bitBase('G') == 3);
  assert(decode2bitBase(3) == 'G');

  //fprintf(stderr, "T -> %2u -> %c\n", encode2bitBase('T'), decode2bitBase(2));
  assert(encode2bitBase('T') == 2);
  assert(decode2bitBase(2) == 'T');

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-a") == 0) {
      seqA = argv[++arg];
    }

    else if (strcmp(argv[arg], "-b") == 0) {
      seqB = argv[++arg];
    }

    if      (strcmp(argv[arg], "-A") == 0) {
      fileA = new dnaSeqFile(argv[++arg]);
      fileA->loadSequence(dseqA);
      seqA = dseqA.bases();
    }

    else if (strcmp(argv[arg], "-B") == 0) {
      fileB = new dnaSeqFile(argv[++arg]);
      fileB->loadSequence(dseqB);
      seqB = dseqB.bases();
    }

    else {
      err++;
    }

    arg++;
  }

  if (seqA == nullptr)  err++;
  if (seqB == nullptr)  err++;

  if (err > 0) {
    fprintf(stderr, "usage: %s -a <bases> -b <bases>\n", argv[0]);
    fprintf(stderr, "       %*s -A <file>  -B <file>\n", (int)strlen(argv[0]), "");
    exit(1);
  }

  ksw2Lib   ksw(2, -4, 1, 1);

  ksw.align(seqA, strlen(seqA), 0, strlen(seqA),
            seqB, strlen(seqB), 0, strlen(seqB), true);

  fprintf(stdout, "identity: %f\n", ksw.percentIdentity());
  fprintf(stdout, "score:    %d\n", ksw.score());

  return(0);
}
