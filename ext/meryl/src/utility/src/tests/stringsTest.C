
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
 
#include "strings.H"

int
main(int argc, char **argv) {
  splitToWords  W;
  splitType     type = splitWords;

  for (uint32 arg=1; arg<argc; arg++) {
    if (strcmp(argv[arg], "-p") == 0) {
      type = splitPaths;
      continue;
    }

    if (strcmp(argv[arg], "-w") == 0) {
      type = splitWords;
      continue;
    }

    W.split(argv[arg], type);

    fprintf(stderr, "'%s'\n", argv[arg]);

    for (uint32 ii=0; ii<W.numWords(); ii++)
      fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
  }

  exit(0);
}
