
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

#include "mt19937ar.H"

int
main(int argc, char **argv) {
  mtRandom  mt;

  if (argc != 4)
    fprintf(stderr, "usage: %s <iterations> <lambda> <rho>\n", argv[0]), exit(1);

  uint32  number  = atoi(argv[1]);
  double  mode    = atof(argv[2]);
  double  scale   = atof(argv[3]);

  for (uint32 ii=0; ii<number; ii++)
    fprintf(stdout, "%f\n", mt.mtRandomExponential(mode, scale));

  exit(0);
}

