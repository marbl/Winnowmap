
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


int
main(int argc, char **argv) {
  char na[32] = { 0 };
  char nb[32] = { 0 };

  strcpy(na, "./timeTestA-XXXXXXXXXXXXXXXX");
  mktemp(na);

  strcpy(nb, "./timeTestB-XXXXXXXXXXXXXXXX");
  mktemp(nb);

  merylutil::muDeltaTime  t(true);

  fprintf(stderr, "Create delta:    %f\n", t.seconds());
  fprintf(stderr, "Create fulltime: %f\n", t.full_seconds());

  fprintf(stderr, "create empty na: '%s'\n", na);
  merylutil::createEmptyFile(na);
  sleep(3);
  fprintf(stderr, "create empty nb: '%s'\n", nb);
  merylutil::createEmptyFile(nb);

  merylutil::muFileTime t1;
  merylutil::muFileTime t2;

  t1.getTimeOfFile(na);
  t2.getTimeOfFile(nb);

  fprintf(stderr, "three second old file: %f\n", t1.fileAge().seconds());
  fprintf(stderr, "zero  second old file: %f\n", t2.fileAge().seconds());

  merylutil::unlink(na);
  merylutil::unlink(nb);

  return 0;
}
