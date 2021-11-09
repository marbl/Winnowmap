
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

#include "files.H"

int32
main(int32 argc, char **argv) {
  uint32   lineMax = 0;
  uint32   lineLen = 0;
  char    *line    = nullptr;
  uint32   nLines  = 0;

  if (argc == 1) {
    fprintf(stderr, "usage: %s inputFile[.gz]\n", argv[0]);
    return(1);
  }

  compressedFileReader  *in = new compressedFileReader(argv[1]);

  while (AS_UTL_readLine(line, lineLen, lineMax, in->file())) {
    nLines++;
  }

  delete    in;
  delete [] line;

  fprintf(stderr, "Found %u lines!  Yay!\n", nLines);

  return(0);
}
