
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

#include "kmers.H"

uint32 kmerTiny::_merSize   = 0;
kmdata kmerTiny::_fullMask  = 0;
kmdata kmerTiny::_leftMask  = 0;
uint32 kmerTiny::_leftShift = 0;



char *
constructBlockName(char   *nameprefix,
                   uint64  outIndex,
                   uint32  numFiles,
                   uint32  iteration,
                   bool    isIndex) {
  char *name = new char [FILENAME_MAX+1];
  char  bits[67] = { 0 };

  bits[0] = '0';
  bits[1] = 'x';

  uint32 bp = 2;

  for (uint32 mask=1; mask < numFiles; mask <<= 1)   //  Count the number of digits we need.
    bp++;

  for (uint32 mask=1; mask < numFiles; mask <<= 1)   //  Then make the name from right to left.
    bits[--bp] = (outIndex & mask) ? '1' : '0';

  if (iteration == 0)
    snprintf(name, FILENAME_MAX, "%s/%s.%s", nameprefix, bits, (isIndex == false) ? "merylData" : "merylIndex");
  else
    snprintf(name, FILENAME_MAX, "%s/%s[%03u].%s", nameprefix, bits, iteration, (isIndex == false) ? "merylData" : "merylIndex");

  return(name);
}



FILE *
openOutputBlock(char   *nameprefix,
                uint64  fileIndex,
                uint32  numFiles,
                uint32  iteration) {
  char    *name = constructBlockName(nameprefix, fileIndex, numFiles, iteration, false);

  FILE *F = AS_UTL_openOutputFile(name);

  delete [] name;

  return(F);
}



FILE *
openInputBlock(char   *nameprefix,
               uint64  fileIndex,
               uint32  numFiles,
               uint32  iteration) {
  char    *name = constructBlockName(nameprefix, fileIndex, numFiles, iteration, false);

  FILE *F = AS_UTL_openInputFile(name);

  delete [] name;

  return(F);
}
