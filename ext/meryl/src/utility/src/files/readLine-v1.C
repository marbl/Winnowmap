
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
#include "arrays.H"
#include "reading-v1.H"

//
//  Two implementations are provided.
//    v0 uses fgets() and is not recently tested.
//    v1 uses getc() and is expected to be faster.
//  This file contains v1.
//

namespace merylutil::inline files::inline v1 {

bool
readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F) {

  if ((L == nullptr) || (Lmax == 0))
    allocateArray(L, Lmax, 1024);

  Llen = 0;
  L[0] = 0;

  if (F == nullptr)
    return(false);

  int32   ch     = getc(F);
  uint32  growth = 1024;

  if (feof(F))
    return false;

  //  Keep reading characters until EOF or a line terminator is encountered.

  while ((feof(F) == 0) && (ferror(F) == 0) && (ch != '\n')) {
    if (Llen + 1 >= Lmax)
      resizeArray(L, Llen, Lmax, Lmax + growth, _raAct::copyData | _raAct::clearNew);  //  Grow the array.

    L[Llen++] = ch;

    ch = getc(F);
  }

  L[Llen] = 0;

  //  Trim trailing whitespace.

  while ((Llen > 0) && (isWhiteSpace(L[Llen-1])))
    L[--Llen] = 0;

  //  Report errors.

  if (ferror(F))
    fprintf(stderr, "ERROR: readLine() got error '%s'.\n", strerror(errno)), exit(1);

  return (ferror(F) == 0);
}

}  //  namespace merylutil::files::v1
