
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
//  This file contains v0.
//

namespace merylutil::inline files::v0 {

bool
readLine(char *&L, uint32 &Llen, uint32 &Lmax, FILE *F) {

  if ((L == nullptr) || (Lmax == 0))
    allocateArray(L, Lmax, 4, _raAct::clearNew);

  L[0]      = 0;
  L[Lmax-2] = 0;
  L[Lmax-1] = 0;

  if (F == nullptr)
    return(false);

  fgets(L, Lmax, F);

  Llen = strlen(L);

  //  fgets() will always NUL-terminate the string.  If the seocnd to last
  //  character exists and is not a newline, we didn't read the whole string.

  while ((L[Lmax-2] != 0) && (L[Lmax-2] != '\n')) {
    uint32   growth = 4;

    assert(Llen == Lmax - 1);

    resizeArray(L, Llen, Lmax, Lmax + growth);  //  Grow the array.
    L[Lmax-2] = 0;
    L[Lmax-1] = 0;

    fgets(L + Llen, 1 + growth, F);             //  Read more bytes.

    Llen += strlen(L + Llen);                   //  How many more?
  }

  //  Trim trailing whitespace.

  while ((Llen > 0) && (isWhiteSpace(L[Llen-1])))
    L[--Llen] = 0;

  return(true);
}

}  //  namespace merylutil::files::v0
