
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

#include "ansi-escape-v1.H"
#include "arrays.H"

//  https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences

static uint64   aeLen[32] = {0};
static uint64   aeMax[32] = {0};
static char    *aeStr[32] = { nullptr };
static uint32   aePos     = 31;

namespace merylutil::inline system::inline v1 {

char const *
makeAnsiEscapeSequence(const ansiCode *codes) {
  uint32  ap = aePos = (aePos + 1) % 32;
  uint64 &al = aeLen[ap];
  uint64 &am = aeMax[ap];
  char  *&as = aeStr[ap];

  //  Count length of output string, then allocate space for it.

  al = 4;   //  For \033, [, ..., m, \000
  for (uint32 cp=0; (codes[cp] != ansiCode::END); cp++) {
    al += 4;
  }

  resizeArray(as, 0, am, al, _raAct::doNothing);

  //  Build string.

  al = 0;

  as[al++] = '\033';
  as[al++] = '[';

  for (uint32 cp=0; (codes[cp] != ansiCode::END); cp++) {
    if (cp != 0)
      as[al++] = ';';

    toDec((uint32)codes[cp], as, al, am);
  }

  as[al++] = 'm';
  as[al] = 0;

  assert(al <= am);

  return as;
}

}

