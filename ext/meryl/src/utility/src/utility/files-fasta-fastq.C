
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
#include "arrays.H"

#include <stdarg.h>



void
outputFASTA(FILE *f, char const *s, uint64 sl, uint64 bl, char const *h, ...) {

  //  Make life easier later.

  if (bl == 0)
    bl = sl;

  //  Allocate space for the output.  We'll grow it later; right now we just
  //  need to fit the header line.

  uint64  outlen = 0;
  uint64  outmax = 4096 + sl + sl/bl;
  char   *outstr = new char [outmax];

  //  Build the name string.

  outstr[0] = '>';

  va_list ap;

  va_start(ap, h);
  outlen = 1 + vsprintf(outstr+1, h, ap);
  va_end(ap);

  outstr[outlen++] = '\n';

  //  Allocate more space if needed.

  resizeArray(outstr, outlen, outmax, outlen + sl + sl/bl + 16);

  //  Copy the sequence to the outstr, adding newlines as desired.  Catch
  //  when we run out of space and get more - it shouldn't be much more.

  if (bl == 0) {                    //  Hopefully the usual case, since
    for (uint64 si=0; si<sl; )      //  multi-line fasta are evil.
      outstr[outlen++] = s[si++];
  }
  else {
    for (uint64 si=0; si<sl; ) {
      outstr[outlen++] = s[si++];

      if ((si % bl) == 0)
        outstr[outlen++] = '\n';
    }
  }

  //  If not a newline already, add one.

  if (outstr[outlen-1] != '\n')
    outstr[outlen++] = '\n';

  outstr[outlen] = 0;

  //  Did we fail?

  assert(outlen <= outmax);

  //  Now just dump the buffer to disk.

  writeToFile(outstr, "outputFASTA", outlen, f);

  delete [] outstr;
}



void
outputFASTQ(FILE *f, char  const *s, char  const *q, uint64 sl, char  const *h, ...) {
  uint64  outlen = 0;
  uint64  outmax = 4096 + sl + sl;
  char   *outstr = new char [outmax];

  //  Build the name string.

  outstr[0] = '@';

  va_list ap;

  va_start(ap, h);
  outlen = 1 + vsprintf(outstr+1, h, ap);
  va_end(ap);

  outstr[outlen++] = '\n';

  //  Allocate more space if needed.

  resizeArray(outstr, outlen, outmax, outlen + sl + 1 + 2 + sl + 1 + 16);

  //  Copy sequence and qualities to the outstr.

  for (uint64 ii=0; ii<sl; ii++)
    outstr[outlen++] = s[ii];
  outstr[outlen++] = '\n';

  outstr[outlen++] = '+';
  outstr[outlen++] = '\n';

  for (uint64 ii=0; ii<sl; ii++)
    outstr[outlen++] = q[ii];
  outstr[outlen++] = '\n';

  outstr[outlen] = 0;

  //  Did we fail?

  assert(outlen <= outmax);

  //  Now just dump the buffer to disk.

  writeToFile(outstr, "outputFASTQ", outlen, f);

  delete [] outstr;
}



void
outputFASTQ(FILE *f, char  const *s, uint8 const *q, uint64 sl, char  const *h, ...) {
  uint64  outlen = 0;
  uint64  outmax = 4096 + sl + sl;
  char   *outstr = new char [outmax];

  //  Build the name string.

  outstr[0] = '@';

  va_list ap;

  va_start(ap, h);
  outlen = 1 + vsprintf(outstr+1, h, ap);
  va_end(ap);

  outstr[outlen++] = '\n';

  //  Allocate more space if needed.

  resizeArray(outstr, outlen, outmax, outlen + sl + 1 + 2 + sl + 1 + 16);

  //  Copy sequence and qualities to the outstr.

  for (uint64 ii=0; ii<sl; ii++)
    outstr[outlen++] = s[ii];
  outstr[outlen++] = '\n';

  outstr[outlen++] = '+';
  outstr[outlen++] = '\n';

  for (uint64 ii=0; ii<sl; ii++)       //  Copy and convert to ASCII
    outstr[outlen++] = q[ii] + '!';
  outstr[outlen++] = '\n';

  outstr[outlen] = 0;

  //  Did we fail?

  assert(outlen <= outmax);

  //  Now just dump the buffer to disk.

  writeToFile(outstr, "outputFASTQ", outlen, f);

  delete [] outstr;
}



void
outputFASTQ(FILE *f, char  const *s, uint8 qv, uint64 sl, char  const *h, ...) {
  uint64  outlen = 0;
  uint64  outmax = 4096 + sl + sl;
  char   *outstr = new char [outmax];

  //  Build the name string.

  outstr[0] = '@';

  va_list ap;

  va_start(ap, h);
  outlen = 1 + vsprintf(outstr+1, h, ap);
  va_end(ap);

  outstr[outlen++] = '\n';

  //  Allocate more space if needed.

  resizeArray(outstr, outlen, outmax, outlen + sl + 1 + 2 + sl + 1 + 16);

  //  Copy sequence and qualities to the outstr.

  for (uint64 ii=0; ii<sl; ii++)
    outstr[outlen++] = s[ii];
  outstr[outlen++] = '\n';

  outstr[outlen++] = '+';
  outstr[outlen++] = '\n';

  for (uint64 ii=0; ii<sl; ii++)       //  Copy and convert to ASCII
    outstr[outlen++] = qv + '!';
  outstr[outlen++] = '\n';

  outstr[outlen] = 0;

  //  Did we fail?

  assert(outlen <= outmax);

  //  Now just dump the buffer to disk.

  writeToFile(outstr, "outputFASTQ", outlen, f);

  delete [] outstr;
}



void
outputSequence(FILE        *f,
               char  const *outputName,
               char  const *outputBases,
               uint8 const *outputQuals,  uint32  outputLen,
               bool         hasQuals,
               bool         asFASTA,
               bool         asFASTQ,
               uint8        QV) {

  if (asFASTA == asFASTQ)  {
    asFASTA = (hasQuals == false);
    asFASTQ = (hasQuals == true);
  }

  if      ((asFASTQ == true) && (hasQuals == true))
    outputFASTQ(f, outputBases, outputQuals, outputLen, outputName);

  else if ((asFASTQ == true) && (hasQuals == false))
    outputFASTQ(f, outputBases, QV, outputLen, outputName);

  else
    outputFASTA(f, outputBases, outputLen, 0, outputName);
}
