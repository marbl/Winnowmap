
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
#include "strings.H"

using merylutil::readBuffer;
using merylutil::displayLetter;


int
main(int argc, char **argv) {

  readBuffer *B = new readBuffer("readBufferTest.txt");

  uint64      Llen = 0;
  uint64      Lmax = 1024;
  char       *L    = new char [Lmax];

  fprintf(stderr, "skip white\n");
  B->skipWhitespace();
  fprintf(stderr, "peek %5s @ %5lu -- expect 'E'\n", displayLetter(B->peek()), B->tell());

  B->copyVisible(L, Llen, Lmax);
  fprintf(stderr, "read '%s' -- expect 'EOL'\n", L);  Llen=0;
  fprintf(stderr, "peek %5s @ %5lu -- expect LF\n", displayLetter(B->peek()), B->tell());

  fprintf(stderr, "skip white\n");
  B->skipWhitespace();
  fprintf(stderr, "peek %5s @ %5lu -- expect 'c'\n", displayLetter(B->peek()), B->tell());

  fprintf(stderr, "skip visible\n");
  B->skipVisible();
  fprintf(stderr, "peek %5s @ %5lu -- expect CR\n", displayLetter(B->peek()), B->tell());

  //-- so really want 'skipLine and all CR/LF's after it.
  //-- or go back to the original find first \n and assume we only get \r\n?

  fprintf(stderr, "skip eol\n");
  B->skipLine();
  B->skipLine();
  B->skipLine();
  B->skipLine();
  B->skipLine();
  B->skipLine();
  B->skipLine();
  fprintf(stderr, "peek %5s @ %5lu -- expect space\n", displayLetter(B->peek()), B->tell());

  fprintf(stderr, "skip white\n");
  B->skipWhitespace();
  fprintf(stderr, "peek %5s @ %5lu -- expect 'a'\n", displayLetter(B->peek()), B->tell());

  B->copyVisible(L, Llen, Lmax);
  fprintf(stderr, "read '%s'\n", L);  Llen=0;
  fprintf(stderr, "peek %5s @ %5lu\n", displayLetter(B->peek()), B->tell());

  fprintf(stderr, "skip line, copy visible\n");
  B->skipLine();
  B->copyVisible(L, Llen, Lmax);
  fprintf(stderr, "read '%s' -- expect 'EOL'\n", L);  Llen=0;
  fprintf(stderr, "peek %5s @ %5lu\n", displayLetter(B->peek()), B->tell());

  delete B;
}


