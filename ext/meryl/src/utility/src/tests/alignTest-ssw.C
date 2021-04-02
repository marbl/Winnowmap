
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
#include "types.H"
#include "files.H"
#include "sequence.H"

#include "align-ssw.H"
#include "align-ssw-driver.H"


void
native_interface(char *seqA,
                 char *seqB) {

  //  Convert the sequences to integers.

  uint32   lenA = strlen(seqA);
  uint32   lenB = strlen(seqB);

  int8    *intA = new int8 [lenA];
  int8    *intB = new int8 [lenB];

  for (uint32 ii=0; ii<lenA; ii++)
    intA[ii] = encode2bitBase(seqA[ii]);

  for (uint32 ii=0; ii<lenB; ii++)
    intB[ii] = encode2bitBase(seqB[ii]);

  //

  //                        A   C   T   G
  int8 scoreMatrix[16] = {  1, -1, -1, -1,    // A
                           -1,  1, -1, -1,    // C
                           -1, -1,  1, -1,    // T
                           -1, -1, -1,  1 };  // G

  // flags
  //
  // bit 3:
  // function ssw_align will return the best percent ident alignment
  // beginning position
  //
  // bit 2:
  // if (ref_end1 - ref_begin1 < filterd && read_end1 - read_begin1 < filterd)
  // the function will return the best alignment beginning position and cigar
  //
  // bit 1:
  // if the best alignment score >= filters, the function will return the
  // best alignment beginning position and cigar
  //
  // bit 0:
  // the function will always return the best alignment beginning position
  // and cigar.
  //
  // When no bits set, only the optimal and sub-optimal scores and the
  // optimal alignment ending position will be returned.

  //  Init on the 'read'
  //  Align on the 'ref'

  double      tBgn = 0;
  double      tEnd = 0;

  s_profile *profile = nullptr;
  s_align   *result  = nullptr;

  tBgn = getTime();
  for (uint32 ii=0; ii<20; ii++) {
    profile = ssw_init(intA, lenA, scoreMatrix, 4, 1);
    result  = ssw_align(profile,
                        intB, lenB, 
                        2,     //  gap open (absolute value)
                        1,     //  gap extend (absolute value)
                        1,     //  flags
                        0,     //  filter score    if flag1 and not flag0
                        0,     //  filter distance if flag2 and not flag0
                        1000); //  mask len
  }
  tEnd = getTime();

  fprintf(stdout, "time %f\n", (tEnd - tBgn) / 20);

  fprintf(stdout, "A: 0 - %5d - %5d - %5d\n",                 result->read_begin1, result->read_end1, lenA);
  fprintf(stdout, "B: 0 - %5d - %5d - %5d  best score %5d\n", result->ref_begin1,  result->ref_end1,  lenB, result->score1);
  fprintf(stdout, "B: 0 - %5s - %5d - %5d   2nd score %5d\n", "",                  result->ref_end2,  lenB, result->score2);


  fprintf(stdout, "cigarLen    %d\n", result->cigarLen);

  for (int32 i = 0; i < result->cigarLen; ++i)
    fprintf(stdout, "%d%c", result->cigar[i]>>4, "MID"[result->cigar[i]&0xf]);
  fprintf(stdout, "\n");
}




int
main(int argc, char **argv) {
  char    *seqA = nullptr;
  char    *seqB = nullptr;

#if 0
  fprintf(stderr, "A -> %2u -> %c\n", encode2bitBase('A'), decode2bitBase(0));
  assert(encode2bitBase('A') == 0);
  assert(decode2bitBase(0) == 'A');

  fprintf(stderr, "C -> %2u -> %c\n", encode2bitBase('C'), decode2bitBase(1));
  assert(encode2bitBase('C') == 1);
  assert(decode2bitBase(1) == 'C');

  fprintf(stderr, "T -> %2u -> %c\n", encode2bitBase('T'), decode2bitBase(2));
  assert(encode2bitBase('T') == 2);
  assert(decode2bitBase(2) == 'T');

  fprintf(stderr, "G -> %2u -> %c\n", encode2bitBase('G'), decode2bitBase(3));
  assert(encode2bitBase('G') == 3);
  assert(decode2bitBase(3) == 'G');

  fprintf(stderr, "N -> %2u -> %c\n", encode2bitBase('N'), decode2bitBase(4));
  assert(encode2bitBase('N') == 4);
  assert(decode2bitBase(4) == 'N');
#endif

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-a") == 0) {
      seqA = argv[++arg];
    }

    else if (strcmp(argv[arg], "-b") == 0) {
      seqB = argv[++arg];
    }

    else {
      err++;
    }

    arg++;
  }

  if (seqA == nullptr)  err++;
  if (seqB == nullptr)  err++;

  if (err > 0) {
    fprintf(stderr, "usage: %s -a <bases> -b <bases>\n", argv[0]);
    exit(1);
  }



  sswLib  *ssw = new sswLib(1, -2, -2, -2);

  ssw->align(seqA, strlen(seqA),
             seqB, strlen(seqB));


  delete ssw;
  return(0);
}



