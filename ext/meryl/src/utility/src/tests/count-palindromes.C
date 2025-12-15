
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
#include "kmers.H"

using namespace merylutil;
using namespace merylutil::kmers::v2;


void
showScanner(uint64 nKmers, uint64 nKmersExpected, uint64 nPalindrome) {
  uint64  S = nKmersExpected / 77;

  if (S < 1000000)   //  Don't show if small.
    return;

  uint64  w = nKmers / S;
  uint64  m = nKmersExpected / S;

  char    s[128] = {0};

  for (uint32 ii=0; ii<m; ii++)
    s[ii] = (ii < w) ? '+' : '-';

  if (nKmers < nKmersExpected)
    fprintf(stderr, "Scanning: [%s] %6.2f%%\r", s, 100.0 * nKmers / nKmersExpected);
  else
    fprintf(stderr, "Scanned!  [%s] %6.2f%%\n", s, 100.0 * nKmers / nKmersExpected);
}


int
main(int argc, char **argv) {
  uint32   k = 0;
  kmer     m;
  kmer     c;

  std::vector<char const *>  err;
  for (int32 ii=1; ii<argc; ii++) {
    if (strcmp(argv[ii], "-k") == 0) {
      k = strtouint32(argv[++ii]);
    }
    else {
      sprintf(err, "Unknown option '%s'.\n", argv[ii]);
    }
  }
  if ((err.size() > 0) || (k == 0)) {
    fprintf(stderr, "usage: %s -k K\n", argv[0]);
    fprintf(stderr, "  Explicitly counts the number of palindromic k-mers as a test\n");
    fprintf(stderr, "  of kmer iteration, reverse-complement, and progress reporting.\n");
    fprintf(stderr, "  Otherwise not interesting.  Performance is around 75 Mkmer/sec.\n");
    fprintf(stderr, "  K more than 17 is not recommended.\n");
    for (char const *e : err)
      fprintf(stderr, "ERROR: %s", e);
    exit(1);
  }

  kmer::setSize(k);

  uint64   nKmers      = 0, nKmersExpected      = uint64one << (2*k);
  uint64   nPalindrome = 0, nPalindromeExpected = ((k % 2) == 0) ? (uint64one << k) : 0;
  uint64   nSmaller    = 0, nSmallerExpected    = nKmersExpected/2 - nPalindromeExpected/2;

  showScanner(nKmers, nKmersExpected, nPalindrome);

  do {
    c = m;
    c.reverseComplement();

    nKmers++;

    if (c == m)
      nPalindrome++;

    if (m < c)
      nSmaller++;

    if ((nKmers & 0x1ffffffllu) == 0x00)
      showScanner(nKmers, nKmersExpected, nPalindrome);
  } while (m++.isLast() == false);

  fprintf(stdout, "\n");
  fprintf(stdout, "K            %12u\n", kmer::merSize());
  fprintf(stdout, "nKmers       %12lu expecting %12lu\n", nKmers,      nKmersExpected);
  fprintf(stdout, "nPalindrome  %12lu expecting %12lu\n", nPalindrome, nPalindromeExpected);
  fprintf(stdout, "nSmaller     %12lu expecting %12lu\n", nSmaller,    nSmallerExpected);

  assert(nKmers      == nKmersExpected);
  assert(nPalindrome == nPalindromeExpected);
  assert(nSmaller    == nSmallerExpected);

  fprintf(stdout, "\n");
  fprintf(stdout, "Pass!\n");

  return 0;
}
