
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

namespace merylutil::inline kmers::v2 {

//  Like loadBlock, but just reports all blocks in the file, ignoring
//  the kmer data.
//
void
dumpMerylDataFile(char *name) {
  FILE            *F = NULL;
  merylFileIndex   I;
  stuffedBits     *D = NULL;

  //  Dump the merylIndex for this block.

  if (fileExists(name, '.', "merylIndex") == false)
    fprintf(stderr, "ERROR: '%s.merylIndex' doesn't exist.  Can't dump it.\n",
            name), exit(1);

  F = merylutil::openInputFile(name, '.', "merylIndex");

  fprintf(stdout, "\n");
  fprintf(stdout, "    prefix    blkPos    nKmers\n");
  fprintf(stdout, "---------- --------- ---------\n");

  while (loadFromFile(I, "merylFileIndex", F, false) != 0) {
    fprintf(stdout, "0x%08x %9lu %9lu\n", I.blockPrefix(), I.blockPosition(), I.numKmers());
  }

  merylutil::closeFile(F);

  //  Read each block, sequentially, and report the header.

  if (fileExists(name, '.', "merylData") == false)
    fprintf(stderr, "ERROR: '%s.merylData' doesn't exist.  Can't dump it.\n",
            name), exit(1);

  F = merylutil::openInputFile(name, '.', "merylData");
  D = new stuffedBits;

  fprintf(stdout, "\n");
  fprintf(stdout, "            prefix   nKmers kCode uBits bBits                 k1 cCode                 c1                 c2\n");
  fprintf(stdout, "------------------ -------- ----- ----- ----- ------------------ ----- ------------------ ------------------\n");

  while (D->loadFromFile(F)) {
    uint64 position   = D->getPosition();

    uint64 m1         = D->getBinary(64), M1 = 0x7461446c7972656dllu;
    uint64 m2         = D->getBinary(64), M2 = 0x0a3130656c694661llu;

    uint64 prefix     = D->getBinary(64);
    uint64 nKmers     = D->getBinary(64);

    uint8  kCode      = D->getBinary(8);
    uint32 unaryBits  = D->getBinary(32);
    uint32 binaryBits = D->getBinary(32);
    uint64 k1         = D->getBinary(64);

    uint8  cCode      = D->getBinary(8);
    uint64 c1         = D->getBinary(64);
    uint64 c2         = D->getBinary(64);

    if ((m1 != M1) ||
        (m2 != M2)) {
      fprintf(stderr, "dumpMerylDataFile()-- Magic number mismatch at position " F_U64 ".\n", position);
      fprintf(stderr, "dumpMerylDataFile()-- Expected 0x%016" F_X64P " got\n", M1);
      fprintf(stderr, "dumpMerylDataFile()            0x%016" F_X64P "\n", m1);
      fprintf(stderr, "dumpMerylDataFile()-- Expected 0x%016" F_X64P " got\n", M2);
      fprintf(stderr, "dumpMerylDataFile()            0x%016" F_X64P "\n", m2);
      exit(1);
    }

    fprintf(stdout, "0x%016lx %8lu %5u %5u %5u 0x%016lx %5u 0x%016lx 0x%016lx\n",
            prefix, nKmers, kCode, unaryBits, binaryBits, k1, cCode, c1, c2);
  }

  delete D;

  merylutil::closeFile(F);

  //  Read each block again, dump the kmers in the block.

  F = merylutil::openInputFile(name, '.', "merylData");
  D = new stuffedBits;

  while (D->loadFromFile(F)) {
    uint64 position   = D->getPosition();

    uint64 m1         = D->getBinary(64);
    uint64 m2         = D->getBinary(64);

    uint64 prefix     = D->getBinary(64);
    uint64 nKmers     = D->getBinary(64);

    uint8  kCode      = D->getBinary(8);
    uint32 unaryBits  = D->getBinary(32);
    uint32 binaryBits = D->getBinary(32);
    uint64 k1         = D->getBinary(64);

    uint8  cCode      = D->getBinary(8);
    uint64 c1         = D->getBinary(64);
    uint64 c2         = D->getBinary(64);

    uint8  lCode      = D->getBinary(8);    //  Only in merylDataFile01!
    uint32 labelBits  = D->getBinary(6);
    uint64 l1         = D->getBinary(58);
    uint64 l2         = D->getBinary(64);

    fprintf(stdout, "\n");
    fprintf(stdout, " kmerIdx prefixDelta      prefix |--- suffix-size and both suffixes ---|    value\n");
    fprintf(stdout, "-------- ----------- ----------- -- ---------------- -- ---------------- --------\n");

    uint64   *pd = new uint64 [nKmers];
    uint64   *s1 = new uint64 [nKmers];
    uint64   *s2 = new uint64 [nKmers];
    uint64   *va = new uint64 [nKmers];
    uint64   *la = new uint64 [nKmers];

    uint32    ls = (binaryBits <= 64) ? (0)          : (binaryBits - 64);
    uint32    rs = (binaryBits <= 64) ? (binaryBits) : (64);

    uint64    tp = 0;

    //  Get all the kmers.
    if (kCode == 1)
      for (uint32 kk=0; kk<nKmers; kk++) {
        pd[kk] = D->getUnary();
        s1[kk] = D->getBinary(ls);
        s2[kk] = D->getBinary(rs);
      }

    else {
      fprintf(stderr, "ERROR: unknown kCode %u\n", kCode), exit(1);
    }

    //  Get all the values.
    if      (cCode == 1) {
      for (uint32 kk=0; kk<nKmers; kk++)
        va[kk] = D->getBinary(32);
    }

    else if (cCode == 2) {
      for (uint32 kk=0; kk<nKmers; kk++)
        va[kk] = D->getBinary(64);
    }

    else {
      fprintf(stderr, "ERROR: unknown cCode %u\n", cCode), exit(1);
    }

    //  Get all the labels.
    if      (lCode == 0) {
      for (uint32 kk=0; kk<nKmers; kk++)
        la[kk] = 0;
    }

    else if (lCode == 1) {
      for (uint32 kk=0; kk<nKmers; kk++)
        la[kk] = D->getBinary(labelBits);
    }

    else {
      fprintf(stderr, "ERROR: unknown lCode 0x%02x\n", lCode), exit(1);
    }

    //  Dump.
    for (uint32 kk=0; kk<nKmers; kk++) {
      tp += pd[kk];

      fprintf(stdout, "%8u %11lu %011lx %2u %016lx %2u %016lx %8lx\n",
              kk, pd[kk], tp, ls, s1[kk], rs, s2[kk], va[kk]);
    }

    delete [] la;
    delete [] va;
    delete [] s2;
    delete [] s1;
    delete [] pd;
  }

  delete D;

  merylutil::closeFile(F);
}

}  //  namespace merylutil::kmers::v2
