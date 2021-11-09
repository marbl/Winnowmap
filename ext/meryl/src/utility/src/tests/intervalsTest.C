
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

#include "runtime.H"
#include "intervals.H"

#include "mt19937ar.H"

void
boringTest(void) {
  bool              errors = false;
  intervals<int32>  t1;

  t1.add_span(11, -4);
  t1.add_position(0, 10);
  t1.add_span(8, 12);

  errors |= ((t1.size() != 3) ||
             (t1.bgn(0) != 7) || (t1.end(0) != 11) ||
             (t1.bgn(1) != 0) || (t1.end(1) != 10) ||
             (t1.bgn(2) != 8) || (t1.end(2) != 20));

  if (errors) {
    fprintf(stderr, "BEFORE:\n");
    for (uint32 ii=0; ii<t1.size(); ii++)
      fprintf(stderr, "%2d %3d-%3d\n", ii, t1.bgn(ii), t1.end(ii));
  }

  t1.squash(-1);

  errors |= ((t1.size() != 1) ||
             (t1.bgn(0) != 0) || (t1.end(0) != 20));

  if (errors) {
    fprintf(stderr, "AFTER:\n");
    for (uint32 ii=0; ii<t1.size(); ii++)
      fprintf(stderr, "%2d %3d-%3d\n", ii, t1.bgn(ii), t1.end(ii));
  }

  if (errors)
    fprintf(stderr, "FAIL.\n");
  else
    fprintf(stderr, "Success!\n");
}



void
invertTest(void) {
  bool              errors = false;
  intervals<int32>  t1;
  intervals<int32>  t2;

  t1.add_position(-30, -10);
  t1.add_position( -5,  5);
  t1.add_position( 10,  30);

  t2.setToInversion(-20, 20, t1);

  errors |= ((t2.size() != 4) ||
             (t2.bgn(0) != -10) || (t2.end(0) !=  20) ||
             (t2.bgn(1) != -20) || (t2.end(1) != -5) ||
             (t2.bgn(2) !=   5) || (t2.end(2) !=  20) ||
             (t2.bgn(3) != -20) || (t2.end(3) !=  10));

  if (errors) {
    fprintf(stderr, "BEFORE:\n");
    for (uint32 ii=0; ii<t1.size(); ii++)
      fprintf(stderr, "%2d %3d-%3d\n", ii, t1.bgn(ii), t1.end(ii));

    fprintf(stderr, "AFTER:\n");
    for (uint32 ii=0; ii<t2.size(); ii++)
      fprintf(stderr, "%2d %3d-%3d\n", ii, t2.bgn(ii), t2.end(ii));
  }

  if (errors)
    fprintf(stderr, "FAIL.\n");
  else
    fprintf(stderr, "Success!\n");
}



void
expensiveTest(uint32 seed) {
  mtRandom  mt(seed);

  //  About 6.5 minutes per million, so this should be about an hour.
  uint32  iterMax = 935000;

  for (uint32 iter=0; iter<iterMax; iter++) {
    uint32  numIntervals =     mt.mtRandom32() % 5000;
    uint32  maxLen       = 1 + mt.mtRandom32() % 1000;
    uint32  maxBgn       = 1 + mt.mtRandom32() % 50000;
    uint32 *depth        = new uint32 [maxBgn + maxLen];

    memset(depth, 0, sizeof(uint32) * (maxBgn + maxLen));

    if (iter % 1000 == 0)
      fprintf(stderr, "%10u/%10u: %3u intervals, each up to %4u long, coords up to %4u\n",
              iter, iterMax,
              numIntervals, maxLen, maxBgn);

    intervals<uint32>  il;

    //  Add intervals to the list.
    //  Sum depths explicitly.
    for (uint32 ii=0; ii<numIntervals; ii++) {
      uint32  bgn = mt.mtRandom32() % maxBgn;   //  bgn between 0 and maxBgn
      uint32  len = mt.mtRandom32() % maxLen;   //  len between 0 and maxLen
      uint32  end = bgn + len;

      if (mt.mtRandom32() < uint32max / 2)
        il.add_span(bgn, len);
      else
        il.add_position(bgn, end);

      for (uint32 xx=bgn; xx<end; xx++)
        depth[xx]++;

      //fprintf(stderr, "IL %u - %u\n", bgn, bgn+len);
    }

    //  Convert intervals to depths.

    intervalsDepth<uint32>  de(il);

    //  Over all the depth regions, subtract the computed depth from
    //  the explicit depth.
    for (uint32 xx=0; xx<de.size(); xx++) {
      uint32  bgn = de.bgn(xx);
      uint32  end = de.end(xx);
      uint32  dpt = de.depth(xx);

      //fprintf(stderr, "ID %u - %u depth %u\n", bgn, end, dpt);

      for (uint32 cc=bgn; cc<end; cc++) {
        //if (cc < 30)
        //  fprintf(stderr, "depth[%u] = %u -> %u\n", cc, depth[cc], depth[cc] - dpt);
        depth[cc] -= dpt;
      }
    }

    //  Every explicit depth should now be zero, even the ones
    //  not covered.
    for (uint32 cc=0; cc<maxBgn + maxLen; cc++) {
      if (depth[cc] != 0)
        fprintf(stderr, "ERROR: depth[%u] = %u in iter %u\n", cc, depth[cc], iter);
      assert(depth[cc] == 0);
    }


    delete [] depth;
  }

  fprintf(stderr, "Success!\n");
}



int
main(int argc, char **argv) {
  bool    doBoring      = false;
  bool    doInvert      = false;
  bool    doExpensive   = false;
  uint32  expensiveSeed = 9;

  AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-boring") == 0) {
      doBoring = true;
    }

    else if (strcmp(argv[arg], "-invert") == 0) {
      doInvert = true;
    }

    else if (strcmp(argv[arg], "-expensive") == 0) {
      doExpensive = true;

      if (++arg < argc)
        expensiveSeed = strtouint32(argv[arg]);
    }

    else {
      err++;
    }

    arg++;
  }

  if ((doBoring    == false) &&
      (doInvert    == false) &&
      (doExpensive == false))
    err++;

  if (err) {
    fprintf(stderr, "usage: %s [-boring] [-expensive seed]\n", argv[0]);
    fprintf(stderr, "  -boring\n");
    fprintf(stderr, "  -invert\n");
    fprintf(stderr, "  -expensive [seed]\n");
  }


  if (doBoring)     boringTest();
  if (doInvert)     invertTest();
  if (doExpensive)  expensiveTest(expensiveSeed);

  exit(0);
}
