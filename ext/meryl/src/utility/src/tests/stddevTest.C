
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

#include "stddev.H"
#include "mt19937ar.H"

//  g++ -Wall -o stddevTest -I. -I.. stddevTest.C

void
testInsert(void) {
  stdDev<uint32>   sdu;
  stdDev<int32>    sdi;
  stdDev<double>   sdd;

  sdu.insert((uint32)2);
  sdu.insert((uint32)4);
  sdu.insert((uint32)4);
  sdu.insert((uint32)4);
  sdu.insert((uint32)5);
  sdu.insert((uint32)5);
  sdu.insert((uint32)7);
  sdu.insert((uint32)9);

  sdi.insert((uint32)2);
  sdi.insert((uint32)4);
  sdi.insert((uint32)4);
  sdi.insert((uint32)4);
  sdi.insert((uint32)5);
  sdi.insert((uint32)5);
  sdi.insert((uint32)7);
  sdi.insert((uint32)9);

  sdd.insert((uint32)2);
  sdd.insert((uint32)4);
  sdd.insert((uint32)4);
  sdd.insert((uint32)4);
  sdd.insert((uint32)5);
  sdd.insert((uint32)5);
  sdd.insert((uint32)7);
  sdd.insert((uint32)9);

  fprintf(stderr, "Expect mean=5, variance=%f, stddev=%f\n", 32.0 / 7.0, sqrt(32.0 / 7.0));

  fprintf(stderr, "  uint32  size %u mean %f variance %f stddev %f\n",
          sdu.size(), sdu.mean(), sdu.variance(), sdu.stddev());
  fprintf(stderr, "  int32   size %u mean %f variance %f stddev %f\n",
          sdi.size(), sdi.mean(), sdi.variance(), sdi.stddev());
  fprintf(stderr, "  double  size %u mean %f variance %f stddev %f\n",
          sdd.size(), sdd.mean(), sdd.variance(), sdd.stddev());

  assert(sdu.variance() == 32.0 / 7.0);
  assert(sdi.variance() == 32.0 / 7.0);
  assert(sdd.variance() == 32.0 / 7.0);

  fprintf(stderr, "\n\n");
}



void
testRemove(void) {
  double   values[10] = { 1, 2, 3, 4, 9, 8, 7, 6, 20, 30 };

  stdDev<double>  sd;

  fprintf(stderr, "Expect final to be zero, and insert() == remove().\n");

  for (int ii=0; ii<10; ii++) {
    sd.insert(values[ii]);
    fprintf(stderr, "insert[%2d] mean %8.4f stddev %8.4f\n", ii+1, sd.mean(), sd.stddev());
  }

  assert(sd.mean() == 9.0);

  fprintf(stderr, "\n");

  for (int ii=9; ii>=0; ii--) {
    sd.remove(values[ii]);
    fprintf(stderr, "remove[%2d] mean %8.4f stddev %8.4f\n", ii, sd.mean(), sd.stddev());
  }

  assert(sd.mean()   == 0.0);
  assert(sd.stddev() == 0.0);
}



void
testBig(uint32 nSamples) {
  histogramStatistics   hist;
  mtRandom              mt(10 + nSamples);

  fprintf(stderr, "\n");
  fprintf(stderr, "testBig for nSamples %u\n", nSamples);

  for (uint32 ii=0; ii<nSamples; ii++) {
    uint32  val = (mt.mtRandom32() % 10);
    uint32  num = (mt.mtRandom32() % 1000);

    //fprintf(stderr, "INSERT val %u num %u\n", val, num);

    hist.add(val, num);
  }

  hist.finalizeData();

  fprintf(stderr, "  size:   %lu\n", hist.numberOfObjects());
  fprintf(stderr, "  mean:   %f +- %f\n", hist.mean(), hist.stddev());
  fprintf(stderr, "  median: %lu +- %lu\n", hist.median(), hist.mad());
}



//  This is testing for an odd bit of apparent numerical instability where
//  the second to last remove() resulted in a negative variance.  d=0.0019
//  was the original case, but many others failed too.
void
testStability(void) {
  double sum = 0.0;

  fprintf(stderr, "\n");
  fprintf(stderr, "testStability (shouldn't crash)\n");

  for (double d = 0.0000; d < 0.5000; d += 0.0001) {
    stdDev<double>  sd;

    sd.insert(0.000000);
    sd.insert(d);
    sd.insert(0.000000);

    sd.remove(0.000000);
    sd.remove(0.000000);  //  Fails here; the two if's in remove() resolve.

    sum += sd.mean();     //  Add d.

    assert(d - 0.00001 <= sd.mean());
    assert(sd.mean() <= d + 0.00001);
    assert(sd.variance() == 0.0);

    sd.remove(d);

    sum += sd.mean();     //  Add zero.

    assert(sd.mean()     == 0.0);
    assert(sd.variance() == 0.0);
  }

  fprintf(stderr, "  %18.16f\n", sum);
}



//  Same idea, but this one fails before we hit the
//  reset for one item.  Grrrr!
void
testStability2(uint32 n) {
  double sum = 0.0;
  stdDev<double>  sd;

  if (n == 1) {
    fprintf(stderr, "\n");
    fprintf(stderr, "testStability2 (values should be positive zero)\n");
  }

  for (uint32 ii=0; ii<n; ii++)
    sd.insert(0.000000);

  sd.insert(0.000190);
  sd.remove(0.000190);
  fprintf(stderr, "%2u  %26.24f\n", n, sd.variance());
  assert(sd.variance() >= 0.0);

  sd.insert(0.000220);
  sd.remove(0.000220);
  fprintf(stderr, "%2u  %26.24f\n", n, sd.variance());
  assert(sd.variance() >= 0.0);

  for (uint32 ii=0; ii<n; ii++)
    sd.remove(0.000000);
}



int
main(int argc, char **argv) {

  testInsert();
  testRemove();

  testBig(1);
  testBig(2);
  testBig(3);
  testBig(4);
  testBig(5);
  testBig(10);
  testBig(100);
  testBig(1000);

  testStability();

  testStability2(1);
  testStability2(2);
  testStability2(3);
  testStability2(4);

  fprintf(stderr, "\n");
  fprintf(stderr, "Success!\n");

  exit(0);
}
