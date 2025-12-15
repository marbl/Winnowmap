
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

using namespace merylutil;


int32
main(int32 argc, char **argv) {
  uint32   Alen=0, Amax=0;   uint32  *A=nullptr;
  uint32   Blen=0, Bmax=0;   uint32  *B=nullptr;
  uint32   Clen=0, Cmax=0;   uint64  *C=nullptr;
  uint32   Dlen=0, Dmax=0;   uint64  *D=nullptr;

  fprintf(stderr, "Testing increaseArray() and related.\n");
  fprintf(stderr, "For best results, run in valgrind.\n");
  fprintf(stderr, "\n");

  //

  fprintf(stderr, "allocateArray()\n");

  allocateArray(A, Amax, 4000);       //  Allocate new arrays, both will
  allocateArray(B, Bmax = 4000);      //  set elements to zero by default.

  fprintf(stderr, "  setting\n");

  for (uint32 ii=0; ii<Amax; ii++)     A[ii] = ii << 8;
  for (uint32 ii=0; ii<Bmax; ii++)     B[ii] = ii << 16;

  fprintf(stderr, "  checking\n");

  for (uint32 ii=0; ii<Amax; ii++)     assert(A[ii] == ii << 8);
  for (uint32 ii=0; ii<Bmax; ii++)     assert(B[ii] == ii << 16);

  fprintf(stderr, "  releasing\n");

  delete [] A;  Amax = Alen = 0;  A = nullptr;
  delete [] B;  Bmax = Blen = 0;  B = nullptr;
  delete [] C;  Cmax = Clen = 0;  C = nullptr;
  delete [] D;  Dmax = Dlen = 0;  D = nullptr;

  //

  fprintf(stderr, "increaseArray() from nullptr\n");

  Alen = 3332;  increaseArray(A,    Alen, Amax, 3333, _raAct::copyDataClearNew);  assert(Amax == 3333);
  Alen = 3333;  increaseArray(A,    Alen, Amax, 3333, _raAct::copyDataClearNew);  assert(Amax == 6666);
  Blen = 3334;  increaseArray(B,    Blen, Bmax, 3333, _raAct::copyDataClearNew);  assert(Bmax == 6666);

  Clen = 3332;  increaseArray(C, D, Clen, Cmax, 3333, _raAct::copyDataClearNew);  assert(Cmax == 3333);
  Clen = 6667;  increaseArray(C, D, Clen, Cmax, 3333, _raAct::copyDataClearNew);  assert(Cmax == 9999);

  fprintf(stderr, "  setting\n");

  for (uint32 ii=0; ii<Alen; ii++)      A[ii] = ii << 8;
  for (uint32 ii=0; ii<Blen; ii++)      B[ii] = ii << 16;
  for (uint64 ii=0; ii<Clen; ii++)    { C[ii] = ii << 24;  D[ii] = ii << 32; }

  fprintf(stderr, "  checking\n");

  for (uint32 ii=0; ii<Alen; ii++)      assert(A[ii] == ii << 8);
  for (uint32 ii=0; ii<Blen; ii++)      assert(B[ii] == ii << 16);
  for (uint64 ii=0; ii<Clen; ii++)    { assert(C[ii] == ii << 24);
                                        assert(D[ii] == ii << 32); }

  for (uint32 ii=Alen; ii<Amax; ii++)   assert(A[ii] == 0);
  for (uint32 ii=Blen; ii<Bmax; ii++)   assert(B[ii] == 0);
  for (uint64 ii=Clen; ii<Cmax; ii++) { assert(C[ii] == 0);
                                        assert(D[ii] == 0); }

  fprintf(stderr, "  releasing\n");

  //lete [] A;  Amax = Alen = 0;  A = nullptr;
  delete [] B;  Bmax = Blen = 0;  B = nullptr;
  //lete [] C;  Cmax = Clen = 0;  C = nullptr;
  delete [] D;  Dmax = Dlen = 0;  D = nullptr;

  //

  fprintf(stderr, "duplicateArray() (lengths %u %u %u %u)\n", Alen, Blen, Clen, Dlen);

  duplicateArray(B, Blen, Bmax, A, Alen, Amax);        //  Duplicate the arrays,
  duplicateArray(D, Dlen, Dmax, C, Clen, Cmax);        //  the last call will force
  duplicateArray(D, Dlen, Dmax, C, Clen, Cmax, true);  //  a reallocation.

  assert(Blen == Alen);
  assert(Bmax == Alen);

  assert(Dlen == Clen);
  assert(Dmax == Clen);

  fprintf(stderr, "  checking\n");

  for (uint32 ii=0; ii<Blen; ii++)  assert(B[ii] == A[ii]);
  for (uint64 ii=0; ii<Dlen; ii++)  assert(D[ii] == C[ii]);

  //

  fprintf(stderr, "Testing increaseArray(single).\n");

  Alen *= 2;
  for (uint32 ii=0; ii<Alen; ii++) {
    increaseArray(A, ii, Amax, 3333, _raAct::copyDataClearNew);
    A[ii] = ii;
  }

  fprintf(stderr, "  Reallocated A from Alen=%u to Alen=%u, now with Amax %lu\n", Alen/2, Alen, Amax);

  for (uint32 ii=0; ii<Alen; ii++)
    assert(A[ii] == ii);
  for (uint32 ii=Alen; ii<Amax; ii++)
    assert(A[ii] == 0);

  fprintf(stderr, "  releasing\n");

  delete [] A;  Amax = Alen = 0;  A = nullptr;
  delete [] B;  Bmax = Blen = 0;  B = nullptr;
  delete [] C;  Cmax = Clen = 0;  C = nullptr;
  delete [] D;  Dmax = Dlen = 0;  D = nullptr;

  //

  fprintf(stderr, "Testing increaseArray(pair).\n");

  Clen *= 2;
  for (uint64 ii=0; ii<Clen; ii++) {
    increaseArrayPair(C, D, ii, Cmax, 3333, _raAct::copyDataClearNew);
    C[ii] = ii << 8;
    D[ii] = ii << 16;

    //fprintf(stderr, "ii = %u %u %u -- %lu %lu\n", ii, ii << 8, ii << 16, C[ii], D[ii]);

    //if (ii > 0) {
    //  fprintf(stderr, "     %u %u %u -- %lu %lu\n", ii-1, (ii-1) << 8, (ii-1) << 16, C[ii-1], D[ii-1]);
    //  assert(C[ii-1] == (ii-1) << 8);
    //  assert(D[ii-1] == (ii-1) << 16);
    //}
  }

  fprintf(stderr, "  Reallocated C and D from Clen=%u to Clen=%u, now with Cmax %lu\n", Clen/2, Clen, Cmax);

  for (uint64 ii=0; ii<Clen; ii++) {
    assert(C[ii] == ii << 8);
    assert(D[ii] == ii << 16);
  }

  fprintf(stderr, "  testing zero\n");

  for (uint64 ii=Clen; ii<Cmax; ii++) {
    assert(C[ii] == 0);
    assert(D[ii] == 0);
  }

  delete [] A;  Amax = Alen = 0;  A = nullptr;
  delete [] B;  Bmax = Blen = 0;  B = nullptr;
  delete [] C;  Cmax = Clen = 0;  C = nullptr;
  delete [] D;  Dmax = Dlen = 0;  D = nullptr;

  //

  fprintf(stderr, "Testing increaseArray(pair) zero.\n");

  increaseArrayPair(C, D, 0, Cmax, 3333, _raAct::copyDataClearNew);
  increaseArrayPair(C, D, 1, Cmax, 3333, _raAct::copyDataClearNew);
  increaseArrayPair(C, D, 0, Cmax, 3333, _raAct::copyDataClearNew);

  fprintf(stderr, "Testing increaseArray(pair) stepped.\n");

  for (Clen=10000; Clen < 110000; Clen += 10000) {
    fprintf(stderr, "  step from %lu to %lu\n", Clen-10000, Clen);

    increaseArrayPair(C, D, Clen, Cmax, 3333, _raAct::copyDataClearNew);

    fprintf(stderr, "  step from %lu to %lu cmax now %u\n", Clen-10000, Clen, Cmax);

    assert(Clen <= Cmax);

    fprintf(stderr, "  setting (%u .. %u)\n", Clen-10000, Clen);

    for (uint64 ii=Clen-10000; ii<Clen; ii++) {
      C[ii] = ii << 8;
      D[ii] = ii << 16;
    }

    fprintf(stderr, "  testing\n");

    for (uint64 ii=0; ii<Clen; ii++) {
      assert(C[ii] == ii << 8);
      assert(D[ii] == ii << 16);
    }

    fprintf(stderr, "  testing zero\n");

    for (uint64 ii=Clen; ii<Cmax; ii++) {
      assert(C[ii] == 0);
      assert(D[ii] == 0);
    }
  }

  delete [] A;  Amax = Alen = 0;  A = nullptr;
  delete [] B;  Bmax = Blen = 0;  B = nullptr;
  delete [] C;  Cmax = Clen = 0;  C = nullptr;
  delete [] D;  Dmax = Dlen = 0;  D = nullptr;

  return(0);
}
