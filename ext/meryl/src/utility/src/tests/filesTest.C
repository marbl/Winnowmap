
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

//  g++6 -o filesTest -I.. -I. filesTest.C files.C

#include "files.H"

typedef  uint8   TYPE;


int32
main(int32 argc, char **argv) {
  uint64   nObj  = (uint64)16 * 1024 * 1024;
  TYPE    *array = new TYPE [nObj];
  TYPE     value = 0;

  if (1) {
    compressedFileWriter *out = new compressedFileWriter("test.gz");

    delete out;
  }

  if (1) {
    fprintf(stderr, "Initializing.\n");

    for (uint64 ii=0; ii<nObj; ii++)
      array[ii] = ii;

    fprintf(stderr, "Writing.\n");

    FILE *OUT = AS_UTL_openOutputFile("./filesTest.dat");

    writeToFile(array, "array", nObj, OUT);

    AS_UTL_closeFile(OUT);
  }


  if (1) {
    fprintf(stderr, "Reading - as one block.\n");

    FILE *IN = AS_UTL_openInputFile("./filesTest.dat");
    loadFromFile(array, "array", nObj, IN);
    AS_UTL_closeFile(IN);

    for (uint64 ii=0; ii<nObj; ii++)
      assert(array[ii] == (TYPE)ii);
  }


  if (1) {
    fprintf(stderr, "Reading.\n");

    FILE *IN = AS_UTL_openInputFile("./filesTest.dat");

    for (uint64 ii=0; ii<nObj; ii++) {
      loadFromFile(value, "value", IN);

      assert(value == (TYPE)ii);
    }

    fprintf(stderr, "Reading - one after eof.\n");
    loadFromFile(value, "value", IN, false);
    loadFromFile(value, "value", IN, true);

    AS_UTL_closeFile(IN);
  }


  exit(0);
}

