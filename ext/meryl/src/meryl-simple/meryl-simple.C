
/******************************************************************************
 *
 *  This file is part of meryl, a genomic k-kmer counter with nice features.
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
#include "sequence.H"
#include "kmers.H"
//#include "strings.H"
//#include "system.H"



int
main(int argc, char **argv) {
  uint32  kSize    = 0;
  char   *sInput   = NULL;   //  input sequences
  char   *mOutput  = NULL;   //  meryl-format output database
  char   *dOutput  = NULL;   //  dump output
  char   *hOutput  = NULL;   //  histogram output
  uint64  memLimit = UINT64_MAX;
  uint64  memUsed  = 0;

  argc = AS_configure(argc, argv, 1);

  std::vector<char *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-k") == 0) {
      kSize = strtouint32(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-S") == 0) {
      sInput = argv[++arg];
    }

    else if (strcmp(argv[arg], "-M") == 0) {
      mOutput = argv[++arg];
    }

    else if (strcmp(argv[arg], "-D") == 0) {
      dOutput = argv[++arg];
    }

    else if (strcmp(argv[arg], "-H") == 0) {
      hOutput = argv[++arg];
    }

    else if (strcmp(argv[arg], "-m") == 0) {
      memLimit   = strtouint32(argv[++arg]);
      memLimit <<= 20;
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  //  If any errors, fail.

  if ((argc == 1) ||        //  No commands
      (err.size() > 0)) {   //  Errors
    fprintf(stderr, "usage: %s -k kmerSize -S input.fasta ...\n", argv[0]);
    fprintf(stderr, "  -k kmerSize\n");
    fprintf(stderr, "  -S input.fasta\n");
    fprintf(stderr, "  -M output.meryl\n");
    fprintf(stderr, "  -D output.dump\n");
    fprintf(stderr, "  -H output.histogram\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -m memLimit_in_MB\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    exit(1);
  }

  kmer::setSize(kSize);

  //  We need to know how many kmers are in the input so we can
  //  allocate a gigantic array.  We could just reallocate and copy,
  //  but then we're limited to 1/2 the memory.

  dnaSeqFile   *seqFile  = new dnaSeqFile(sInput, false);
  dnaSeq       *seq      = new dnaSeq;
  uint64        kmersLen = 0;
  uint64        kmersMax = 0;
  kmdata       *kmers    = NULL;

  fprintf(stderr, "-- Scanning input sequences to determine the number of kmers.\n");

  while (seqFile->loadSequence(*seq) == true)
    kmersMax += seq->length();

  delete seq;
  delete seqFile;

  //  Allocate space for the kmers list, reopen the file, build the list.

  memUsed  = sizeof(kmdata) * kmersMax;

  if (memUsed > memLimit) {
    fprintf(stderr, "-- Need %lu MB for kmers data, limited to %lu MB.  Fail.\n",
            memUsed >> 20, memLimit >> 20);
    exit(1);
  }

  fprintf(stderr, "-- Allocating space for %lu kmers (%lu MB).\n", kmersMax, memUsed >> 20);

  kmers   = new kmdata [kmersMax];
  seqFile = new dnaSeqFile(sInput, false);
  seq     = new dnaSeq;

  fprintf(stderr, "-- Loading kmers.\n");

  while (seqFile->loadSequence(*seq) == true) {
    kmerIterator  ki(seq->bases(), seq->length());

    //fprintf(stderr, "--   '%s'\n", seq->name());

    while (ki.nextMer())
      if (ki.fmer() < ki.rmer())
        kmers[kmersLen++] = ki.fmer();
      else
        kmers[kmersLen++] = ki.rmer();

    //fprintf(stderr, "--   '%s' %lu\n", seq->name(), kmersLen);

    assert(kmersLen < kmersMax);
  }

  delete seq;
  delete seqFile;

  //  Sort.

  fprintf(stderr, "-- Sorting %lu kmers.\n", kmersLen);

  std::sort(kmers, kmers + kmersLen);

  //  Scan, count and output stuff.

  uint32    histLen = 0;
  uint32    histMax = 16 * 1024 * 1024;
  uint32   *hist    = new uint32 [histMax];
  kmer      k;
  uint32    c;
  char      kstr[65];

  for (uint32 ii=0; ii<histMax; ii++)
    hist[ii] = 0;

  //LE *M = AS_UTL_openOutputFile(mOutput);
  FILE *D = AS_UTL_openOutputFile(dOutput);

  fprintf(stderr, "-- Output.\n");

  for (uint64 ii=0, jj=1; ii<kmersLen; ) {
    while ((kmers[jj] == kmers[ii]) && (jj < kmersLen))
      jj++;

    k.setPrefixSuffix(0, kmers[ii], kSize);
    c = jj - ii;

    if (mOutput) {
    }

    if (dOutput) {
      fprintf(D, "%s\t%u\n", k.toString(kstr), c);
    }

    if (c < histMax)
      hist[c]++;
    if (histLen < c)
      histLen = c;

    ii = jj;
    jj = jj + 1;
  }

  AS_UTL_closeFile(D);


  if (hOutput) {
    FILE *H = AS_UTL_openOutputFile(hOutput);

    for (uint32 ii=0; ii<=histLen; ii++)
      if (hist[ii] > 0)
        fprintf(H, "%u\t%u\n", ii, hist[ii]);

    AS_UTL_closeFile(H);
  }

  //  Cleanup.

  delete [] hist;
  delete [] kmers;

  fprintf(stderr, "Bye.\n");

  return(0);
}
