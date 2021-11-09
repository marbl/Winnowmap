
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

#include "kmers.H"
#include "sequence.H"
#include "bits.H"



void
loadLookup(char const         *inputDBname,
           uint64              minV,
           uint64              maxV,
           merylExactLookup   &lookup) {

  fprintf(stderr, "==\n");
  fprintf(stderr, "==  Create merylExactLookup from '%s'.\n", inputDBname);
  fprintf(stderr, "==\n");

  merylFileReader   *merylDB = new merylFileReader(inputDBname);

  lookup.load(merylDB, 16.0, 0, minV, maxV);

  fprintf(stderr, "\n");

  delete merylDB;
}


void
loadMap(char const             *inputDBname,
        uint64                  minV,
        uint64                  maxV,
        std::map<kmer, kmvalu> &lookup) {

  fprintf(stderr, "==\n");
  fprintf(stderr, "==  Create map<kmer,kmvalu> from '%s'.\n", inputDBname);

  merylFileReader   *merylDB = new merylFileReader(inputDBname);

  uint64     nKmers = 0;
  uint64     nSkips = 0;

  while (merylDB->nextMer() == true) {
    kmer    kmer  = merylDB->theFMer();
    uint32  value = merylDB->theValue();

    if ((minV <= value) &&
        (value <= maxV)) {
      lookup[kmer] = value;
      nKmers++;
    } else {
      nSkips++;
    }

    if (((nKmers + nSkips) % 100000) == 0)
      fprintf(stderr, "==    Loaded %lu kmers; ignored %lu.\r", nKmers, nSkips);
  }

  fprintf(stderr, "==    Loaded %lu kmers; ignored %lu; map size %lu.\n", nKmers, nSkips, lookup.size());
  fprintf(stderr, "==\n");
  fprintf(stderr, "\n");

  delete merylDB;
}



int
main(int argc, char **argv) {
  char   *inputSeqName   = nullptr;
  char   *inputDBname    = nullptr;
  uint64  minV           = 0;
  uint64  maxV           = uint64max;
  uint32  threads        = 1;

  bool    skipCountCheck = false;

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  int                        arg = 1;
  while (arg < argc) {
    if      (strcmp(argv[arg], "-sequence") == 0)
      inputSeqName = argv[++arg];

    else if (strcmp(argv[arg], "-mers") == 0)
      inputDBname = argv[++arg];

    else if (strcmp(argv[arg], "-no-count") == 0)
      skipCountCheck = true;

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (inputSeqName == nullptr)   err.push_back("No input sequences (-sequence) supplied.\n");
  if (inputDBname  == nullptr)   err.push_back("No query meryl database (-mers) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s -sequence X.fasta -mers X.meryl\n", argv[0]);
    fprintf(stderr, "  -no-count      don't check the value of kmers\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Loads kmers in X.meryl into a merylExactLookup table and a standard\n");
    fprintf(stderr, "C++ associative map.  Verifies that every kmer present in X.fasta is\n");
    fprintf(stderr, "present in both the merylExactLookup and the associative map, and that\n");
    fprintf(stderr, "the value returned by both is the same.\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  merylExactLookup        kmerLookup;
  std::map<kmer,kmvalu>   kmerValue;
  std::map<kmer,kmvalu>   kmerCheck;

  loadLookup(inputDBname, minV, maxV, kmerLookup);
  loadMap   (inputDBname, minV, maxV, kmerValue);

  //fprintf(stderr, "==\n");
  //fprintf(stderr, "==  Copy kmerValue to kmerCheck.\n");
  //fprintf(stderr, "==\n");

  kmerCheck = kmerValue;

  //

  fprintf(stderr, "\n");
  fprintf(stderr, "==\n");
  fprintf(stderr, "==  Stream kmers from '%s'.\n", inputSeqName);
  fprintf(stderr, "==\n");

  dnaSeqFile  *seqFile    = new dnaSeqFile(inputSeqName);

  dnaSeq   seq;
  char     fString[64];
  char     rString[64];
  uint64   nTest = 0;
  uint64   nFail = 0;

  while (seqFile->loadSequence(seq)) {
    kmerIterator  kiter(seq.bases(), seq.length());

    while (kiter.nextMer()) {
      bool     fail  = false;
      kmer     fMer  = kiter.fmer();
      kmer     rMer  = kiter.rmer();
      kmer     cMer  = (fMer < rMer) ? fMer : rMer;
      kmvalu   value = 0;

      //
      //  Test exists() true/false.
      //

      if ((fail == false) && (kmerLookup.exists(cMer) == false)) {
        fprintf(stdout, "%s\t%s\t%s MISSING from kmerLookup::exists()\n",
                seq.ident(),
                kiter.fmer().toString(fString),
                kiter.rmer().toString(rString));
        fail = true;
      }

      //
      //  Test value().
      //

      if ((fail == false) && (kmerLookup.value(cMer) != kmerValue[cMer])) {
        fprintf(stdout, "%s\t%s\t%s MISSING from kmerLookup::value() -- kmerLookup=%u kmerValue=%u\n",
                seq.ident(),
                kiter.fmer().toString(fString),
                kiter.rmer().toString(rString),
                kmerLookup.value(cMer),
                kmerCheck[cMer]);
        fail = true;
      }

      //
      //  Test exists() true/false AND return the value.
      //

      if ((fail == false) && (kmerLookup.exists(cMer, value) == false)) {
        fprintf(stdout, "%s\t%s\t%s MISSING from kmerLookup::exists(mer, value) - (not found)\n",
                seq.ident(),
                kiter.fmer().toString(fString),
                kiter.rmer().toString(rString));
        fail = true;
      }
      if ((fail == false) && (value != kmerValue[cMer])) {
        fprintf(stdout, "%s\t%s\t%s MISSING from kmerLookup::exists(mer, value) -  kmerLookup=%u kmerValue=%u\n",
                seq.ident(),
                kiter.fmer().toString(fString),
                kiter.rmer().toString(rString),
                kmerLookup.value(cMer),
                kmerCheck[cMer]);
        fail = true;
      }

      //
      //  Subtract one from the kmer check counters.  If this is zero, the
      //  kmerIterator returned too many kmers.
      //
      if (skipCountCheck == false) {
        if (kmerCheck[cMer] == 0) {
          fprintf(stdout, "%s\t%s\t%s ZERO\n",
                  seq.ident(),
                  kiter.fmer().toString(fString),
                  kiter.rmer().toString(rString));
          fail = true;
        }

        --kmerCheck[cMer];
      }

      //  Log.

      if (fail)
        kmerLookup.exists_test(cMer);

      if (fail)
        nFail++;
#if 1
      else
        fprintf(stdout, "%s\t%s\t%s PASSES\n",
                seq.ident(),
                kiter.fmer().toString(fString),
                kiter.rmer().toString(rString));
#endif

      if ((++nTest % 100000) == 0)
        fprintf(stderr, "==    Tested %lu kmers.\r", nTest);
    }
  }

  delete seqFile;

  //  Check that all values are zero.

  if (skipCountCheck == false) {
    fprintf(stderr, "\n");
    fprintf(stderr, "==\n");
    fprintf(stderr, "==  Checking all kmers were seen.\n");
    fprintf(stderr, "==\n");

    for (auto it=kmerCheck.begin(); it != kmerCheck.end(); it++) {
      kmer    k = it->first;
      uint32  v = it->second;

      if (v != 0) {
        char   kmerString[64];

        fprintf(stderr, "%s\t%u\n", k.toString(kmerString), v);
      }
    }
  }

  if (nFail == 0) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Success!\n");
    fprintf(stderr, "\n");
    return(0);
  }

  return(1);
}
