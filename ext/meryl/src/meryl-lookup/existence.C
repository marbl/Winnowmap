
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

#include "meryl-lookup.H"
#include "sweatShop.H"



class existInput {
public:
  existInput() {
  };
  ~existInput() {
    delete [] nFound;
  };

  dnaSeq        seq;

  uint64        nTotal  = 0;
  uint64       *nFound  = nullptr;
};



static    //  (This really came from merfin)
void *
loadSequence(void *G) {
  lookupGlobal *g = (lookupGlobal *)G;
  existInput   *s = new existInput();

  if (g->seqFile1->loadSequence(s->seq) == false) {
    delete s;
    return(nullptr);
  }

  return(s);
}



static
void
processSequence(void *G, void *T, void *S) {
  lookupGlobal *g   = (lookupGlobal *)G;
  existInput   *s   = (existInput   *)S;
  int32         nIn = g->lookupDBs.size();

  //  Allocate and clear outputs.

  s->nTotal = 0;
  s->nFound = new uint64 [nIn];

  for (uint32 dd=0; dd<nIn; dd++)
    s->nFound[dd] = 0;

  //  Zip through the kmers, counting how many kmers we have and how many we
  //  found in each input.

  kmerIterator  kiter(s->seq.bases(), s->seq.length());

  while (kiter.nextMer()) {
    s->nTotal++;

    for (uint32 dd=0; dd<nIn; dd++) {
      if ((g->lookupDBs[dd]->value(kiter.fmer()) > 0) ||
          (g->lookupDBs[dd]->value(kiter.rmer()) > 0))
        s->nFound[dd]++;
    }
  }

  //  Release the memory use for storing the sequence.

  s->seq.releaseBases();
}



static
void
outputSequence(void *G, void *S) {
  lookupGlobal *g   = (lookupGlobal *)G;
  existInput   *s   = (existInput   *)S;
  int32         nIn = g->lookupDBs.size();

  //  Allocate space for the output string.

  resizeArray(g->outstring, 0, g->outstringMax, 16 + 16 * 2 * nIn, _raAct::doNothing);

  //  Create the string.

  char *t = g->outstring;

  *t++ = '\t';
  t = toDec(s->nTotal, t);

  for (uint32 dd=0; dd<nIn; dd++) {
    *t++ = '\t';
    t = toDec(g->lookupDBs[dd]->nKmers(), t);

    *t++ = '\t';
    t = toDec(s->nFound[dd], t);
  }

  *t++ = '\n';
  *t   = 0;

  //  And output it.

  fputs(s->seq.ident(), g->outFile1->file());
  fputs(g->outstring,   g->outFile1->file());

  delete s;
}




void
reportExistence(lookupGlobal *g) {
  sweatShop     *ss = new sweatShop(loadSequence, processSequence, outputSequence);

  ss->setLoaderQueueSize(4096);
  ss->setNumberOfWorkers(omp_get_max_threads());
  ss->setWriterQueueSize(4096);

  ss->run(g, g->showProgress);

  delete ss;
}
