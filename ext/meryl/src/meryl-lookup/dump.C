
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


class dumpInput {
public:
  dumpInput() {
  };
  ~dumpInput() {
    delete [] fwd;
    delete [] rev;
  };

  dnaSeq        seq;
  uint64        seqIdx;

  kmvalu       *fwd    = nullptr;
  kmvalu       *rev    = nullptr;

  uint64        maxP;
};



static    //  (This really came from merfin)
void *
loadSequence(void *G) {
  lookupGlobal *g = (lookupGlobal *)G;
  dumpInput    *s = new dumpInput;

  if (g->seqFile1->loadSequence(s->seq) == false) {
    delete s;
    return(nullptr);
  }

  s->seqIdx = g->seqFile1->seqIdx();

  return(s);
}



static
void
processSequence(void *G, void *T, void *S) {
  lookupGlobal     *g   = (lookupGlobal *)G;
  dumpInput        *s   = (dumpInput    *)S;
  merylExactLookup *L = g->lookupDBs[0];

  //  Allocate and clear outputs.

  s->fwd = new kmvalu [s->seq.length()];
  s->rev = new kmvalu [s->seq.length()];

  for (uint32 ii=0; ii<s->seq.length(); ii++)
    s->fwd[ii] = s->rev[ii] = 0;

  //  Zip down all the kmers, saving the value of each.

  kmerIterator  kiter(s->seq.bases(), s->seq.length());

  while (kiter.nextMer()) {
    uint64  p = kiter.bgnPosition();

    s->fwd[p] = L->value(kiter.fmer());
    s->rev[p] = L->value(kiter.rmer());

    s->maxP = p+1;
  }

  //  Release the memory use for storing the sequence.

  s->seq.releaseBases();
}



static
void
outputSequence(void *G, void *S) {
  lookupGlobal     *g      = (lookupGlobal *)G;
  dumpInput        *s      = (dumpInput    *)S;

  //  Allocate space for the output string.

  resizeArray(g->outstring, 0, g->outstringMax, strlen(s->seq.ident()) + 16 + 16 + 16, _raAct::doNothing);

  //  Copy the sequence ident into the output strig.

  char *outptr = g->outstring;

  for (char const *x = s->seq.ident(); *x; )
    *outptr++ = *x++;

  *outptr++ = '\t';

  //  'outptr' is now where we start adding new info for each kmer,
  //  and we output the string from 'outroot'.

  for (uint64 p=0; p<s->maxP; p++) {
    char *t;

    if (s->fwd[p] + s->rev[p] == 0)
      continue;

    t = toDec(s->seqIdx, outptr);   *t++ = '\t';
    t = toDec(p, t);                *t++ = '\t';
    t = toDec(s->fwd[p], t);        *t++ = '\t';
    t = toDec(s->rev[p], t);        *t++ = '\n';   *t = 0;

    fputs(g->outstring, g->outFile1->file());
  }

  delete s;
}



void
dumpExistence(lookupGlobal *g) {
  sweatShop     *ss = new sweatShop(loadSequence, processSequence, outputSequence);

  ss->setLoaderQueueSize(4096);
  ss->setNumberOfWorkers(omp_get_max_threads());
  ss->setWriterQueueSize(4096);

  ss->run(g, g->showProgress);

  delete ss;
}
