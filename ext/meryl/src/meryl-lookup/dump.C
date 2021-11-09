
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
    for (uint32 ii=0; ii<existLen; ii++)
      delete exist[ii];
    delete [] exist;
    delete [] count;
    delete [] depth;
  };

  void     outputBED    (lookupGlobal *g);
  void     outputBEDruns(lookupGlobal *g);
  void     outputWIG    (lookupGlobal *g);

  dnaSeq        seq;          //  Sequence we're processing.
  uint64        seqIdx = 0;   //  Index of that sequence in the input file.

  uint64        maxP   = 0;   //  Maximum position set in the following data.

  //  BED format output is:
  //    seqName bgn end [label [score [...]]]
  //
  //  This will be stored as a bitvector, one vector per input DB.

  uint32        existLen = 0;
  bitArray    **exist    = nullptr;

  //  WIGGLE format output is:
  //    track type=<name>             <wiggle track definition line>
  //    variableStep chrom=<name>
  //    position value
  //    position value
  //
  //  'value' is either:
  //    the sum of the fwd and rev kmer counts,
  //    the number of kmers touching this base (aka, depth)

  kmvalu       *count = nullptr;
  uint8        *depth = nullptr;
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

  //
  //  In all the lookups below, we ask for both the F and the R mer, instead
  //  of just the canonical mer, so we can support non-canonical databases.
  //

  kmerIterator      kiter(s->seq.bases(), s->seq.length());

  //  For BED, remember which kmers are found in each database.

  if (g->reportType == lookupOp::opBED) {
    s->existLen = g->lookupDBs.size();
    s->exist    = new bitArray * [s->existLen];

    for (uint32 dd=0; dd<s->existLen; dd++)
      s->exist[dd] = new bitArray(s->seq.length());

    while (kiter.nextMer()) {
      uint64  p  = kiter.bgnPosition();
      kmer    f  = kiter.fmer(), c;
      kmer    r  = kiter.rmer(), n;

      if (f < r) {   //  A small optimization; since (hopefully)
        c = f;       //  the usual use case will be with a
        n = r;       //  canonical DB, we can test for the
      } else {       //  canonical and non-canonical kmers,
        c = r;       //  instead of the f and r kmers.  If we
        n = f;       //  test the canonical first, we'll skip the
      }              //  non-canonical lookup.

      for (uint32 dd=0; dd<g->lookupDBs.size(); dd++) {
        merylExactLookup *L = g->lookupDBs[dd];

        if ((L->exists(c) == true) ||
            (L->exists(n) == true)) {
          if (g->lookupDBlabelLen > 0)       //  If labels are present, remember which
            s->exist[dd]->setBit(p, true);   //  database the kmer was found in.  If not,
          else                               //  remove duplicate outputs by flagging the
            s->exist[0]->setBit(p, true);    //  kmer as present only in the first db.

          s->maxP = p+1;
        }
      }
    }
  }

  //  For WIGgle, remember the count of the kmer in all databases.

  if (g->reportType == lookupOp::opWIGcount) {
    s->count = new kmvalu [s->seq.length()];

    for (uint32 ii=0; ii<s->seq.length(); ii++)
      s->count[ii] = 0;

    while (kiter.nextMer()) {
      uint64  p  = kiter.bgnPosition();
      kmer    f  = kiter.fmer();
      kmer    r  = kiter.rmer();

      for (uint32 dd=0; dd<g->lookupDBs.size(); dd++) {
        merylExactLookup *L = g->lookupDBs[dd];

        kmvalu  fv = L->value(f);   //  Ask for each individually, instead of the
        kmvalu  rv = L->value(r);   //  canonical, to support non-canonical DBs.

        if (f == r)                 //  Don't double count palindromes.
          s->count[p] += fv;
        else
          s->count[p] += fv + rv;

        s->maxP = p+1;
      }
    }
  }

  //  For WIGgle, compute the depth of coverage of kmers that exist in the DB.
  //
  //  BPW expects that the obvious implementation used below will be faster
  //  than the clever implementation (build list of increment, decrement,
  //  sort, scan list to find depth at each position) but didn't test.
  //
  //  It will certainly be smaller memory.

  if (g->reportType == lookupOp::opWIGdepth) {
    merylExactLookup *L = g->lookupDBs[0];

    s->depth = new uint8 [s->seq.length()];

    for (uint32 ii=0; ii<s->seq.length(); ii++)
      s->depth[ii] = 0;

    //  Confused.  The SIMPLE_DEPTH method _should_ be slower; it does
    //  'kmersize' work per kmer, compared to N work.  But if the number of
    //  kmers found is small, it turns out to be (slightly) faster.
    //
    //  On a large database of 'the wrong' kmers (homopoly compressed chm13 reads)
    //  a human genome assembly (hg0733) takes:
    //    1080 user seconds, 7:25 wall  SIMPLE_DEPTH
    //    1123 user seconds, 7:33 wall
    //
    //  But on a large database of 'the correct' kmers (chm13 reads, uncompressed)
    //  a human genome assembly (chm13) takes:
    //    3699 user seconds, 30:29 wall  SIMPLE_DEPTH
    //    3740 user seconds, 27:45 wall
    //  A second run:
    //    3608 user seconds, 27:56 wall  SIMPLE_DEPTH
    //    3622 user seconds, 26:49 wall
    //  So largely not different.  I'm leaving SIMPLE_DEPTH _disabled_;
    //  suspecting that a larger k (k=22 here) will matter.

#define SIMPLE_DEPTH
#undef  SIMPLE_DEPTH

#ifdef SIMPLE_DEPTH

    while (kiter.nextMer()) {
      kmer    f = kiter.fmer();
      kmer    r = kiter.rmer();

      if ((L->exists(f) == true) ||
          (L->exists(r) == true)) {
        for (uint64 p=kiter.bgnPosition(); p<kiter.endPosition(); p++)
          s->depth[p]++;

        s->maxP = kiter.endPosition();
      }
    }

#else

    while (kiter.nextMer()) {
      kmer    f = kiter.fmer();
      kmer    r = kiter.rmer();

      if ((L->exists(f) == true) ||
          (L->exists(r) == true)) {
        s->depth[kiter.bgnPosition()] += 1;
        s->depth[kiter.endPosition()] -= 1;

        s->maxP = kiter.endPosition();
      }
    }
  
    //  Convert that change-in-depth into depth, and rewrite the array.

    uint32 d = 0;
    for (uint64 pp=0; pp<=s->maxP; pp++) {
      d += s->depth[pp];
      s->depth[pp] = d;
    }

#endif
  }


  //  Release the memory use for storing the sequence.

  s->seq.releaseBases();
}




void
dumpInput::outputBED(lookupGlobal *g) {

  //  Allocate space for the output string.

  uint32  maxIdent = strlen(seq.ident());
  uint32  maxLabel = g->lookupDBlabelLen;

  resizeArray(g->outstring, 0, g->outstringMax, maxIdent + 16 + 16 + maxLabel + 1, _raAct::doNothing);

  //  Copy the sequence ident into the output string.

  char *outptr = g->outstring;

  for (char const *x = seq.ident(); *x; )
    *outptr++ = *x++;

  *outptr++ = '\t';

  //  'outptr' is now where we start adding new info for each kmer,
  //  and we output the string from 'outroot'.

  uint32  k = kmer::merSize();

  for (uint64 p=0; p<maxP; p++) {
    for (uint32 dd=0; dd<g->lookupDBs.size(); dd++) {
      char *t;

      if (exist[dd]->getBit(p) == false)                   //  Not set?  No kmer here.
        continue;

      t = toDec(p,     outptr);                            //  Append the begin position.

      *t++ = '\t';
      t = toDec(p + k, t);                                 //  Append the end position.

      if (dd < g->lookupDBlabel.size()) {                  //  If a label exists,
        *t++ = '\t';                                       //
        for (char const *x = g->lookupDBlabel[dd]; *x; )   //  append the label.
          *t++ = *x++;
      }

      *t++ = '\n';                                         //  Terminate the string.
      *t   = 0;

      fputs(g->outstring, g->outFile1->file());
    }
  }
}


void
dumpInput::outputBEDruns(lookupGlobal *g) {

  //  Allocate space for the output string.

  uint32  maxIdent = strlen(seq.ident());
  uint32  maxLabel = g->lookupDBlabelLen;

  resizeArray(g->outstring, 0, g->outstringMax, maxIdent + 16 + 16 + maxLabel + 1, _raAct::doNothing);

  //  Copy the sequence ident into the output string.

  char *outptr = g->outstring;

  for (char const *x = seq.ident(); *x; )
    *outptr++ = *x++;

  *outptr++ = '\t';

  //  'outptr' is now where we start adding new info for each kmer,
  //  and we output the string from 'outroot'.

  uint32  k = kmer::merSize();

  uint64  bgn[g->lookupDBs.size()];
  for (uint32 ii=0; ii < g->lookupDBs.size(); ii++)
    bgn[ii] = uint32max;

  for (uint64 p=0; p <= maxP; p++) {
    for (uint32 dd=0; dd<g->lookupDBs.size(); dd++) {
      bool   bit = (p < maxP) ? (exist[dd]->getBit(p)) : false;
      char  *t;

      if (bit == true) {            //  If a true bit, and we aren't in a run
        if (bgn[dd] == uint32max)   //  already, remember the start of the run.
          bgn[dd] = p;              //
        continue;                   //  Then keep on trucking.
      }

      if (bgn[dd] == uint32max)     //  And if we aren't in a run,
        continue;                   //  keep on searching.

      //  Otherwise, we've just found the end of a run, so output it.

      t = toDec(bgn[dd], outptr);                          //  Append the begin position.

      bgn[dd] = uint32max;                                 //  Declare that we aren't in a run.

      *t++ = '\t';
      t = toDec(p + k, t);                                 //  Append the end position.

      if (dd < g->lookupDBlabel.size()) {                  //  If a label exists,
        *t++ = '\t';                                       //
        for (char const *x = g->lookupDBlabel[dd]; *x; )   //  append the label.
          *t++ = *x++;
      }

      *t++ = '\n';                                         //  Terminate the string.
      *t   = 0;

      fputs(g->outstring, g->outFile1->file());
    }
  }
}


void
dumpInput::outputWIG(lookupGlobal *g) {

  //  Allocate space for the output string - this is fixed length, just an
  //  integer position and an integer count/depth.

  resizeArray(g->outstring, 0, g->outstringMax, 16 + 16 + 1, _raAct::doNothing);

  //  Output a track definition line.

  fprintf(g->outFile1->file(), "variableStep chrom=%s\n", seq.ident());

  //  Then all the data.

  bool   isCount = (g->reportType == lookupOp::opWIGcount);
  bool   isDepth = (g->reportType == lookupOp::opWIGdepth);

  for (uint64 p=0; p<maxP; p++) {
    char *t;

    if (((isCount) && (count[p] == 0)) ||
        ((isDepth) && (depth[p] == 0)))
      continue;

    t = toDec(p + 1, g->outstring);

    *t++ = '\t';

    if (isCount)
      t = toDec(count[p], t);

    if (isDepth)
      t = toDec(depth[p], t);

    *t++ = '\n';
    *t   = 0;

    fputs(g->outstring, g->outFile1->file());
  }
}


static
void
outputSequence(void *G, void *S) {
  lookupGlobal     *g      = (lookupGlobal *)G;
  dumpInput        *s      = (dumpInput    *)S;

  if ((g->reportType == lookupOp::opBED) && (g->mergeBedRuns == false))
    s->outputBED(g);

  if ((g->reportType == lookupOp::opBED) && (g->mergeBedRuns == true))
    s->outputBEDruns(g);

  if (g->reportType == lookupOp::opWIGcount)
    s->outputWIG(g);

  if (g->reportType == lookupOp::opWIGdepth)
    s->outputWIG(g);

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
