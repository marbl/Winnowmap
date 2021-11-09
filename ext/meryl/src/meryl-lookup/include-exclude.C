
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



class filterInput {
public:
  filterInput() {
  };
  ~filterInput() {
  };

  dnaSeq        seq1;
  dnaSeq        seq2;

  uint64        nTotal;
  uint64        nFound;
};




static
void *
loadSequence(void *G) {
  lookupGlobal *g = (lookupGlobal *)G;
  filterInput  *s = new filterInput;

  bool   load1 = (g->seqFile1 != nullptr) && (g->seqFile1->loadSequence(s->seq1) == true);
  bool   load2 = (g->seqFile2 != nullptr) && (g->seqFile2->loadSequence(s->seq2) == true);

  if ((load1 == false) &&
      (load2 == false)) {
    delete s;
    return(nullptr);
  }

  return(s);
}



static
uint64
processSequence(merylExactLookup *L, dnaSeq &seq, bool is10x) {
  kmerIterator kiter(seq.bases(), seq.length());
  uint64       found = 0;

  //  Ignore the first 23 kmers in seq

  while (kiter.nextMer()) {
    if (is10x && kiter.bgnPosition() < 23)  continue;
    if ((L->value(kiter.fmer()) > 0) ||
        (L->value(kiter.rmer()) > 0))
      found++;
  }

  return(found);
}


static
void
processSequence(void *G, void *T, void *S) {
  lookupGlobal *g = (lookupGlobal *)G;
  filterInput  *s = (filterInput  *)S;

  //  Count the number of kmers found in the database from either
  //  seq1 or seq2.
  //
  //  If this is 10X Genomics reads, ignore counting in the first 23 bp of the seq1.

  s->nFound  = processSequence(g->lookupDBs[0], s->seq1, g->is10x);
  s->nFound += processSequence(g->lookupDBs[0], s->seq2, false);
}



static
void
outputSequence(compressedFileWriter  *O,
               dnaSeq                &seq,
               uint64                 nFound) {

  if (O == nullptr)
    return;

  if (seq.quals()[0] == 0)   outputFASTA(O->file(), seq.bases(),              seq.length(), 0, "%s nKmers=%lu", seq.ident(), nFound);
  else                       outputFASTQ(O->file(), seq.bases(), seq.quals(), seq.length(),    "%s nKmers=%lu", seq.ident(), nFound);
}


static
void
outputSequence(void *G, void *S) {
  lookupGlobal *g = (lookupGlobal *)G;
  filterInput  *s = (filterInput  *)S;

  g->nReadsTotal++;

  //  Write output if:
  //    'include' and    mers found.
  //    'exclude' and no mers found.

  if (((s->nFound  > 0) && (g->reportType == lookupOp::opInclude)) ||
      ((s->nFound == 0) && (g->reportType == lookupOp::opExclude))) {
    g->nReadsFound++;

    outputSequence(g->outFile1, s->seq1, s->nFound);
    outputSequence(g->outFile2, s->seq2, s->nFound);
  }

  delete s;
}



void
filter(lookupGlobal *g) {

  if (g->is10x)
    fprintf(stderr, "\nRunning in 10x mode. The first 23 bp of every sequence in %s will be ignored while looking up.\n", g->seqFile1->filename());
  
  sweatShop     *ss = new sweatShop(loadSequence, processSequence, outputSequence);

  ss->setLoaderQueueSize(10240);
  ss->setNumberOfWorkers(omp_get_max_threads());
  ss->setWriterQueueSize(omp_get_max_threads());

  ss->run(g, g->showProgress);

  delete ss;

  fprintf(stderr, "\nIncluding %lu reads (or read pairs) out of %lu.\n", g->nReadsTotal, g->nReadsFound);
}

