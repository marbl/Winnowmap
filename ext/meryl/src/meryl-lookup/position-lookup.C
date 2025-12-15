
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

#include <algorithm>

#include "meryl-lookup.H"
#include "system.H"

using namespace merylutil;

class hit_s {
public:
  hit_s() {
  }
  hit_s(uint32 rid, uint32 rpos, uint32 qid, uint32 qpos) :
    refID(rid), refPos(rpos), qryID(qid), qryPos(qpos) {
  }

  uint32  refID  = 0;   //  refID is relative to the reference file.
  uint32  refPos = 0;   //  refPos is the index of the kmer in the lookup table.
  uint32  qryID  = 0;   //  qryID is relative to the batch.
  uint32  qryPos = 0;   //  qryPos is the position of the kmer in the query sequence.
};


class seqBatch {
public:
  seqBatch();
  ~seqBatch();

  uint64                 batchid = 0;

  std::vector<dnaSeq *>  sequences;
  std::vector<hit_s>     hits;
};

seqBatch::seqBatch() {
}
seqBatch::~seqBatch() {
  for (dnaSeq *seq : sequences)
    delete seq;
}


class posLookGlobal {
public:
  posLookGlobal();
  ~posLookGlobal();

  void      addInputFile (char const *in)  { inputNames.push_back(in); };
  void      setRefMerName(char const *in)  { refMerName = in; };
  void      setRefSeqName(char const *in)  { refSeqName = in; };

  void      initialize(void);

  seqBatch *loadBatch(uint64 nSeqs, uint64 nBases);
  void      lookupBatch(seqBatch *g);
  void      writeHits(uint64 idx, bool doWrite);
  void      writeBatch(seqBatch  *g);

  void      savePainting(char const *filename);

private:
  bool      openInputFile(void) {
    if (inputFile != nullptr)                 return(true);   //  we have an input file already
    if (inputNamesPos >= inputNames.size())   return(false);  //  no more input files

    inputFile = openSequenceFile(inputNames[inputNamesPos]);  //  open next file
    inputNamesPos++;                                          //  point to the next file name

    return(true);
  }

private:
  char const                *refMerName = nullptr;
  char const                *refSeqName = nullptr;

  std::vector<char const *>  inputNames;
  uint32                     inputNamesPos = 0;

  dnaSeqFile                *inputFile  = nullptr;

  uint64                     batchid = 0;

  merylExactLookup          *refLookup  = nullptr;
  uint32                    *nQmerPer   = nullptr;   //  Number of query kmers with a hit here
  uint32                    *nQseqPer   = nullptr;   //  Number of query sequence with a hit here

public:
  compressedFileWriter      *hitsPerQuery      = nullptr;
  compressedFileWriter      *paintMerPerBase   = nullptr;
  compressedFileWriter      *paintQueryPerBase = nullptr;

  void                       saveMerPerBase(uint32 windowSize);
  void                       saveQueryPerBase(uint32 windowSize);

public:
  uint32                     paintSave  = 0;
};



posLookGlobal::posLookGlobal() {
}
posLookGlobal::~posLookGlobal() {
  delete    inputFile;
  delete    refLookup;
  delete [] nQseqPer;
  delete [] nQmerPer;

  delete    hitsPerQuery;
  delete    paintMerPerBase;
  delete    paintQueryPerBase;
}



void
posLookGlobal::initialize(void) {
  merylFileReader   *refKmers    = new merylFileReader(refMerName);
  dnaSeqFile        *refSequence = openSequenceFile(refSeqName);

  refLookup = new merylExactLookup();

  refLookup->load(refKmers, 32, 18);      //  Load kmers, allowing 32GB and 18 bits of prefix.
  refLookup->loadPositions(refSequence);  //  Load positions of those kmers.

  delete refSequence;
  delete refKmers;

  nQseqPer = new uint32 [256 * 1024 * 1024];
  nQmerPer = new uint32 [256 * 1024 * 1024];

  for (uint32 ii=0; ii<256 * 1024 * 1024; ii++)
    nQseqPer[ii] = nQmerPer[ii] = 0;
}


seqBatch *
posLookGlobal::loadBatch(uint64 nSeqs, uint64 nBases) {
  uint64    nS  = 0;
  uint64    nB  = 0;
  seqBatch *s   = nullptr;
  dnaSeq   *seq = nullptr;

  while ((openInputFile() == true) &&   //  While there is an open input file and
         (nS < nSeqs) &&                //  the batch size is small...
         (nB < nBases)) {
    seq = new dnaSeq;

    if (inputFile->loadSequence(*seq) == true) {   //  If we load an input sequence
      if (s == nullptr) {                          //  add it to the seqBatch,
        s = new seqBatch();                        //  creating a seqBatch if necessary.
        s->batchid = batchid++;
      }

      nS += 1;
      nB += seq->length();

      s->sequences.push_back(seq);  seq = nullptr;
    }
    else {                                         //  If nothing loaded, remove the
      delete inputFile;  inputFile = nullptr;      //  empty input file and the
      delete seq;        seq       = nullptr;      //  unused dnaSeq.
    }
  }

  if (s)
    fprintf(stderr, "Loaded batch %4lu with %5lu sequences and %8lu bases.\n", s->batchid, nS, nB);
  return(s);
}

void
posLookGlobal::lookupBatch(seqBatch *s) {

  for (uint32 ii=0; ii<s->sequences.size(); ii++) {
    dnaSeq       *seq = s->sequences[ii];
    kmerIterator  kiter(seq->bases(), seq->length());

    //fprintf(stderr, "  SEQ: %s\r", seq->ident());

    while (kiter.nextMer()) {
      kmer    cmer  = std::min(kiter.fmer(), kiter.rmer());
      uint64  idx   = refLookup->index(cmer);

      if (idx != uint64max) {
        s->hits.push_back( hit_s(uint32max, idx, ii, kiter.bgnPosition()) );
      }
      else {
      }
    }
  } 
}

#if 0
void
posLookGlobal::writeHits(uint32 *nPer, uint64 idx, bool doFullWrite) {
  uint64  base = refLookup->_posStart->get(idx);
  uint32  nmax = refLookup->valueAtIndex(idx);

  for (uint32 nn=0; nn<nmax; nn++) {
    uint64 spp  = refLookup->_posData->get(base + nn);
    uint64 sID  = refLookup->decodeID(spp);
    uint64 sPos = refLookup->decodePos(spp);

    if (nQseqPer)   nQseqPer[sPos]++;

    if (doFullWrite)
      fprintf(stdout, "%lu %lu\n", sID, sPos);
  }
}
#endif

void
posLookGlobal::writeBatch(seqBatch *s) {
  bool   doFullWrite = false;
  uint32 hh=0;

  //if (s->hits.size() == 0)
  //  return;

  fprintf(stderr, "Write  batch %4lu with ----- sequences and -------- bases -- %9lu hits\n", s->batchid, s->hits.size());

  //  Scan the list of hits, count the number of hits to each contig.

  if (hitsPerQuery) {
    uint32 *tCov = new uint32 [s->sequences.size()];
    uint32 *nPer = new uint32 [s->sequences.size()];

    for (uint32 ii=0; ii<s->sequences.size(); ii++)
      tCov[ii] = nPer[ii] = 0;

    for (uint32 hh=0; hh < s->hits.size(); hh++) {
      uint32  qryid = s->hits[hh].qryID;
      uint64  idx   = s->hits[hh].refPos;
      uint64  base  = refLookup->_posStart->get(idx);
      uint32  nmax  = refLookup->valueAtIndex(idx);

      //  Each hit is one kmer in one query to the reference.  That kmer
      //  occurs nmax times in the reference.

      tCov[qryid] += 1;
      nPer[qryid] += nmax;
    }

    for (uint32 ii=0; ii<s->sequences.size(); ii++)
      fprintf(hitsPerQuery->file(), "%u\t%u\t%lu\t%s\n", nPer[ii], tCov[ii], s->sequences[ii]->length(), s->sequences[ii]->ident());

    delete [] nPer;
    delete [] tCov;
  }

  //  Scan the list of hits, count the number of kmers matches to each
  //  reference position.

  if (paintMerPerBase) {
    for (uint32 hh=0; hh < s->hits.size(); hh++) {
      uint32  qryid = s->hits[hh].qryID;
      uint64  idx   = s->hits[hh].refPos;
      uint64  base  = refLookup->_posStart->get(idx);
      uint32  nmax  = refLookup->valueAtIndex(idx);

      //  Over all the nmax positions in the reference, add one
      //  for the hit to this kmer in this contig.

      for (uint32 nn=0; nn<nmax; nn++) {
        uint64 rpp    = refLookup->_posData->get(base + nn);
        uint64 refID  = refLookup->decodeID(rpp);
        uint64 refPos = refLookup->decodePos(rpp);

        nQmerPer[refPos]++;
      }
    }
  }

  //  Scan the list of hits, count the number of distinct contigs that hit
  //  each reference position.

  if (paintQueryPerBase) {
    uint64  lq = uint64max;

    //  Sort by refPos, then by qryID.

    std::sort(s->hits.begin(), s->hits.end(), [](hit_s &a,
                                                 hit_s &b)   { return(((a.refPos  < b.refPos)) ||
                                                                      ((a.refPos == b.refPos) &&
                                                                       (a.qryID   < b.qryID))); });

    //  Scan the hits, incrementing the count of contigs-per-base, skipping
    //  multiple hits from the same contig to the same refPos.

    for (uint32 pp=0; pp<s->hits.size(); pp++) {
      if ((pp != 0) && (s->hits[pp-1].refPos == s->hits[pp].refPos) && (s->hits[pp-1].qryID == s->hits[pp].qryID))
        continue;

      //fprintf(stderr, "");

      //  Over all the nmax positions in the reference, add one
      //  for the hit to this kmer in this contig.

      uint64  idx   = s->hits[pp].refPos;
      uint64  base  = refLookup->_posStart->get(idx);
      uint32  nmax  = refLookup->valueAtIndex(idx);

      for (uint32 nn=0; nn<nmax; nn++) {
        uint64 rpp    = refLookup->_posData->get(base + nn);
        uint64 refID  = refLookup->decodeID(rpp);
        uint64 refPos = refLookup->decodePos(rpp);

        nQseqPer[refPos]++;
      }
    }
  }

  //  Scan the list of hits, output a map from query to reference position.

  {
  }

  fflush(stdout);
}


void
posLookGlobal::saveMerPerBase(uint32 windowSize) {

  if (paintMerPerBase == nullptr)
    return;

  for (uint32 ii=0; ii<256*1024*1024; ii++) {
    if (nQmerPer[ii] > 0)
      fprintf(paintMerPerBase->file(), "%u %u\n", ii, nQmerPer[ii]);
  }
}


void
posLookGlobal::saveQueryPerBase(uint32 windowSize) {
  uint32 winb = 0;
  uint32 wine = windowSize/2;

  if (paintQueryPerBase == nullptr)
    return;

  for (uint32 ii=0; ii<256*1024*1024; ii++) {
    if (nQseqPer[ii] > 0)
      fprintf(paintQueryPerBase->file(), "%u %u\n", ii, nQseqPer[ii]);
  }
}


//
//  sweatShop stubs.
//

void *loadSeqs(void *G) {
  posLookGlobal *g = (posLookGlobal *)G;
  return(g->loadBatch(4 * 1024, 16 * 1048576));
}

void  lookupSeqs(void *G, void *T, void *S) {
  posLookGlobal *g = (posLookGlobal *)G;
  seqBatch      *s = (seqBatch      *)S;
  g->lookupBatch(s);
}

void  writeSeqs(void *G, void *S) {
  posLookGlobal *g = (posLookGlobal *)G;
  seqBatch      *s = (seqBatch      *)S;
  g->writeBatch(s);
  delete s;
}


//
//  Main!
//

int
main(int argc, char **argv) {
  posLookGlobal   *g  = new posLookGlobal;

  for (int32 arg=1; arg<argc; arg++) {
    if      (strcmp(argv[arg], "-m") == 0) {
      g->setRefMerName(argv[++arg]);
    }
    else if (strcmp(argv[arg], "-s") == 0) {
      g->setRefSeqName(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-hpq") == 0) {
      g->hitsPerQuery = new compressedFileWriter(argv[++arg]);
    }
    else if (strcmp(argv[arg], "-mpb") == 0) {
      g->paintMerPerBase = new compressedFileWriter(argv[++arg]);
    }
    else if (strcmp(argv[arg], "-qpb") == 0) {
      g->paintQueryPerBase = new compressedFileWriter(argv[++arg]);
    }

    else if (fileExists(argv[arg])) {
      g->addInputFile(argv[arg]);
    }
    else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      exit(1);
    }
  }

  g->initialize();

  sweatShop      *ss = new sweatShop(loadSeqs, lookupSeqs, writeSeqs);

  ss->setLoaderQueueSize(16);
  ss->setNumberOfWorkers(14);
  ss->setWriterQueueSize(128);
  ss->setInOrderOutput(true);

  //ss->run(g, true);
  ss->run(g, false);

  delete ss;

  g->saveMerPerBase(10000);
  g->saveQueryPerBase(10000);

  delete g;

  return(0);
}
