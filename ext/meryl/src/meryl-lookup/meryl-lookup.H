
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


enum class lookupOp {
  opNone,
  opEstimate,
  opBED,
  opWIGcount,
  opWIGdepth,
  opExistence,
  opInclude,
  opExclude
};

inline
const char *
toString(lookupOp op) {
  switch (op) {
    case lookupOp::opNone:       return("(not supplied)");  break;
    case lookupOp::opEstimate:   return("(estimate)");      break;
    case lookupOp::opBED:        return("-bed");            break;
    case lookupOp::opWIGcount:   return("-wig-count");      break;
    case lookupOp::opWIGdepth:   return("-wig-depth");      break;
    case lookupOp::opExistence:  return("-existence");      break;
    case lookupOp::opInclude:    return("-include");        break;
    case lookupOp::opExclude:    return("-exclude");        break;
  }
  return("(not supplied)");
}



class lookupGlobal {
public:
  lookupGlobal() {
  };

  ~lookupGlobal() {
    delete [] outstring;

    for (uint32 ii=0; ii<lookupDBs.size(); ii++)
      delete lookupDBs[ii];

    delete seqFile1;
    delete seqFile2;

    delete outFile1;
    delete outFile2;
  };

  void checkInvalid(std::vector<char const *> &err);

  void initialize(void);
  void loadLookupTables(void);
  void openInputs(void);
  void openOutputs(void);

  //  Inputs.

  double                            maxMemory    = getMaxMemoryAllowed() / 1024.0 / 1024.0 / 1024.0;

  char const                       *seqName1     = nullptr;
  char const                       *seqName2     = nullptr;

  dnaSeqFile                       *seqFile1     = nullptr;
  dnaSeqFile                       *seqFile2     = nullptr;

  char const                       *outName1     = "-";
  char const                       *outName2     = nullptr;

  compressedFileWriter             *outFile1     = nullptr;
  compressedFileWriter             *outFile2     = nullptr;


  std::vector<const char *>         lookupDBname;
  std::vector<const char *>         lookupDBlabel;
  uint32                            lookupDBlabelLen = 0;
  std::vector<merylExactLookup *>   lookupDBs;   //  Kmer lookup table.

  kmvalu                            minV         = 0;
  kmvalu                            maxV         = kmvalumax;

  lookupOp                          reportType   = lookupOp::opNone;

  bool                              is10x        = false;
  bool                              mergeBedRuns = false;

  bool                              doEstimate   = false;
  bool                              showProgress = false;

  //  Outputs for existence.

  char                             *outstring    = nullptr;
  uint32                            outstringMax = 0;

  //  Outputs for include/exclude.

  uint64                            nReadsTotal  = 0;
  uint64                            nReadsFound  = 0;
};



void helpBED           (char const *progname=nullptr);
void helpWIGcount      (char const *progname=nullptr);
void helpWIGdepth      (char const *progname=nullptr);
void helpExistence     (char const *progname=nullptr);
void helpIncludeExclude(char const *progname=nullptr);
void help              (char const *progname=nullptr);

void dumpExistence(lookupGlobal *G);
void reportExistence(lookupGlobal *G);
void filter(lookupGlobal *G);
