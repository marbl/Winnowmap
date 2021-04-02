
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
  opDump,
  opExistence,
  opInclude,
  opExclude
};



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

  void initialize(void);
  void loadLookupTables(void);
  void openInputs(void);
  void openOutputs(void);

  //  Inputs.

  uint32                            nThreads     = getMaxThreadsAllowed();
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
  std::vector<merylExactLookup *>   lookupDBs;   //  Kmer lookup table.

  kmvalu                            minV         = 0;
  kmvalu                            maxV         = kmvalumax;

  lookupOp                          reportType   = lookupOp::opNone;

  bool                              doEstimate   = false;
  bool                              showProgress = false;

  //  Outputs for existence.

  char                             *outstring    = nullptr;
  uint32                            outstringMax = 0;

  //  Outputs for include/exclude.

  uint64                            nReadsTotal  = 0;
  uint64                            nReadsFound  = 0;

};



