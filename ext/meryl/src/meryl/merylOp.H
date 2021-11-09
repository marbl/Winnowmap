
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

#ifndef MERYLOP_H
#define MERYLOP_H

#include "meryl.H"

//  These would (maybe) be handy, but they change the order of kmers, so we need
//  to sort again.
//
//  opReverse,             //  Reverse the kmer.     Accepts a single input only.
//  opComplement,          //  Complement the kmer.  Accepts a single input only.
//  opReverseComplement,   //  Rev-Compl the kmer.   Accepts a single input only.


enum merylOp {
  opCount,
  opCountForward,
  opCountReverse,
  opPassThrough,         //  Internal operation; behaves just like opUnion, but only a single input.

  opLessThan,
  opGreaterThan,
  opAtLeast,
  opAtMost,
  opEqualTo,
  opNotEqualTo,

  opIncrease,
  opDecrease,
  opMultiply,
  opDivide,
  opDivideRound,
  opModulo,

  opUnion,
  opUnionMin,
  opUnionMax,
  opUnionSum,

  opIntersect,
  opIntersectMin,
  opIntersectMax,
  opIntersectSum,

  opDifference,
  opSymmetricDifference,

  opSubtract,	// if count(a) >= count(b), keep count(a)-count(b). else, discard.

  opHistogram,
  opStatistics,

  opCompare,

  opNothing
};


enum merylVerbosity {
  sayNothing,
  sayStandard,
  sayConstruction,
  sayDetails,
  sayEverything
};


//  For the simple counting algorithm, how bit of a word to use for initial counts.
typedef uint16 lowBits_t;


class merylOperation {
public:
  merylOperation();
  merylOperation(merylOperation *op,
                 uint32 fileNum,
                 uint32 nInputs,
                 uint32 threads, uint64 memory);
  ~merylOperation();

private:
  void    clearInputs(void);
  void    checkInputs(const char *name);

  uint64  guesstimateNumberOfkmersInInput_dnaSeqFile(dnaSeqFile *sequence);
  uint64  guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID);
  uint64  guesstimateNumberOfkmersInInput(void);

  //  Returns true
  void    configureCounting(uint64   memoryAllowed,      //  Input:  Maximum allowed memory in bytes
                            bool    &useSimple_,         //  Output: algorithm to use
                            uint32  &wPrefix_,           //  Output: Number of bits in the prefix (== bucket address)
                            uint64  &nPrefix_,           //  Output: Number of prefixes there are (== number of buckets)
                            uint32  &wData_,             //  Output: Number of bits in kmer data
                            kmdata  &wDataMask_);        //  Output: A mask to return just the data of the mer);

public:
  void    addInputFromOp  (merylOperation *operation);
  void    addInputFromDB  (char *dbName);
  void    addInputFromCanu(char *stName, uint32 segment, uint32 segmentMax);
  void    addInputFromSeq (char *sqName, bool doCompression);

  void    addOutput(char *wrName);
  void    addPrinter(char *prName, bool ACGTorder);



  void    finalize(void);

  void    setOperation(merylOp op) { _operation = op;    };
  merylOp getOperation(void)       { return(_operation); };

  void    setConstant(uint64 p)                { _mathConstant = p; };
  void    setThreshold(uint64 p)               { _threshold    = p; };
  void    setFractionDistinct(double p)        { _fracDist     = p; };
  void    setWordFrequency(double p)           { _wordFreq     = p; };

  void    setExpectedNumberOfKmers(uint64 n)   { _expNumKmers  = n; };

  void    setCountSuffix(char *s)              {
    _countSuffixLength = strlen(s);

    for (uint32 ii=0; ii<_countSuffixLength; ii++) {
      _countSuffixString[ii] = s[ii];

      _countSuffix.addR(s[ii]);
    }
  };

  void    setMemoryLimit(uint64 m)             { _maxMemory    = m; };
  void    setThreadLimit(uint32 t)             { _maxThreads   = t; };

  uint64  getMemoryLimit(void)                 { return(_maxMemory);  };
  uint32  getThreadLimit(void)                 { return(_maxThreads); };

  bool    isCounting(void) {
    return((_operation == opCount)        ||
           (_operation == opCountForward) ||
           (_operation == opCountReverse));
  };

  bool    isNormal(void) {
    return(isCounting() == false);
  };

  bool    needsThreshold(void) {
    return((_operation == opLessThan)     ||
           (_operation == opGreaterThan)  ||
           (_operation == opAtLeast)      ||
           (_operation == opAtMost)       ||
           (_operation == opEqualTo)      ||
           (_operation == opNotEqualTo));
  };

  bool    needsConstant(void) {
    return((_operation == opIncrease)     ||
           (_operation == opDecrease)     ||
           (_operation == opMultiply)     ||
           (_operation == opDivide)       ||
           (_operation == opDivideRound)  ||
           (_operation == opModulo));
  };

  //bool    needsRealParameter(void) {
  //  return((_operation == opMostFrequent));
  //};

  bool    isMultiSet(void)         { return(_isMultiSet); };

  kmer   &theFMer(void)            { return(_kmer);       };
  kmvalu  theValue(void)           { return(_value);      };

public:
  void    doCounting(void);

public:
  void    initializeThreshold(void);
  bool    initialize(bool isMasterTree=false);

private:
  void    nextMer_findSmallestNormal(void);
  void    nextMer_findSmallestMultiSet(void);
  bool    nextMer_finish(void);

public:
  bool    nextMer(void);

  bool    validMer(void)           { return(_valid);  };

  void    countSimple(void);
  void    countThreads(uint32  wPrefix,
                       uint64  nPrefix,
                       uint32  wData,
                       kmdata  wDataMask);
  void    count(uint32  wPrefix,
                uint64  nPrefix,
                uint32  wData,
                kmdata  wDataMask);

  void    reportHistogram(void);
  void    reportStatistics(void);

  static
  void    onlyConfigure(void)               { _onlyConfig   = true;       };
  static
  void    showProgress(void)                { _showProgress = true;       };
  static
  void    increaseVerbosity(void) {
    if      (_verbosity == sayNothing)        _verbosity    = sayStandard;
    else if (_verbosity == sayStandard)       _verbosity    = sayConstruction;
    else if (_verbosity == sayConstruction)   _verbosity    = sayDetails;
    else if (_verbosity == sayDetails)        _verbosity    = sayEverything;
    else                                      _verbosity    = sayEverything;
  };
  static
  void    beQuiet(void)                     { _verbosity    = sayNothing; };


private:
  void    findMinCount(void);
  void    findMaxCount(void);
  void    findSumCount(void);
  void    subtractCount(void);

  std::vector<merylInput *>      _inputs;
  bool                           _isMultiSet = false;

  merylOp                        _operation = opNothing;

  uint64                         _mathConstant = 0;
  uint64                         _threshold    = UINT64_MAX;
  double                         _fracDist     = DBL_MAX;
  double                         _wordFreq     = DBL_MAX;

  char                           _countSuffixString[65] = {0};
  uint32                         _countSuffixLength     =  0;
  kmer                           _countSuffix;

  uint64                         _expNumKmers = 0;

  uint32                         _maxThreads = 0;
  uint64                         _maxMemory  = 0;

  merylHistogram                *_stats = nullptr;

  merylFileWriter               *_outputO = nullptr;    //  Main output object (this pointer owns the object)
  merylFileWriter               *_outputP = nullptr;    //  Main output object
  merylStreamWriter             *_writer  = nullptr;    //  Per-thread output

  FILE                          *_printer        = nullptr;
  char                          *_printerName    = nullptr;
  bool                           _printACGTorder = false;

  char                           _kmerString[65] = {0};

  uint32                         _fileNumber = UINT32_MAX;

  uint32                         _actLen   = 0;
  kmvalu                        *_actCount = nullptr;
  uint32                        *_actIndex = nullptr;

  kmer                           _kmer;
  kmvalu                         _value = 0;
  bool                           _valid = true;

  static bool                    _onlyConfig;
  static bool                    _showProgress;
  static merylVerbosity          _verbosity;

  friend class merylCommandBuilder;
};


char const *toString(merylOp op);


#endif  //  MERYLOP_H
