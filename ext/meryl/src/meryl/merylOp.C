
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

#include "meryl.H"
#ifdef CANU
#include "sqStore.H"
#endif

bool            merylOperation::_onlyConfig   = false;
bool            merylOperation::_showProgress = false;
merylVerbosity  merylOperation::_verbosity    = sayStandard;



merylOperation::merylOperation(void) {
}



FILE *
openPerThreadOutput(FILE *prFile, char *prName, uint32 fileNum) {
  char  T[FILENAME_MAX+1] = { 0 };
  char  N[FILENAME_MAX+1] = { 0 };

  if (prFile != nullptr)   //  If we have a file already, it's stdout.
    return(prFile);

  if (prName == nullptr)   //  If no name, no print requested.
    return(nullptr);

  strncpy(T, prName, FILENAME_MAX);

  char   *pre = T;
  char   *suf = strchr(T, '#');
  uint32  len = 0;

  while ((suf) && (*suf == '#')) {
    *suf = 0;
    len++;
    suf++;
  }

  if (len == 0)
    snprintf(N, FILENAME_MAX, "%s.%02d", prName, fileNum);
  else
    snprintf(N, FILENAME_MAX, "%s%0*d%s", pre, len, fileNum, suf);

  return(AS_UTL_openOutputFile(N));
}


merylOperation::merylOperation(merylOperation *op,
                               uint32 fileNum,
                               uint32 nInputs,
                               uint32 threads, uint64 memory) {

  //  Set our operation and basic parameters.
  _operation     = op->_operation;

  _mathConstant  = op->_mathConstant;
  _threshold     = op->_threshold;
  _fracDist      = op->_fracDist;
  _wordFreq      = op->_wordFreq;

  //  Limit resource usage.
  _maxThreads    = threads;
  _maxMemory     = memory;

  //  Remember which piece we're processing.
  _fileNumber    = fileNum;

  //  Allocate space for the input buffers.
  _actCount      = new kmvalu [nInputs];
  _actIndex      = new uint32 [nInputs];

  //  Set pointers to the database output object.
  _outputO = nullptr;
  _outputP = op->_outputO;
  _writer  = nullptr;

  //  Open per-thread printing output.
  _printer = openPerThreadOutput(op->_printer, op->_printerName, _fileNumber);
}



merylOperation::~merylOperation() {

  clearInputs();

  delete    _stats;

  delete    _outputO;

  assert(_writer == NULL);

  if (_printer != stdout)
    AS_UTL_closeFile(_printer);

  delete [] _printerName;

  delete [] _actCount;
  delete [] _actIndex;
}




void
merylOperation::clearInputs(void) {

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  _inputs.clear();

  _actLen = 0;
}



void
merylOperation::checkInputs(const char *name) {

  if ((_actLen > 1) && ((_operation == opPassThrough) ||
                        (_operation == opLessThan)    ||
                        (_operation == opGreaterThan) ||
                        (_operation == opEqualTo)     ||
                        (_operation == opHistogram))) {
    fprintf(stderr, "merylOp::addInput()-- ERROR: can't add input '%s' to operation '%s': only one input supported.\n",
            name, toString(_operation));
    exit(1);
  }
}



void
merylOperation::addInputFromOp(merylOperation *operation) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from operation '%s' to operation '%s'\n",
            toString(operation->_operation), toString(_operation));

  _inputs.push_back(new merylInput(operation));

  if (_actIndex)
    _actIndex[_actLen++] = _inputs.size() - 1;

  if (operation->_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't be used as an input: it doesn't supply kmers.\n", toString(operation->_operation)), exit(1);

  if (_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't take input from '%s': it can only accept databases.\n", toString(_operation), toString(operation->_operation)), exit(1);

  checkInputs(toString(operation->getOperation()));
}


void
merylOperation::addInputFromDB(char *dbName) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
            dbName, toString(_operation));

  merylFileReader *db = new merylFileReader(dbName);

  _inputs.push_back(new merylInput(db->filename(), db, _fileNumber));

  if (_actIndex)
    _actIndex[_actLen++] = _inputs.size() - 1;

  checkInputs(db->filename());
}


void
merylOperation::addInputFromSeq(char *sqName, bool doCompression) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from file '%s' to operation '%s'\n",
            sqName, toString(_operation));

  dnaSeqFile *sq = new dnaSeqFile(sqName);

  _inputs.push_back(new merylInput(sq->filename(), sq, doCompression));

  if (_actIndex)
    _actIndex[_actLen++] = _inputs.size() - 1;

  if (isCounting() == false)
    fprintf(stderr, "ERROR: operation '%s' cannot use sequence files as inputs.\n", toString(_operation)), exit(1);
}



void
merylOperation::addInputFromCanu(char *sqName, uint32 segment, uint32 segmentMax) {

#ifdef CANU
  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding input from sqStore '%s' to operation '%s'\n",
            sqName, toString(_operation));

  sqStore *store = new sqStore(sqName);

  _inputs.push_back(new merylInput(store->sqStore_path(), store, segment, segmentMax));

  if (_actIndex)
    _actIndex[_actLen++] = _inputs.size() - 1;

  if (isCounting() == false)
    fprintf(stderr, "ERROR: operation '%s' cannot use sqStore as inputs.\n", toString(_operation)), exit(1);
#else
#endif
}



void
merylOperation::addOutput(char *wrName) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding output to file '%s' from operation '%s'\n",
            wrName, toString(_operation));

  if (_outputO)
    fprintf(stderr, "ERROR: already have an output set!\n"), exit(1);

  if (_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't use 'output' modifier.\n", toString(_operation));

  _outputO = new merylFileWriter(wrName);
}



void
merylOperation::addPrinter(char *prName, bool ACGTorder) {

  if (_verbosity >= sayConstruction)
    fprintf(stderr, "Adding printer to %s from operation '%s'\n",
            (prName == nullptr) ? "(stdout)" : prName,
            toString(_operation));

  if (_printerName)
    fprintf(stderr, "ERROR: already have a printer set!\n"), exit(1);

  if (_operation == opHistogram)
    fprintf(stderr, "ERROR: operation '%s' can't use 'output' modifier.\n", toString(_operation));

  //  For stdout, prFile is defined.  For files, it is nullptr.

  if (prName == nullptr) {
    _printer        = stdout;
    _printerName    = duplicateString("(stdout)");
    _printACGTorder = ACGTorder;
  } else {
    _printer        = nullptr;
    _printerName    = duplicateString(prName);
    _printACGTorder = ACGTorder;
  }
}




//  We're all done processing this operation.  Clean up what we can.
//  The _output CANNOT be deleted until all operations are done with it.
//  Yes, I should be pointer counting or something smart like that.
void
merylOperation::finalize(void) {

  clearInputs();
}







char const *
toString(merylOp op) {
  switch (op) {
    case opCount:                return("opCount");                break;
    case opCountForward:         return("opCountForward");         break;
    case opCountReverse:         return("opCountReverse");         break;
    case opPassThrough:          return("opPassThrough");          break;

    case opLessThan:             return("opLessThan");             break;
    case opGreaterThan:          return("opGreaterThan");          break;
    case opAtLeast:              return("opAtLeast");              break;
    case opAtMost:               return("opAtMost");               break;
    case opEqualTo:              return("opEqualTo");              break;
    case opNotEqualTo:           return("opNotEqualTo");           break;

    case opIncrease:             return("opIncrease");             break;
    case opDecrease:             return("opDecrease");             break;
    case opMultiply:             return("opMultiply");             break;
    case opDivide:               return("opDivide");               break;
    case opDivideRound:          return("opDivideRound");          break;
    case opModulo:               return("opModulo");               break;

    case opUnion:                return("opUnion");                break;
    case opUnionMin:             return("opUnionMin");             break;
    case opUnionMax:             return("opUnionMax");             break;
    case opUnionSum:             return("opUnionSum");             break;

    case opIntersect:            return("opIntersect");            break;
    case opIntersectMin:         return("opIntersectMin");         break;
    case opIntersectMax:         return("opIntersectMax");         break;
    case opIntersectSum:         return("opIntersectSum");         break;

    case opSubtract:             return("opSubtract");             break;

    case opDifference:           return("opDifference");           break;
    case opSymmetricDifference:  return("opSymmetricDifference");  break;

    case opHistogram:            return("opHistogram");            break;
    case opStatistics:           return("opStatistics");           break;

    case opCompare:              return("opCompare");              break;

    case opNothing:              return("opNothing");              break;
  }

  assert(0);
  return(NULL);
}

