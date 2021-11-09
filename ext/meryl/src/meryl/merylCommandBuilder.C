
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

#include "runtime.H"
#include "strings.H"
#include "system.H"



//  In meryOp-count.C
uint64
findMaxInputSizeForMemorySize(uint32 kMerSize, uint64 memorySize);



bool
isDigit(char c) {
  return(('0' <= c) && (c <= '9'));
}

bool
isNumber(char *s, char dot='.') {

  if ((s    == NULL) ||
      (s[0] == 0))
    return(false);

  for (uint32 ii=0; s[ii] != 0; ii++)
    if ((isDigit(s[ii]) == false) &&
        (s[ii] != dot))
      return(false);

  return(true);
}






//  Everything is initialized in the declaration.  Nothing really to do here.
merylCommandBuilder::merylCommandBuilder() {
  _allowedThreads = getMaxThreadsAllowed();   //  Absolute maximum limits on
  _allowedMemory  = getMaxMemoryAllowed();    //  memory= and threads= values
}




merylCommandBuilder::~merylCommandBuilder() {

  for (uint32 ii=0; ii<_opRoot.size(); ii++) {
    uint32  rr = _opRoot[ii];

    for (uint32 tt=0; tt<64; tt++)   //  Destroy threads first.
      delete _thList[tt][rr];

    delete _opList[rr];              //  Then destroy the master.
  }

  for (uint32 tt=0; tt<64; tt++)
    delete [] _thList[tt];
}





void
merylCommandBuilder::terminateOperation(void) {

  for (; _terminating > 0; _terminating--) {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "terminate()-  Pop operation '%s' from top of stack.\n", toString(_opStack.top()->getOperation()));
    _opStack.pop();
  }
}





//  Initialize, either for a new operation or just for a new option string.
//
//  Strip off any leading '['s, and count the number of closing ']'s.
//
//  Save some copies of the stripped string, converted to file paths,
//  for later use.
//
//  Make sure there is an operation on the stack.
//
void
merylCommandBuilder::initialize(char *opt) {

  if (merylOperation::_verbosity >= sayConstruction)
    fprintf(stderr, "initialize()-\n");

  //  Process any existing left-over termination requests.

  terminateOperation();

  //  Save a copy of the string.

  _optStringLen = 0;

  while ((_optStringLen < FILENAME_MAX) && (opt[_optStringLen] != 0)) {
    _optString[_optStringLen] = opt[_optStringLen];
    _optStringLen++;
  }

  _optString[_optStringLen] = 0;

  //  Ignore '[' at the start of the string.  Their purpose is to visually
  //  match the ']' which tells us to stop adding inputs to the current
  //  command.  There should only be one opening bracket.

  if (_optString[0] == '[') {
    for (uint32 ii=0; ii<_optStringLen; ii++)
      _optString[ii] = _optString[ii+1];

    _optStringLen--;
  }

  //  If we have a ']' as the last character, strip it off and remember that
  //  we need to close the command on the stack after we process this arg.
  //  We can get any number of closing brackets.

  while ((_optStringLen > 0) &&
         (_optString[_optStringLen-1] == ']')) {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "initialize()- Found terminator.\n");

    _optString[_optStringLen-1] = 0;
    _optStringLen--;

    _terminating++;
  }

  //  Save a few copies of the command line word.

  strncpy(_inoutName, _optString, FILENAME_MAX + 1);

  snprintf(_indexName, FILENAME_MAX, "%s/merylIndex", _optString);
  snprintf(_sqInfName, FILENAME_MAX, "%s/info",       _optString);
  snprintf(_sqRdsName, FILENAME_MAX, "%s/reads",      _optString);

  //  If the stack is empty, push on a new operation.

  if (_opStack.size() == 0) {
    merylOperation  *op = new merylOperation;

    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "initialize()- Add new empty root operation at position %zu.\n", _opList.size());

    _opStack.push(op);
    _opList.push_back(op);
    _opRoot.push_back(_opList.size() - 1);
  }

  if (merylOperation::_verbosity >= sayConstruction)
    fprintf(stderr, "initialize()- Got option '%s' for '%s' at stack level %zu.\n",
            _optString,
            toString(_opStack.top()->getOperation()),
            _opStack.size());
}





bool
merylCommandBuilder::processOptions(void) {
  KeyAndValue   kv(_optString);
  char         *key = kv.key();
  char         *val = kv.value();

  //  Check for global options.

  if (strncmp(_optString, "-V", 2) == 0) {           //  Anything that starts with -V
    for (uint32 vv=1; vv<strlen(_optString); vv++)   //  increases verbosity by the
      merylOperation::increaseVerbosity();           //  number of letters.
    return(true);
  }

  if (strncmp(_optString, "-Q", 3) == 0) {
    merylOperation::beQuiet();
    return(true);
  }

  if (strncmp(_optString, "-P", 3) == 0) {
    merylOperation::showProgress();
    return(true);
  }

  if (strncmp(_optString, "-C", 3) == 0) {
    merylOperation::onlyConfigure();
    return(true);
  }

  //  If the string is entirely a number, treat it as either a threshold or a
  //  constant, depending on the operation.  This is used for things like
  //  "greater-than 45" and "divide 2".
  //
  //  If there is no operation, or it doesn't want a number, we fall trhough
  //  and return 'false' when key/val is checked below.

  bool  isNum = isNumber(_optString, 0);

  if ((_opStack.top()->needsThreshold() == true) && (isNum == true)) {
    _opStack.top()->setThreshold(strtouint64(_optString));
    return(true);
  }

  if ((_opStack.top()->needsConstant() == true) && (isNum == true)) {
    _opStack.top()->setConstant(strtouint64(_optString));
    return(true);
  }

  //  Some options have no value, unfortunately.

  if (strcmp(_optString, "compress") == 0) {
    _doCompression = true;
    return(true);
  }

  //  The rest should be in a key=value structure.  If we don't find that,
  //  just return.

  if ((key == NULL) ||
      (val == NULL))
    return(false);

  uint32 val32 = strtouint32(val);
  uint64 val64 = strtouint64(val);
  double valDB = strtodouble(val);

  //  Kmer size.
  if (strcmp(key, "k") == 0) {
    if ((kmerTiny::merSize() != 0) &&
        (kmerTiny::merSize() != val32)) {
      fprintf(stderr, "ERROR: kmer size mismatch: %u != %u\n", kmerTiny::merSize(), val32);
      exit(1);
    }
    kmerTiny::setSize(val32);
    return(true);
  }

  //  Number of kmers expected for counting.
  if (strcmp(key, "n") == 0) {
    _opStack.top()->setExpectedNumberOfKmers(val64);
    return(true);
  }

  //  A suffix to filter kmers by when counting.
  if (strcmp(key, "count-suffix") == 0) {
    _opStack.top()->setCountSuffix(val);
    return(true);
  }

  //  Threshold values for less-than, etc, specifed as a fraction of the
  //  total distinct kmers, or as a word-frequency, or as an absolute count.

  if ((strcmp(key, "d") == 0) ||
      (strcmp(key, "distinct") == 0)) {
    _opStack.top()->setFractionDistinct(valDB);
    return(true);
  }

  if ((strcmp(key, "f") == 0) ||
      (strcmp(key, "word-frequency") == 0)) {
    _opStack.top()->setWordFrequency(valDB);
    return(true);
  }

  if ((strcmp(key, "t") == 0) ||            //  See above for special case of this!
      (strcmp(key, "threshold") == 0)) {
    _opStack.top()->setThreshold(val64);
    return(true);
  }

  //  Memory limit, in GB, either global or per-task.

  if (strcmp(key, "memory") == 0) {
    _allowedMemory = (uint64)(valDB * 1024 * 1024 * 1024);
    return(true);
  }

  //  Thread limit, either global or per-task.

  if (strcmp(key, "threads") == 0) {
    _allowedThreads = val32;
    setNumThreads(_allowedThreads);
    return(true);
  }

  //  Segment of input, for counting from seqStore.  Useless otherwise.
  if ((strcmp(key, "segment") == 0) &&
      (isNumber(val, '/'))) {
    decodeRange(val, _segment, _segmentMax);
#ifndef CANU
    fprintf(stderr, "WARNING: option '%s' ignored, available only with Canu support.\n", _optString);
#endif
    return(true);
  }

  //  If nothing triggered, we don't recognize it as an option.  Maybe the
  //  filename had an '=' in it?

  return(false);
}





bool
merylCommandBuilder::processOperation(void) {
  merylOp     non = opNothing;   //  The new op name.

  assert(_opStack.size() > 0);

  //  If the string is of length zero, explicitly do nothing and return
  //  success.  This was (probably) just an isolated terminating bracket.

  if      (0 == _optStringLen)
    return(true);

  //  Check for an operation string.

  if      (0 == strcmp(_optString, "count"))                  non = opCount;
  else if (0 == strcmp(_optString, "count-forward"))          non = opCountForward;
  else if (0 == strcmp(_optString, "count-reverse"))          non = opCountReverse;

  else if (0 == strcmp(_optString, "less-than"))              non = opLessThan;
  else if (0 == strcmp(_optString, "greater-than"))           non = opGreaterThan;
  else if (0 == strcmp(_optString, "at-least"))               non = opAtLeast;
  else if (0 == strcmp(_optString, "at-most"))                non = opAtMost;
  else if (0 == strcmp(_optString, "equal-to"))               non = opEqualTo;
  else if (0 == strcmp(_optString, "not-equal-to"))           non = opNotEqualTo;

  else if (0 == strcmp(_optString, "increase"))               non = opIncrease;
  else if (0 == strcmp(_optString, "decrease"))               non = opDecrease;
  else if (0 == strcmp(_optString, "multiply"))               non = opMultiply;
  else if (0 == strcmp(_optString, "divide"))                 non = opDivide;
  else if (0 == strcmp(_optString, "divide-round"))           non = opDivideRound;
  else if (0 == strcmp(_optString, "modulo"))                 non = opModulo;

  else if (0 == strcmp(_optString, "union"))                  non = opUnion;
  else if (0 == strcmp(_optString, "union-min"))              non = opUnionMin;
  else if (0 == strcmp(_optString, "union-max"))              non = opUnionMax;
  else if (0 == strcmp(_optString, "union-sum"))              non = opUnionSum;

  else if (0 == strcmp(_optString, "intersect"))              non = opIntersect;
  else if (0 == strcmp(_optString, "intersect-min"))          non = opIntersectMin;
  else if (0 == strcmp(_optString, "intersect-max"))          non = opIntersectMax;
  else if (0 == strcmp(_optString, "intersect-sum"))          non = opIntersectSum;

  else if (0 == strcmp(_optString, "subtract"))               non = opSubtract;
  
  else if (0 == strcmp(_optString, "difference"))             non = opDifference;
  else if (0 == strcmp(_optString, "symmetric-difference"))   non = opSymmetricDifference;

  else if (0 == strcmp(_optString, "histogram"))              non = opHistogram;
  else if (0 == strcmp(_optString, "statistics"))             non = opStatistics;

  else if (0 == strcmp(_optString, "compare"))                non = opCompare;

  else return(false);   //  optString is not an operation.

  //  If the top-of-stack command is counting, pop it off and possibly make a
  //  new command.  Counting operations cannot take input from a command.
  //  We're also guaranteed to have either an empty stack, or a valid
  //  non-counting operation on the stack after this.

  if (_opStack.top()->isCounting() == true) {
    _opStack.pop();

    if (_opStack.size() == 0) {
      merylOperation   *op = new merylOperation();

      if (merylOperation::_verbosity >= sayConstruction)
        fprintf(stderr, "processOp()-  Add new empty root operation at position %zu.\n", _opList.size());

      _opStack.push(op);
      _opList.push_back(op);
      _opRoot.push_back(_opList.size() - 1);
    }

    assert(_opStack.top()->isCounting() == false);
  }

  //  If there is a valid command on the stack, push a new command onto the
  //  stack, and add it to the inputs list.

  if (_opStack.top()->getOperation() != opNothing) {
    merylOperation   *op = new merylOperation();

    _opStack.top()->addInputFromOp(op);

    _opStack.push(op);
    _opList.push_back(op);

    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "processOp()-  Operation '%s' added to stack at level %zu\n", toString(non), _opStack.size());
  } else {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "processOp()-  Operation '%s' replaces '%s' at level %zu\n",
              toString(non), toString(_opStack.top()->getOperation()), _opStack.size());
  }

  //  Set the type of this operation.

  _opStack.top()->setOperation(non);

  return(true);  //  optString was an operation.
}




bool
merylCommandBuilder::isOutput(void) {

  //  If we see 'output', flag the next arg as the output name.
  if (strcmp(_optString, "output") == 0) {
    _isOutput = true;
    return(true);
  }

  //  If the flag isn't set, this wasn't 'output' nor is it an output name.
  if (_isOutput == false)
    return(false);

  //  Must be an output name.  Reset the flag, and add the output
  //  to whatever operation is current.

  _isOutput   = false;

  _opStack.top()->addOutput(_inoutName);

  return(true);
}



bool
merylCommandBuilder::isPrinter(void) {

  //  If we see 'print' or 'printACGT', flag the next arg as the output name.
  if (strcmp(_optString, "print") == 0) {
    _printACGTorder = false;
    _isPrint        = true;
    return(true);
  }

  if (strcmp(_optString, "printACGT") == 0) {
    _printACGTorder = true;
    _isPrint        = true;
    return(true);
  }

  //  If the flag isn't set, this isn't a printer output name.
  if (_isPrint == false)
    return(false);

  //  Must be a printer name.  Reset the flag.  Add the name, unless it looks
  //  like a meryl database; in this case, the user requested 'print
  //  db.meryl' and output should go to stdout.  We need to add the printer
  //  to this operation, and return false ("this is not a printer option") so
  //  we can properly handle the database.

  _isPrint = false;

  if (fileExists(_indexName) == true) {
    _opStack.top()->addPrinter(nullptr, _printACGTorder);
    return(false);
  }

  _opStack.top()->addPrinter(_inoutName, _printACGTorder);
  return(true);
}




bool
merylCommandBuilder::isMerylInput(void) {

  if (fileExists(_indexName) == false)
    return(false);

  _opStack.top()->addInputFromDB(_inoutName);

  return(true);
}

bool
merylCommandBuilder::isCanuInput(std::vector<char *> &err) {

  if ((fileExists(_sqInfName) == false) ||
      (fileExists(_sqRdsName) == false))
    return(false);

#ifndef CANU
  char *s = new char [FILENAME_MAX + 129];
  snprintf(s, FILENAME_MAX + 129, "Detected seqStore input '%s', but no Canu support available.", _inoutName);
  err.push_back(s);
#endif

  _opStack.top()->addInputFromCanu(_inoutName, _segment, _segmentMax);

  _segment    = 1;
  _segmentMax = 1;

  return(true);
}

bool
merylCommandBuilder::isSequenceInput(void) {

  if (fileExists(_inoutName) == false)
    return(false);

  _opStack.top()->addInputFromSeq(_inoutName, _doCompression);

  _doCompression = false;

  return(true);
}





void
merylCommandBuilder::finalize(void) {

  //  Finish processing the last option string.
  terminateOperation();

  //  If no operation supplied, assume it's a print, and make it print all kmers.
  if ((_opStack.size() > 0) &&
      (_opStack.top()->getOperation() == opNothing)) {
    if (merylOperation::_verbosity >= sayConstruction)
      fprintf(stderr, "finalize()- Change opNothing to opLessThan at stack level %zu.\n", _opStack.size());
    _opStack.top()->setOperation(opLessThan);
  }

  //  Clear the stack.
  while (_opStack.size() > 0)
    _opStack.pop();

  //  Update memory and threads for everything.
  for (uint32 oo=0; oo<_opList.size(); oo++) {
    _opList[oo]->_maxMemory  = _allowedMemory;
    _opList[oo]->_maxThreads = _allowedThreads;
  }
}



void
merylCommandBuilder::printTree(merylOperation *op, uint32 indent) {

  fprintf(stderr, "%*s%-s\n", indent, "", toString(op->getOperation()));

  if (op->_mathConstant > 0)
    fprintf(stderr, "%*sconstant=%lu\n", indent+2, "", op->_mathConstant);
  if (op->_threshold != UINT64_MAX)
    fprintf(stderr, "%*sthreshold=%lu\n", indent+2, "", op->_threshold);
  if (op->_fracDist != DBL_MAX)
    fprintf(stderr, "%*sfraction-distinct=%f\n", indent+2, "", op->_fracDist);
  if (op->_wordFreq != DBL_MAX)
    fprintf(stderr, "%*sword-frequenct=%f\n", indent+2, "", op->_wordFreq);

  for (uint32 ii=0; ii<op->_inputs.size(); ii++) {
    merylInput  *in = op->_inputs[ii];

    if (in->isFromOperation() == true) {
      printTree(in->_operation, indent+2);
    }

    if (in->isFromDatabase() == true) {
      fprintf(stderr, "%*s%s\n", indent+2, "", in->_name);
    }

    if (in->isFromSequence() == true) {
      fprintf(stderr, "%*s%s%s\n", indent+2, "", in->_name, in->_homopolyCompress ? " (homopoly compressed)" : "");
    }

    if (in->isFromStore() == true) {
      fprintf(stderr, "%*s%s (reads %u through %u)\n", indent+2, "", in->_name, in->_sqBgn, in->_sqEnd);
    }
  }

  if (op->_outputO) {
    fprintf(stderr, "%*soutput to %s\n", indent+2, "", op->_outputO->filename());
  }

  if (op->_printerName) {
    fprintf(stderr, "%*sprint to %s\n", indent+2, "", op->_printerName);
  }
}



//  Clone the command tree(s) into thread-specific copies, one tree per thread.
//
//
void
merylCommandBuilder::spawnThreads(void) {
  uint32  indent = 0;

  setNumThreads(_allowedThreads);

  for (uint32 tt=0; tt<64; tt++) {

    //  Construct a list of operations for each thread.
    _thList[tt] = new merylOperation * [_opList.size()];

    //  Copy operations from the main list to our thread list.
    for (uint32 oo=0; oo<_opList.size(); oo++)
      _thList[tt][oo] = new merylOperation(_opList[oo],
                                           tt,
                                           _opList[oo]->_inputs.size(),
                                           _allowedThreads, _allowedMemory);

    //  Update all the input/output files to be per-thread.
    for (uint32 oo=0; oo<_opList.size(); oo++) {
      merylOperation  *op = _thList[tt][oo];   //  The per-thread operation we're fixing up.
      merylOperation  *OP = _opList[oo];       //  The master operation we're using as a template.

      for (uint32 ii=0; ii<OP->_inputs.size(); ii++) {
        merylInput  *IN = OP->_inputs[ii];     //  The template input for this operation.

        //  If the template input is from an operation, we need to search for
        //  that operation in the master list of operations, then set the
        //  per-thread input to be the corresponding operation.
        if (IN->isFromOperation() == true) {
          uint32  inop = UINT32_MAX;

          for (uint32 xx=0; xx<_opList.size(); xx++)   //  Search all template operations for
            if (IN->_operation == _opList[xx])         //  the one that is our input.
              inop = xx;

          if (inop == UINT32_MAX)
            fprintf(stderr, "Failed to find corresponding operation.\n"), exit(1);

          op->addInputFromOp(_thList[tt][inop]);       //  Add input from the same op in our list.
        }

        //  If the template input is from a database, make a new input for
        //  just the piece we're processing in this thread (done implicitly
        //  in addInputFromDB()).
        if (IN->isFromDatabase() == true) {
          op->addInputFromDB(IN->_name);
        }

        //  We should never get inputs from a sequence file.
        if (IN->isFromSequence() == true) {
          assert(0);
          continue;
        }

        //  We should never get inputs from a Canu seqStore file.
        if (IN->isFromStore() == true) {
          assert(0);
          continue;
        }
      }
    }
  }
}
