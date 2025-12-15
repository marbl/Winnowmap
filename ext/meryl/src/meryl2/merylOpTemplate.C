
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


merylOpTemplate::merylOpTemplate(uint32 ident) {
  _ident = ident;
}

merylOpTemplate::~merylOpTemplate() {

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  _inputs.clear();

  for (uint32 ii=0; ii<64; ii++)
    delete _computes[ii];

  delete    _counting;

  delete [] _outDbseName;   delete _outDbse;
  delete [] _outListName;   delete _outList;
  delete [] _outShowName;   delete _outShow;
  delete [] _outPipeName;   //lete _outPipe;

  delete [] _outStatsName;  delete _outStats;
  delete [] _outHistoName;  delete _outHisto;

  delete [] _valueString;
  delete [] _labelString;
}




//  Add a 'database' output to the template command.
//
void
merylOpTemplate::addOutputToDB(char const *wrName, std::vector<char const *> &err) {

  if (globals.showConstruction() == true)
    fprintf(stderr, "addOutputToDB()-- action #%u -> output to database '%s'\n",
            _ident, wrName);

  if (_outDbse != nullptr) {
    sprintf(err, "Operation #%u already writes output to '%s',", _ident, _outDbseName);
    sprintf(err, "  can't add another output to '%s'!", wrName);
    sprintf(err, "");
    return;
  }

  _outDbseName = duplicateString(wrName);
  _outDbse     = new merylFileWriter(wrName);
}



//  Add a 'list' output to the template command.
//
//   If the name is '-':
//    - output will go to stdout.
//
//   If the name does not contain two or more '#' symbols
//    - and output will go to a single file with that name.
//
//   If the name does contain two or more '#' symbols
//    - and output will go to 64 files, one per slice, replacing
//      the ##'s with digits.
//
//  See merylOpCompute::addLister() for how the last two cases are distinguished.
//
void
merylOpTemplate::addOutputToList(char const *prName, bool ACGTorder, std::vector<char const *> &err) {
  bool   isNormal = false;
  bool   isStdout = false;

  //  Catch outputs to "-' and redirect them to 'stdout', by convention,
  //  communicated as prName = nullptr.

  if ((prName != nullptr) && (prName[0] == '-') && (prName[1] == 0))
    prName = nullptr;

  //  Decide if we're outputting to a single normal file, or to parallel files.
  //  Note this only counts the first set of consecutive hashes.

  uint32  nHash = 0;

  if (prName)
    for (char const *suf = strchr(prName, '#'); ((suf) && (*suf == '#')); suf++)
      nHash++;

  //
  //  List kmers to stdout.
  //
  if (prName == nullptr) {
    if (globals.showConstruction() == true)
      fprintf(stderr, "addOutput()-- action #%u -> text list\n", _ident);

    if (_outShow) {
      sprintf(err, "Operation #%u is already showing kmers to '%s',", _ident, _outShowName);
      sprintf(err, "  can't show them twice.");
      sprintf(err, "");
      return;
    }

    _outShowName = duplicateString("(stdout)");
    _outShow     = new compressedFileWriter(nullptr);
  }

  //
  //  List kmers to a single file.
  //
  else if (nHash == 0) {
    if (globals.showConstruction() == true)
      fprintf(stderr, "addOutput()-- action #%u -> text list '%s'\n", _ident, prName);

    if (_outList) {
      sprintf(err, "Operation #%u is already listing kmers to '%s',", _ident, _outListName);
      sprintf(err, "  can't also list them to '%s'.", prName);
      sprintf(err, "");
      return;
    }

    _outListName = duplicateString(prName);
    _outList     = new compressedFileWriter(prName);
    return;
  }

  //
  //  List kmers to parallel files.
  //
  else {
    if (globals.showConstruction() == true)
      fprintf(stderr, "addOutput()-- action #%u -> text list '%s' (parallel mode)\n", _ident, prName);

    if (_outList) {
      sprintf(err, "Operation #%u is already listing kmers to '%s',", _ident, _outListName);
      sprintf(err, "  can't also list them to '%s'.", prName);
      sprintf(err, "");
      return;
    }

    _outListName = duplicateString(prName);
    _outList     = nullptr;
  }
}



//
//  Add 'histogram' or 'statistics' output to this operation.
//  Fails if 

void
merylOpTemplate::addStatsOutput(char const *hiName, std::vector<char const *> &err) {

  if ((hiName == nullptr) || (hiName[0] == 0))
    hiName = "-";

  if (globals.showConstruction() == true)
    fprintf(stderr, "addOutput()-- action #%u -> statistics to '%s'\n", _ident, hiName);

  if (_outStats != nullptr) {
    sprintf(err, "Operation #%u already has 'statistics' output to file '%s',", _ident, _outStats->filename());
    sprintf(err, "  can't add another output to file '%s'.", hiName);
    sprintf(err, "");
    return;
  }

  _outStatsName = duplicateString(hiName);
  _outStats     = new compressedFileWriter(hiName);
}

void
merylOpTemplate::addHistoOutput(char const *hiName, std::vector<char const *> &err) {

  if ((hiName == nullptr) || (hiName[0] == 0))
    hiName = "-";

  if (globals.showConstruction() == true)
    fprintf(stderr, "addOutput()-- action #%u -> histogram to '%s'\n", _ident, hiName);

  if (_outHisto != nullptr) {
    sprintf(err, "Operation #%u already has 'histogram' output to file '%s',", _ident, _outHisto->filename());
    sprintf(err, "  can't add another output to file '%s'.", hiName);
    sprintf(err, "");
    return;
  }

  _outHistoName = duplicateString(hiName);
  _outHisto     = new compressedFileWriter(hiName);
}



//  This is called by merylCommandBuilder after the entire tree has been
//  built, but before any processing threads are created.  It will:
//   - open any input databases (to set the kmer size)
//   - make sure counting actions are from seequence and processing
//     actions are from databases.
//   - let select:input figure out what 'all' means, and fail
//     if an input is not present (asking for input 4, but only
//     3 are supplied).
//   - let the action itself fail if there are too few or too
//     many inputs.
//
//  It is NOT intended to do any processing, like examining histograms to
//  choose thresholds.  That is done in initializeTemplate().
//
void
merylOpTemplate::finalizeTemplateInputs(uint32 oo, std::vector<char const *> &err) {

  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    _inputs[ii]->openInput(err);

    if ((_isCounting == true) &&
        (_inputs[ii]->isFromSequence() == false) &&
        (_inputs[ii]->isFromStore()    == false))
      sprintf(err, "ERROR: counting action at position %u input %u must be from a sequence file or Canu seqStore.", oo, ii);

    if ((_isCounting == false) &&
        (_inputs[ii]->isFromTemplate() == false) &&
        (_inputs[ii]->isFromDatabase() == false))
      sprintf(err, "ERROR: action at position %u input %u must be from a meryl database or another action.", oo, ii);
  }


  for (uint32 f1=0; f1<_select.size(); f1++)
    for (uint32 f2=0; f2<_select[f1].size(); f2++)
      _select[f1][f2].finalizeSelectorInputs(this, err);


  if (_inputs.size() < _inputsMin)
    sprintf(err, "Action at position %u has %u inputs, but requires at least %u.\n",
            oo, _inputs.size(), _inputsMin);

  if (_inputs.size() > _inputsMax)
    sprintf(err, "Action at position %u has %u inputs, but requires at most %u.\n",
            oo, _inputs.size(), _inputsMax);


  if ((_isCounting == true) &&
      (_outDbseName == nullptr))
    sprintf(err, "ERROR: counting operation at position %u must have an output database.\n", oo);
}



//  Counting is done, we're about ready to process.
//   - open any output databases
//   - figure out any last minute parameters (that are found
//     by querying databases, for example).
//
void
merylOpTemplate::finalizeTemplateParameters(void) {

  if (_outDbse)                                    //  Create the master output object.  We'll later
    _outDbse->initialize(0, false);                //  request per-thread writer objects from this.

  for (uint32 f1=0; f1<_select    .size(); f1++)   //  Let selectors query inputs
  for (uint32 f2=0; f2<_select[f1].size(); f2++)   //  for parameters.
    _select[f1][f2].finalizeSelectorParameters(this);
}



void
merylOpTemplate::finishAction(void) {

  //  Forward the request to any inputs that are actions.

  for (uint32 ii=0; ii<_inputs.size(); ii++)
    if (_inputs[ii]->_template != nullptr)
      _inputs[ii]->_template->finishAction();

  //  Gather stats from the compute threads then output.

  if ((_outStats != nullptr) ||
      (_outHisto != nullptr)) {
    merylHistogram  statsTot;

    for (uint32 ss=0; ss<64; ss++)
      statsTot.insert( _computes[ss]->_statsAcc );

    if (_outStats != nullptr)
      statsTot.reportStatistics(_outStats->file());

    if (_outHisto != nullptr)
      statsTot.reportHistogram(_outHisto->file());
  }

  //  Delete the compute objets.

  //  Anything else should be done in the various destructors.
}





void
merylOpTemplate::doCounting(uint64 allowedMemory,
                            uint32 allowedThreads) {

  if ((_counting == nullptr) ||
      (_outDbse == nullptr))
    return;

  _counting->setDefaultLabel(_labelConstant);
  _counting->setDefaultValue(0);
  _counting->doCounting(_inputs, allowedMemory, allowedThreads, _outDbse);

  if (_onlyConfig == true)   //  If only configuring, stop now.
    return;

  //  Convert this op into a pass through
#warning NEED TO RESET COUNTING OPERATION TO PASS THROUGH
  //  Fiddle with the operation.
  //   - remove the output; it's already been written.
  //   - remove all the inputs
  //   - convert the operation to a simple 'pass through'
  //   - add the counted output as an input

  //char *name = duplicateString(_outDbse->filename());

  //char name[FILENAME_MAX + 1];
  //strncpy(name, _outDbse->filename(), FILENAME_MAX + 1);   //  know which input to open later.

  //  Close the inputs and forget about them too.
  for (uint32 ii=0; ii<_inputs.size(); ii++)
    delete _inputs[ii];

  _inputs.clear();

  //  But now remember what that output was and make it an input.
  //    (see addInputFromDB())

  merylInput  *in = new merylInput;

  in->registerMerylDB(_outDbseName);
  addInput(in);

  //  Close the output and forget about it.
  delete [] _outDbseName;   _outDbseName = nullptr;
  delete    _outDbse;       _outDbse     = nullptr;

  delete    _counting;      _counting    = nullptr;
  //_inputs.push_back(new merylInput(new merylFileReader(name)));
};



char const *
merylOpTemplate::displayValueAssignment(char *ts) {
  static char const *formats[26];

  formats[0+2*merylAssignValue::valueNOP]      = "no-operation";
  formats[0+2*merylAssignValue::valueSet]      = "set to constant zero because no constant supplied";
  formats[0+2*merylAssignValue::valueSelected] = "that of the kmer selected by LABEL selector";
  formats[0+2*merylAssignValue::valueFirst]    = "that of the kmer in the first input";
  formats[0+2*merylAssignValue::valueMin]      = "the minimum of all kmer values";
  formats[0+2*merylAssignValue::valueMax]      = "the maximum of all kmer values";
  formats[0+2*merylAssignValue::valueAdd]      = "the sum of all kmer values";
  formats[0+2*merylAssignValue::valueSub]      = "the selected kmer value minus all others";
  formats[0+2*merylAssignValue::valueMul]      = "the product of all kmer values";
  formats[0+2*merylAssignValue::valueDiv]      = "the selected kmer divided by all others";
  formats[0+2*merylAssignValue::valueDivZ]     = "the selected kmer divided by all others, rounding zero up to one";
  formats[0+2*merylAssignValue::valueMod]      = "the remainder of the selected kmer divided by all others";
  formats[0+2*merylAssignValue::valueCount]    = "the number of input databases the kmer is present in";

  formats[1+2*merylAssignValue::valueNOP]      = "no-operation (supplied constant '%s' ignored)";
  formats[1+2*merylAssignValue::valueSet]      = "constant '%s'";
  formats[1+2*merylAssignValue::valueSelected] = "that of the kmer selected by label selector (supplied constant '%s' ignored)";
  formats[1+2*merylAssignValue::valueFirst]    = "that of the kmer in the first input (supplied constant '%s' ignored)";
  formats[1+2*merylAssignValue::valueMin]      = "the minimum of all kmer values and constant '%s'";
  formats[1+2*merylAssignValue::valueMax]      = "the maximum of all kmer values and constant '%s'";
  formats[1+2*merylAssignValue::valueAdd]      = "the sum of all kmer values and constant '%s'";
  formats[1+2*merylAssignValue::valueSub]      = "the selected kmer minus all others and constant '%s'";
  formats[1+2*merylAssignValue::valueMul]      = "the product of all kmer values and constant '%s'";
  formats[1+2*merylAssignValue::valueDiv]      = "the selected kmer divided by all others and constant '%s'";
  formats[1+2*merylAssignValue::valueDivZ]     = "the selected kmer divided by all others and constant '%s', rounding zero up to one";
  formats[1+2*merylAssignValue::valueMod]      = "the remainder of the selected kmer divided by all others and constant '%s'";
  formats[1+2*merylAssignValue::valueCount]    = "the number of input databases the kmer is present in (supplied constant '%s' ignored)";

  uint32      fi = 2 * _valueAssign + (_valueString == nullptr ? 0 : 1);

  sprintf(ts, formats[fi], _valueString);

  return ts;
}


char const *
merylOpTemplate::displayLabelAssignment(char *ts) {
  static char const *formats[34];

  formats[0+2*merylAssignLabel::labelNOP]         = "no-operation";
  formats[0+2*merylAssignLabel::labelSet]         = "set to constant zero because no constant supplied!";
  formats[0+2*merylAssignLabel::labelSelected]    = "that of the kmer selected by VALUE selector";
  formats[0+2*merylAssignLabel::labelFirst]       = "that of the kmer in the first input";
  formats[0+2*merylAssignLabel::labelMin]         = "the minimum of all labels";
  formats[0+2*merylAssignLabel::labelMax]         = "the maximum of all labels";
  formats[0+2*merylAssignLabel::labelAnd]         = "the bitwise AND of all labels";
  formats[0+2*merylAssignLabel::labelOr]          = "the bitwise OR of all labels";
  formats[0+2*merylAssignLabel::labelXor]         = "the bitwise XOR of all labels";
  formats[0+2*merylAssignLabel::labelDifference]  = "the selected kmer label bitwise minus all others";
  formats[0+2*merylAssignLabel::labelLightest]    = "the label with the fewest set bits (the lightest)";
  formats[0+2*merylAssignLabel::labelHeaviest]    = "the label with the most set bits (the heaviest)";
  formats[0+2*merylAssignLabel::labelInvert]      = "the bitwise complement of the selected label";
  formats[0+2*merylAssignLabel::labelShiftLeft]   = "the selected label shifted left one position";
  formats[0+2*merylAssignLabel::labelShiftRight]  = "the selected label shifted right one position";
  formats[0+2*merylAssignLabel::labelRotateLeft]  = "the selected label rotated left one position";
  formats[0+2*merylAssignLabel::labelRotateRight] = "the selected label rotated right one position";

  formats[1+2*merylAssignLabel::labelNOP]         = "no-operation (supplied constant '%s' ignored)";
  formats[1+2*merylAssignLabel::labelSet]         = "set to constant '%s'";
  formats[1+2*merylAssignLabel::labelSelected]    = "that of the kmer selected by VALUE selector (supplied constant '%s' ignored)";
  formats[1+2*merylAssignLabel::labelFirst]       = "that of the kmer in the first input (supplied constant '%s' ignored)";
  formats[1+2*merylAssignLabel::labelMin]         = "the minimum of all labels and constant '%s'";
  formats[1+2*merylAssignLabel::labelMax]         = "the maximum of all labels and constant '%s'";
  formats[1+2*merylAssignLabel::labelAnd]         = "the bitwise AND of all labels and constant '%s'";
  formats[1+2*merylAssignLabel::labelOr]          = "the bitwise OR of all labels and constant '%s'";
  formats[1+2*merylAssignLabel::labelXor]         = "the bitwise XOR of all labels and constant '%s'";
  formats[1+2*merylAssignLabel::labelDifference]  = "the selected kmer label bitwise minus all others and constant '%s'";
  formats[1+2*merylAssignLabel::labelLightest]    = "the label or constant '%s' with the fewest set bits (the lightest)";
  formats[1+2*merylAssignLabel::labelHeaviest]    = "the label or constant '%s' with the most set bits (the heaviest)";
  formats[1+2*merylAssignLabel::labelInvert]      = "the bitwise complement of the selected label (supplied constant '%s' ignored)";
  formats[1+2*merylAssignLabel::labelShiftLeft]   = "the selected label shifted left '%s' positions";
  formats[1+2*merylAssignLabel::labelShiftRight]  = "the selected label shifted right '%s' positions";
  formats[1+2*merylAssignLabel::labelRotateLeft]  = "the selected label rotated left '%s' positions";
  formats[1+2*merylAssignLabel::labelRotateRight] = "the selected label rotated right '%s' positions";

  uint32      fi = 2 * _labelAssign + (_labelString == nullptr ? 0 : 1);

  sprintf(ts, formats[fi], _labelString);

  return ts;
}


