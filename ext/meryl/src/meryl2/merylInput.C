
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


merylInput::~merylInput() {
  //lete    _template;    //  We do not own this.
  delete    _action;
  //lete    _pipe;
  delete    _db;
  delete    _list;
  delete    _sequence;
  delete    _store;
  delete    _read;
  delete [] _templateName;
  delete [] _actionName;
  delete [] _pipeName;
  delete [] _dbName;
  delete [] _listName;
  delete [] _sequenceName;
  delete [] _storeName;
}

//  Lazy creation of input sources.  Save whatever is given to us
//  so we can open the input after we've verified the action tree
//  is sane.

bool
merylInput::registerTemplate(merylOpTemplate *t, std::vector<char const *> &err, bool displayErrors) {
  _type = merylInputType::inTemplate;
  _template = t;
  return true;
}

bool
merylInput::registerAction(merylOpCompute *a, std::vector<char const *> &err, bool displayErrors) {
  _type = merylInputType::inAction;
  _action = a;
  return true;
}

bool
merylInput::registerActionPipe(char const *name, std::vector<char const *> &err, bool displayErrors) {
  _type = merylInputType::inPipe;
  _pipeName = duplicateString(name);
  return true;
}

bool
merylInput::registerMerylDB(char const *path, std::vector<char const *> &err, bool displayErrors) {
  if (isMerylDatabase(path) == false) {
    if (displayErrors)
      sprintf(err, "'%s' does not appear to be a meryl database\n", path);
    return false;
  }
  _type = merylInputType::inDB;
  _dbName = duplicateString(path);
  return true;
}

merylInput *
merylInput::registerMerylDB(char const *path) {
  _type = merylInputType::inDB;
  _dbName = duplicateString(path);
  return this;
}

bool
merylInput::registerMerylList(char const *file, std::vector<char const *> &err, bool displayErrors) {
  if (fileExists(file) == false) {
    if (displayErrors)
      sprintf(err, "list file '%s' does not exist\n", file);
    return false;
  }
  _type = merylInputType::inList;
  _listName = duplicateString(file);
  return true;
}

bool
merylInput::registerSeqFile(char const *file, bool doCompression, std::vector<char const *> &err, bool displayErrors) {
  if (fileExists(file) == false) {
    if (displayErrors)
      sprintf(err, "sequence file '%s' does not exist\n", file);
    return false;
  }
  _type = merylInputType::inSequence;
  _sequenceName = duplicateString(file);
  _squish       = doCompression;
  return true;
}

bool
merylInput::registerSeqStore(char const *path, uint32 seg, uint32 segMax, std::vector<char const *> &err, bool displayErrors) {
  if (isCanuSeqStore(path) == false) {
    if (displayErrors)
      sprintf(err, "'%s' does not appear to be a Canu seqStore\n", path);
    return false;
  }
  _type = merylInputType::inCanu;
  _storeName   = duplicateString(path);
  _storeSeg    = seg;
  _storeSegMax = segMax;
  return true;
}



#ifndef CANU

void
merylInput::openInputSeqStore(void) {}

bool
merylInput::loadBasesFromCanu(char    *seq,
                              uint64   maxLength,
                              uint64  &seqLength,
                              bool    &endOfSequence)  { return(false); }

#endif


void
merylInput::openInput(std::vector<char const *> &err) {

  switch (_type) {
    case merylInputType::inNowhere:
      assert(0);
      break;
    case merylInputType::inTemplate:
      //  Do nothing.
      break;
    case merylInputType::inAction:
      //  Do nothing.
      break;
    case merylInputType::inPipe:
#warning not opening pipe inputs
      break;
    case merylInputType::inDB:
      _db = new merylFileReader(_dbName, false, &err);
      break;
    case merylInputType::inList:
#warning not opening list inputs
      assert(0);
      break;
    case merylInputType::inSequence:
      _sequence = openSequenceFile(_sequenceName);
      break;
    case merylInputType::inCanu:
      openInputSeqStore();
      break;
    default:
      assert(0);
      break;
  }
}



void
merylInput::nextMer(void) {

  if (_db) {
    _valid = _db->nextMer();
    _kmer  = _db->theFMer();
  }

  if (_action) {
    _valid = _action->nextMer();
    _kmer  = _action->theFMer();
  }
}



bool
merylInput::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength,
                      bool    &endOfSequence) {
  bool  gotBases = false;

  seqLength     = 0;
  endOfSequence = true;

  if (_sequence)    gotBases = _sequence->loadBases(seq, maxLength, seqLength, endOfSequence);
  if (_store)       gotBases = loadBasesFromCanu(seq, maxLength, seqLength, endOfSequence);

  //  Homopoly compress if there are bases.
  if ((gotBases) && (_squish))
    seqLength = homopolyCompress(seq, seqLength, seq, NULL, _lastByte);

  //  Save the last byte of the buffer.
  if ((seqLength > 0) && (endOfSequence == false))
    _lastByte = seq[seqLength - 1];
  else
    _lastByte = 0;

  return(gotBases);
}
