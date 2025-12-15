
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

#include "meryl.H"


merylSelector::merylSelector(merylSelectorQuantity type,
                             merylSelectorRelation rela,
                             bool                  invert,
                             char const           *str)   {
  _q   = type;      //  Type and Relation are obvious. The desired
  _r   = rela;      //  result _t is true when 'invert' is false, and
  _t   = !invert;   //  false when 'invert' is true.
  _str = duplicateString(str);
}


merylSelector::merylSelector(const merylSelector &that) {

  *this = that;

  _str = duplicateString(that._str);   //  Replace with our own copy.

  assert(_presentInNum  == nullptr);
  assert(_presentInIdx  == nullptr);
  assert(_presentInList == nullptr);
}



merylSelector::~merylSelector()  {
  delete [] _str;
  delete [] _presentInNum;
  delete [] _presentInIdx;
  delete [] _presentInList;
}



//  Returns true if the selector test is true, meaning the kmer should be
//  output.
//
//  _t is the desired result for a 'true' test here; normally it is 'true',
//  but when 'not' is encountered on the command line, it is set to 'false'.
//
//  This flow control doesn't seem to be a bottleneck.
//   - Entirely hardcoding the test in nextMer() and skipping ALL of the
//     isKmerSelectedOut() code, takes 42u.
//   - Calling isKmerSelectedOut() (implemented as a switch in nextMer.C) and
//    using a set of unrelated if tests here takes 46u.  Using a switch here
//    did not change times.
//

bool
merylSelector::isTrue(kmer const &k, uint32 actLen, merylActList *act, merylActList *inp) const {
  bool   result = _t;


  if (_q == merylSelectorQuantity::isNOP) {
    assert(0);
  }


  if (_q == merylSelectorQuantity::isValue) {
    kmvalu  rhs;
    kmvalu  lhs;

    //  Get the left hand value
    if      (_vIndex1 == uint32max)            lhs = _vValue1;
    else if (_vIndex1 == 0)                    lhs = k._val;
    else if (inp[_vIndex1-1]._idx == 0)        lhs = inp[_vIndex1-1]._val;
    else                                       return(false);

    //  Get the right hand value
    if      (_vIndex2 == uint32max)            rhs = _vValue2;
    else if (_vIndex2 == 0)                    rhs = k._val;
    else if (inp[_vIndex2-1]._idx == 0)        rhs = inp[_vIndex2-1]._val;
    else                                       return(false);

    result = compare(lhs, rhs);
  }


  if (_q == merylSelectorQuantity::isLabel) {
    kmlabl  rhs;
    kmlabl  lhs;

    //  Get the left hand value
    if      (_vIndex1 == uint32max)            lhs = _vLabel1;
    else if (_vIndex1 == 0)                    lhs = k._lab;
    else if (inp[_vIndex1-1]._idx == 0)        lhs = inp[_vIndex1-1]._lab;
    else                                       return(false);

    //  Get the right hand value
    if      (_vIndex2 == uint32max)            rhs = _vLabel2;
    else if (_vIndex2 == 0)                    rhs = k._lab;
    else if (inp[_vIndex2-1]._idx == 0)        rhs = inp[_vIndex2-1]._lab;
    else                                       return(false);

    result = compare(lhs, rhs);
  }


  if (_q == merylSelectorQuantity::isBases) {             //  The two nonsense cases below
    uint32 c = (((_countA == false) ? 0 : countA(k)) +  //  should be caught by isBasesSelector().
                ((_countC == false) ? 0 : countC(k)) +  //  We'll compute the result anyway, so
                ((_countG == false) ? 0 : countG(k)) +  //  if something does change we still
                ((_countT == false) ? 0 : countT(k)));  //  provide some sensible result.

    if      ((_vIndex1 == 0) && (_vIndex2 == 0))    //  Nonsense, just comparing  ! see comment
      result = compare(_vBases1, _vBases2);         //  the two input constants!  ! above

    else if ((_vIndex1 == 0) && (_vIndex2 != 0))    //  Useful: kmer OP constant
      result = compare(c, _vBases2);                //

    else if ((_vIndex1 != 0) && (_vIndex2 == 0))    //  Useful: constant OP kmer
      result = compare(_vBases1, c);                //

    else                                            //  Nonsense, comparing self  ! see comment
      result = compare(c, c);                       //  to self!                  ! above
  }


  if (_q == merylSelectorQuantity::isIndex) {        //  The index selector checks two things:
    result = _presentInNum[actLen];                //    does the kmer appear in the correct
                                                   //    number of inputs (an array lookup).
    for (uint32 pp=0; pp<_presentInLen; pp++) {    //
      uint32  ai = inp[ _presentInList[pp] ]._idx; //    does the kmer appear in ALL specified
                                                   //    inputs (a loop).  presentInList is
      if (ai == uint32max)                         //    the list of inputs that must supply a
        result = false;                            //    kmer, actRdx maps an input to the act[]
    }                                              //    index it is associated with, or uint32max
  }                                                //    if the input didn't supply a kmer.


  //  If the result is the desired result, return 'true'.

  return(result == _t);
}



void
merylSelector::finalizeSelectorInputs(merylOpTemplate *mot, std::vector<char const *> &err) {
  uint32  nInputs = mot->_inputs.size();

  //
  //  Check that any distinct= or word-frequency= value selectors have
  //  exactly one database input.
  //

  if ((_vValue2Distinct >= 0) ||
      (_vValue2WordFreq >= 0)) {

    if (nInputs != 1)
      sprintf(err, "ERROR: selector '%s' invalid; exactly one input required but %u supplied.\n", _str, nInputs);

    if ((nInputs > 0) &&
        (mot->_inputs[0]->_db == nullptr))
      sprintf(err, "ERROR: selector '%s' invalid; database input expected, but type '%s' supplied.\n", _str, mot->_inputs[0]->inputType());

    if ((mot->_valueAssign   != merylAssignValue::valueFirst) ||
        (mot->_valueConstant != 0)) {
      sprintf(err, "ERROR: selector '%s' invalid; the value cannot be modified, remove '%s' modifier\n", _str, toString(mot->_valueAssign));
    }
  }

  //
  //  Build the index: lookup tables.
  //

  _presentInNum = new bool [nInputs + 1];   //  +1 because we actually access [nInputs].
  _presentInIdx = new bool [nInputs];

  _presentInList = new uint32 [nInputs];

  //  Initialize defaults.  If nothing was specified, default to allowing
  //  'any' number of inputs, then set the state of the lookup table to
  //  'true' if 'any' is specified or 'false' otherwise.  Finally, if 'all'
  //  was specified, now that we know the number of inputs, we can set that
  //  to 'true'.

  if ((_input_num.size() == 0) && (_input_num_all == false))
    _input_num_any = true;

  for (uint32 ii=0; ii<=nInputs; ii++)
    _presentInNum[ii] = (_input_num_any == true) ? true : false;

  if (_input_num_all == true)
    _presentInNum[nInputs] = true;

  //  Check each of the entries in _input_num and set the flag for each to true.

  for (uint32 ii=0; ii<_input_num.size(); ii++) {
    uint32  a = _input_num[ii];

    if (a == 0)
      sprintf(err, "selector '%s' invalid; there is no 0th input database.\n", _str);
    else if (a > nInputs)
      sprintf(err, "selector '%s' invalid: cannot occur in %u inputs; there are only %u inputs.\n", _str, a, nInputs);
    else
      _presentInNum[a] = true;
  }

  //  Build a lookup table for the databases that the kmer must be present in.
  //   - initialize everything to false.
  //   - if '@a-all' was supplied, set all those to true.
  //   - set any explicitly specified index to true.

  for (uint32 ii=0; ii<nInputs; ii++)
    _presentInIdx[ii] = false;

  for (uint32 ii=0; ii<_input_idx.size(); ii++) {
    uint32 a = _input_idx[ii];

    if (a == 0)
      sprintf(err, "selector '%s' invalid; there is no 0th input database.\n", _str);
    else if (a > nInputs)
      sprintf(err, "selector '%s' invalid: input %u does not exist; there are only %u inputs.\n", _str, a, nInputs);
    else {
      _presentInIdx[a-1] = true;
      _presentInList[_presentInLen++] = a-1;
    }
  }

  for (uint32 ii=_input_num_at_least; ii<=nInputs; ii++) {
    _presentInIdx[ii-1] = true;
    _presentInList[_presentInLen++] = ii-1;
  }

  std::sort(_presentInList, _presentInList + _presentInLen);
}



//  If either distinct= or word-frequency= were supplied, we need to
//  load the histogram from the input database and compute the
//  threhsold value to use.
//
void
merylSelector::finalizeSelectorParameters(merylOpTemplate *mot) {
  bool  verbose = globals.showConstruction();

  if ((_vValue2Distinct < 0) &&
      (_vValue2WordFreq < 0))
    return;
  
  assert(mot->_inputs.size() == 1);
  assert(mot->_inputs[0]->_db != nullptr);

  merylFileReader         *input = mot->_inputs[0]->_db;
  merylHistogram          *stats = input->stats();   //  Calls input->loadStatistics().
  merylHistogramIterator  *histo = new merylHistogramIterator(stats);

  if (_vValue2Distinct >= 0) {
    uint64  nKmers       = 0;
    uint64  nKmersTarget = _vValue2Distinct * stats->numDistinct();

    if (verbose)
      fprintf(stderr, "finalizeSelectorParameters()-- database '%s' has %lu distinct kmers -> target %lu kmers\n",
              mot->_inputs[0]->_db->filename(), stats->numDistinct(), nKmersTarget);

    for (uint64 ii=0; ii<histo->histogramLength(); ii++) {
      nKmers += histo->histogramOccurrences(ii);

      if (verbose)
        fprintf(stderr, "finalizeSelectorParameters()--   threshold %lu -> %lu cumulative kmers\n",
                histo->histogramValue(ii), nKmers);

      if (nKmers >= nKmersTarget) {
        _vValue2 = histo->histogramValue(ii);
        break;
      }
    }

    if (verbose)
      fprintf(stderr, "finalizeSelectorParameters()-- distinct %f -> threshold %u\n", _vValue2Distinct, _vValue2);
  }

  if (_vValue2WordFreq >= 0) {
    if (verbose)
      fprintf(stderr, "finalizeSelectorParameters()-- database '%s' has %lu total kmers\n",
              mot->_inputs[0]->_db->filename(), stats->numTotal());

    _vValue2 = _vValue2WordFreq * stats->numTotal();

    if (verbose)
      fprintf(stderr, "finalizeSelectorParameters()-- word-frequency %f -> threshold %u\n", _vValue2WordFreq, _vValue2);
  }

  delete histo;

  input->dropStatistics();
}



char *
merylSelector::describeValueSelector(char *str) {
  char  lhs[256] = {0};
  char  rhs[256] = {0};

  if      (_vIndex1 == uint32max) sprintf(lhs, "constant value %s", toDec(_vValue1));
  else if (_vIndex1 == 0)         sprintf(lhs, "output kmer value");
  else                            sprintf(lhs, "kmer value from input %u", _vIndex1);

  if      (_vIndex2 == uint32max) sprintf(rhs, "constant value %s", toDec(_vValue2));
  else if (_vIndex2 == 0)         sprintf(rhs, "output kmer value");
  else                            sprintf(rhs, "kmer value from input %u", _vIndex2);

  sprintf(str, "%s %s %s %s\n", lhs, (_t == false) ? "not" : "   ", toString(_r), rhs);

  return str;
}


char *
merylSelector::describeLabelSelector(char *str) {
  char  lhs[256] = {0};
  char  rhs[256] = {0};

  if      (_vIndex1 == uint32max) sprintf(lhs, "constant label %s", toDec(_vLabel1));
  else if (_vIndex1 == 0)         sprintf(lhs, "output kmer label");
  else                            sprintf(lhs, "kmer label from input %u", _vIndex1);

  if      (_vIndex2 == uint32max) sprintf(rhs, "constant label %s", toDec(_vLabel2));
  else if (_vIndex2 == 0)         sprintf(rhs, "output kmer label");
  else                            sprintf(rhs, "kmer label from input %u", _vIndex2);

  sprintf(str, "%s %s %s %s\n", lhs, (_t == false) ? "not" : "   ", toString(_r), rhs);

  return str;
}


char *
merylSelector::describeBasesSelector(char *str) {
  char  lhs[256] = {0};
  char  rhs[256] = {0};
  char  acgt[5] = {0};
  int32 acgtlen = 0;

  if (_countA == true)   acgt[acgtlen++] = 'A';
  if (_countC == true)   acgt[acgtlen++] = 'C';
  if (_countT == true)   acgt[acgtlen++] = 'T';
  if (_countG == true)   acgt[acgtlen++] = 'G';

  //  vIndex is 0 if nothing was supplied to the selector, which tells us
  //  to get the value from the output kmer.

  if (_vIndex1 == 0) sprintf(lhs, "number of %s in the kmer", acgt);
  else               sprintf(lhs, "%s", toDec(_vBases1));

  if (_vIndex2 == 0) sprintf(rhs, "number of %s in the kmer", acgt);
  else               sprintf(rhs, "%s", toDec(_vBases2));

  sprintf(str, "%s %s %s %s\n", lhs, (_t == false) ? "not" : "   ", toString(_r), rhs);

  return str;
}


char *
merylSelector::describeIndexSelector(char *str) {
  char  lhs[256] = {0};
  char  rhs[256] = {0};

  sprintf(str, "<index selector not described>\n");

  return str;
}


char *
merylSelector::describe(char *str) {

  //  Lots of duplication with isTrue().

  switch (_q) {
    case merylSelectorQuantity::isValue:   describeValueSelector(str);    break;
    case merylSelectorQuantity::isLabel:   describeLabelSelector(str);    break;
    case merylSelectorQuantity::isBases:   describeBasesSelector(str);    break;
    case merylSelectorQuantity::isIndex:   describeIndexSelector(str);    break;
    default:
      sprintf(str, "<empty selector>\n");
      break;
  }

  return(str);
}
