
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
#include <cmath>



void
printKmer(FILE *f, kmer pk, bool printACGTorder) {
  char   outstr[256];
  char  *outptr = outstr;

  if (printACGTorder == true)              //  If requested, recompute the canonical
    pk.recanonicalizeACGTorder();          //  mer in ACGT order.  Yuck.

  pk.toString(outptr);                     //  Convert the kmer to ASCII, then
  while (*outptr)                          //  advance to the end of the string.
    outptr++;

  *outptr++ = '\t';                        //  Add the value.  There is always
  outptr = toDec(pk._val, outptr);         //  a value to add.

#warning need to print label as binary or hex, user supplied
  if (kmer::labelSize() > 0) {             //  If a label exists, add it too.
    *outptr++ = '\t';
    outptr = toBin(pk._lab, outptr, kmer::labelSize());
  }

  *outptr++ = '\n';                        //  Terminate the string and
  *outptr++ = 0;                           //  emit it.

#pragma omp critical (printLock)           //  fputs() is not thread safe and will
  fputs(outstr, f);                        //  happily intermix on e.g. Linux.
}



//  Return true if a single selector product term (_select[ii]) is true.  A
//  single term is true if all of it's terms (_select[ii][0], _select[ii][1],
//  ...) are true.
//
inline
bool
merylOpCompute::shouldKmerBeOutput(uint32 ii) {

#if 1

  switch (_ot->_select[ii].size()) {
    case 0:
      return(true);
      break;
    case 1:
      return((_ot->_select[ii][0].isTrue(_kmer, _actLen, _acta, _inpa)));
      break;
    case 2:
      return((_ot->_select[ii][0].isTrue(_kmer, _actLen, _acta, _inpa)) &&
             (_ot->_select[ii][1].isTrue(_kmer, _actLen, _acta, _inpa)));
      break;
    case 3:
      return((_ot->_select[ii][0].isTrue(_kmer, _actLen, _acta, _inpa)) &&
             (_ot->_select[ii][1].isTrue(_kmer, _actLen, _acta, _inpa)) &&
             (_ot->_select[ii][2].isTrue(_kmer, _actLen, _acta, _inpa)));
      break;
    case 4:
      return((_ot->_select[ii][0].isTrue(_kmer, _actLen, _acta, _inpa)) &&
             (_ot->_select[ii][1].isTrue(_kmer, _actLen, _acta, _inpa)) &&
             (_ot->_select[ii][2].isTrue(_kmer, _actLen, _acta, _inpa)) &&
             (_ot->_select[ii][3].isTrue(_kmer, _actLen, _acta, _inpa)));
      break;
    default:
      for (uint32 tt=0; tt<_ot->_select[ii].size(); tt++)
        if (_ot->_select[ii][tt].isTrue(_kmer, _actLen, _acta, _inpa) == false)
          return(false);
      return(true);
      break;
  }

#else

  for (uint32 tt=0; tt<_ot->_select[ii].size(); tt++)
    if (_ot->_select[ii][tt].isTrue(_kmer, _actLen, _acta, _inpa) == false)
      return(false);
  return(true);

#endif
}



//  Returns true if the kmer should be output, based on the selectors.
//
//  This function is true of ANY of the selector product terms (_select[0],
//  _select[1], ...) are true.
//
inline
bool
merylOpCompute::shouldKmerBeOutput(void) {

  //  Simple end cases (note that order is important):
  //    the kmer IS     selected out if the value is zero
  //    the kmer is NOT selected out if there are no selectors
  //
  if (_kmer._val == 0)            return(false);
  if (_ot->_select.size() == 0)   return(true);

  //  Test each selector.  A selector will return true if the kmer should
  //  be output; false if the kmer should be skipped.
  //
  //  A kmer should be output if ANY selector[ii] is true.
  //
  //  A single selector[ii] is true if ALL of it's pieces are true.

#if 1

  //  The switches are significantly faster than a for loop!
  //
  //  Looking at user time for a run with three product terms, the first two
  //  with one selector, the last with N selectors:
  //
  //    SWITCH  FOR-LOOP  DIFF        time
  //    ------  --------  ----        /work/meryl-redo/build/bin/meryl -Q
  //     55.22     59.44  4.22          label:eq:9
  //     58.58     62.14  3.56            or
  //     62.73     66.25  3.52          label:eq:9
  //     67.59     71.17  3.58            or
  //     72.99     74.77  1.78          value:ge:0 value:ge:0 value:ge:100000
  //     77.57     78.31  0.74          fArcCen1.k21.gt1.meryl
  //     81.55     82.36  0.81

  switch (_ot->_select.size()) {
    case 1:
      return((shouldKmerBeOutput(0) == true));
      break;
    case 2:
      return((shouldKmerBeOutput(0) == true) ||
             (shouldKmerBeOutput(1) == true));
      break;
    case 3:
      return((shouldKmerBeOutput(0) == true) ||
             (shouldKmerBeOutput(1) == true) ||
             (shouldKmerBeOutput(2) == true));
      break;
    case 4:
      return((shouldKmerBeOutput(0) == true) ||
             (shouldKmerBeOutput(1) == true) ||
             (shouldKmerBeOutput(2) == true) ||
             (shouldKmerBeOutput(3) == true));
      break;
    default:
      for (uint32 ii=0; ii<_ot->_select.size(); ii++)
        if (shouldKmerBeOutput(ii) == true)
          return(true);
      return(false);
      break;
  }

  assert(0);
  return(true);

#else

  bool   r = false;

  for (uint32 ii=0; ii<_ot->_select.size(); ii++) {
    bool t = true;

    for (uint32 tt=0; tt<_ot->_select[ii].size(); tt++)
      t &= _ot->_select[ii][tt].isTrue(_kmer, _actLen, _acta, _inpa);

    r |= t;
  }

  return(r);

#endif
}



bool
merylOpCompute::nextMer(void) {
  bool   isEmpty = false;

 nextMerAgain:

  //  Grab the next mer for every input that was active in the last
  //  iteration.  (on the first call, all inputs were 'active' last time)

  for (uint32 ii=0; ii<_actLen; ii++)
    _inputs[_acta[ii]._idx]->nextMer();

  //  Find the smallest kmer in any input, and remember the values and labels
  //  of the kmer in each input file.

  findOutputKmer();

  //  If no active kmers, we're done.  There are no kmers left in any input.

  if (_actLen == 0)
    return(false);

  //  Figure out the value/label of the output kmer.

  findOutputValue();
  findOutputLabel();

  //  If this kmer is selected, go back and get another kmer from the inputs.

  if (shouldKmerBeOutput() == false)
    goto nextMerAgain;

  //  Analyze the kmer.

  if (_statsAcc)
    _statsAcc->addValue(_kmer._val);

  //  Output the kmer to whatever outputs exist.

  if (_outDbse != nullptr)
    _outDbse->addMer(_kmer);

  if (_outList != nullptr)
    printKmer(_outList->file(), _kmer, _ot->_printACGTorder);

  if (_outShow != nullptr)
    printKmer(_outShow->file(), _kmer, _ot->_printACGTorder);

  //  We have now loaded the nextMer; return and let our clients query us.

  return(true);
}
