
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

#include <memory>

merylOpCompute::merylOpCompute(merylOpTemplate *ot, uint32 dbSlice, uint32 nInputs) {
  _ot            = ot;

  //  Allocate space for the active list, or use our built in space.  This
  //  _might_ be solving a performance bottleneck, though benchmarks seem to
  //  show no change between using _acts and heap allocated space.

  if (nInputs > merylActListMax) {
    _actAlloc = new merylActList [nInputs * 2];
    _acta     = _actAlloc + 0;
    _inpa     = _actAlloc + nInputs;
  }
  else {
    _acta = _acts;
    _inpa = _inps;
  }

  //  Initialize so that each input is listed as active.  This forces
  //  merylOpCompute::nextMer() to call nextMer() on the first iteration.

  for (uint32 ii=0; ii<nInputs; ii++) {
    _acta[_actLen]._idx = ii;
    _acta[_actLen]._val = 0;
    _acta[_actLen]._lab = 0;

    _inpa[_actLen]._idx = uint32max;
    _inpa[_actLen]._val = 0;
    _inpa[_actLen]._lab = 0;

    _actLen++;
  }
}



merylOpCompute::~merylOpCompute() {

  for (uint32 ii=0; ii<_inputs.size(); ii++)   //  Destroy the inputs, which will recursively destroy
    delete _inputs[ii];                        //  any merylOpCompute objects we have as inputs.

  _inputs.clear();

  delete    _statsAcc;                         //  Close per-slice stats accumulator.
  delete    _outDbseSlice;                     //  Close per-slice output.
  delete    _outListSlice;                     //  Close per-slice list.
}



////////////////////////////////////////
//
//  COMPUTING the kmer/value/label to output.
//
//   - findOutputKmer() scans all the inputs to find the smallest kmer, then
//     creates a list of the inputs with that kmer.
//
//     _kmer._mer is not valid if _actLen is zero after this function.
//
//   - findOutputValue() and findOutputLabel() are reasonably straight
//     forward, but long, and compute an output value/label based on the
//     action specified.
//
void
merylOpCompute::findOutputKmer(void) {

  _actLen = 0;

  //  This sets:
  //    _kmer._mer to the smallest kmer in the input
  //    _inpa[]    to the value/label of each input
  //    _acta[]    to the value/label of each input with the smallest kmer
  //
  //    _inpa[]._idx is 0 if the kmer is the smallest kmer
  
  for (uint32 ii=0; ii<_inputs.size(); ii++) {
    kmdata kmer = _inputs[ii]->_kmer._mer;

    if (_inputs[ii]->_valid == false)      //  No more kmers in the file,
      continue;                            //  skip it.

    _inpa[ii]._idx = uint32max;
    _inpa[ii]._val = _inputs[ii]->_kmer._val;
    _inpa[ii]._lab = _inputs[ii]->_kmer._lab;

    if ((_actLen > 0) &&                   //  If we've picked a kmer already,
        (kmer > _kmer._mer))               //  and this one is bigger,
      continue;                            //  skip this one.

    if ((_actLen > 0) &&                   //  If we've picked a kmer already,
        (kmer < _kmer._mer))               //  but this one is smaller,
      _actLen = 0;                         //  forget everything we've done.

    if (_actLen == 0)                      //  Pick this one if nothing picked yet.
      _kmer._mer = kmer;

    _acta[_actLen]._idx = ii;
    _acta[_actLen]._val = _inputs[ii]->_kmer._val;
    _acta[_actLen]._lab = _inputs[ii]->_kmer._lab;

    _actLen++;
  }

  //  If we need to use _inpa[] do another pass to set _idx correctly.

  for (uint32 ii=0; ii<_actLen; ii++) {
    uint32  idx = _acta[ii]._idx;

    _inpa[idx]._idx = 0;
  }
}



void
merylOpCompute::findOutputValue(void) {
  kmvalu  q = 0;
  kmvalu  r = 0;

  switch (_ot->_valueAssign) {
    case merylAssignValue::valueNOP:
      break;

    case merylAssignValue::valueSet:
      _kmer._val = _ot->_valueConstant;
      break;

    case merylAssignValue::valueSelected:
#warning wrong - need to figure out which input to select
      _kmer._val = _acta[0]._val;
      break;

    case merylAssignValue::valueFirst:
#warning wrong - do we need to verify that actIdx[0] is 0?
      _kmer._val = _acta[0]._val;
      break;

    case merylAssignValue::valueMin:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        _kmer._val = std::min(_kmer._val, _acta[ii]._val);
      break;

    case merylAssignValue::valueMax:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        _kmer._val = std::max(_kmer._val, _acta[ii]._val);
      break;

    case merylAssignValue::valueAdd:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        if (kmvalumax - _kmer._val < _acta[ii]._val)
          _kmer._val = kmvalumax;
        else
          _kmer._val = _kmer._val + _acta[ii]._val;
      break;

    case merylAssignValue::valueSub:
      _kmer._val = _acta[0]._val;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_kmer._val > _acta[ii]._val)
          _kmer._val -= _acta[ii]._val;
        else
          _kmer._val = 0;

      if (_kmer._val > _ot->_valueConstant)
        _kmer._val -= _ot->_valueConstant;
      else
        _kmer._val = 0;

      break;

    case merylAssignValue::valueMul:
      _kmer._val = _ot->_valueConstant;
      for (uint32 ii=0; ii<_actLen; ii++)
        if (kmvalumax / _kmer._val < _acta[ii]._val)
          _kmer._val = kmvalumax;
        else
          _kmer._val = _kmer._val * _acta[ii]._val;
      break;


    case merylAssignValue::valueDiv:
      _kmer._val = _acta[0]._val;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_acta[ii]._val > 0)
          _kmer._val /= _acta[ii]._val;
        else
          _kmer._val = 0;

      if (_ot->_valueConstant > 0)
        _kmer._val /= _ot->_valueConstant;
      else
        _kmer._val = 0;

      break;


      //  Division, but with rounding instead of truncation.
      //  Additionally, values between 0 and 0.5 are rounded up to 1.
      //
      //  However, division by zero results in a 0 output.
    case merylAssignValue::valueDivZ:
      _kmer._val = _acta[0]._val;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_acta[ii]._val == 0)
          _kmer._val = 0;
        else if (_kmer._val < _acta[ii]._val)
          _kmer._val = 1;
        else
          _kmer._val = round(_kmer._val / (double)_acta[ii]._val);

      if (_ot->_valueConstant == 0)
        _kmer._val = 0;
      else if (_kmer._val < _ot->_valueConstant)
        _kmer._val = 1;
      else
        _kmer._val = round(_kmer._val / (double)_ot->_valueConstant);

      break;


    case merylAssignValue::valueMod:
      q = _acta[0]._val;
      r = 0;

      for (uint32 ii=1; ii<_actLen; ii++)
        if (_acta[ii]._val > 0) {
          kmvalu  qt = q / _acta[ii]._val;

          r += q - qt * _acta[ii]._val;
          q  = qt;
        } else {
          r += q;
          q  = 0;
        }

      if (_ot->_valueConstant > 0) {
        kmvalu  qt = q / _ot->_valueConstant;

        r += q - qt * _ot->_valueConstant;
        q  = qt;
      } else {
        r += q;
        q  = 0;
      }

      _kmer._val = r;

      break;


    case merylAssignValue::valueCount:
      _kmer._val = _actLen;
      break;
  }
}



void
merylOpCompute::findOutputLabel(void) {
  kmlabl  l;
  kmvalu  v;

  switch (_ot->_labelAssign) {
    case merylAssignLabel::labelNOP:
      break;

    case merylAssignLabel::labelSet:
      _kmer._lab = _ot->_labelConstant;
      break;

    case merylAssignLabel::labelSelected:
#warning wrong - need to figure out which input to select
      _kmer._lab = _acta[0]._lab;
      break;

    case merylAssignLabel::labelFirst:
#warning wrong - do we need to verify that actIdx[0] is 0?
      _kmer._lab = _acta[0]._lab;
      break;

    case merylAssignLabel::labelMin:
      l = _ot->_labelConstant;
      v = kmvalumax;

      for (uint32 ll=0; ll<_actLen; ll++)
        if (_acta[ll]._val < v) {
          l = _acta[ll]._lab;
          v = _acta[ll]._val;
        }

      _kmer._lab = l;
      break;

    case merylAssignLabel::labelMax:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab = std::max(_kmer._lab, _acta[ll]._lab);
      break;


    case merylAssignLabel::labelAnd:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab &= _acta[ll]._lab;
      break;

    case merylAssignLabel::labelOr:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab |= _acta[ll]._lab;
      break;

    case merylAssignLabel::labelXor:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        _kmer._lab ^= _acta[ll]._lab;
      break;

    case merylAssignLabel::labelDifference:
#warning check this
      _kmer._lab = _acta[0]._lab & ~_ot->_labelConstant;
      for (uint32 ll=1; ll<_actLen; ll++)
        _kmer._lab &= ~_acta[ll]._lab;
      break;

    case merylAssignLabel::labelLightest:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        if (countNumberOfSetBits64(_acta[ll]._lab) < countNumberOfSetBits64(_kmer._lab))
          _kmer._lab = _acta[ll]._lab;
      break;

    case merylAssignLabel::labelHeaviest:
      _kmer._lab = _ot->_labelConstant;
      for (uint32 ll=0; ll<_actLen; ll++)
        if (countNumberOfSetBits64(_acta[ll]._lab) > countNumberOfSetBits64(_kmer._lab))
          _kmer._lab = _acta[ll]._lab;
      break;

    case merylAssignLabel::labelInvert:
      assert(_actLen == 1);
      _kmer._lab = ~_acta[0]._lab;
      break;

    case merylAssignLabel::labelShiftLeft:
      assert(_actLen == 1);
      _kmer._lab = _acta[0]._lab >> _ot->_labelConstant;
      break;

    case merylAssignLabel::labelShiftRight:
      assert(_actLen == 1);
      _kmer._lab = _acta[0]._lab << _ot->_labelConstant;
      break;

#warning wrong
    case merylAssignLabel::labelRotateLeft:
      assert(_actLen == 1);
      _kmer._lab = _acta[0]._lab >> _ot->_labelConstant;
      break;

#warning wrong
    case merylAssignLabel::labelRotateRight:
      assert(_actLen == 1);
      _kmer._lab = _acta[0]._lab << _ot->_labelConstant;
      break;
  }
}
