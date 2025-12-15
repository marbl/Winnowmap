
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

//
//  Handle value and label assignment methods.
//
//  Unlike 'output', the _curParam MUST be valid in all cases.  The word must
//  be of form 'value=something' or 'label=something'.  It would be trivial
//  to allow 'value= something'....it just looks weird.
//

void
merylOpTemplate::setValueConstant(char const *s, merylAssignValue a, kmvalu v) {
  _valueString   = duplicateString(s);
  _valueAssign   = a;
  _valueConstant = v;
}

void
merylOpTemplate::setLabelConstant(char const *s, merylAssignLabel a, kmlabl l) {
  _labelString   = duplicateString(s);
  _labelAssign   = a;
  _labelConstant = l;
}


bool
merylCommandBuilder::isAssignValue(void) {
  merylOpTemplate  *op = getCurrent();
  kmvalu            cs;

  if (_curPname != opPname::pnValue)
    return false;

  fprintf(stderr, "isAssignValue()-- '%s' with param '%s'\n", _optString, _curParam);

  //  Doesn't work; default value is not valueNOP!
  //if (op->_valueAssign != merylAssignValue::valueNOP) {
  //  sprintf(_errors, "Assignment '%s' would replace existing %s assignment.", _optString, toString(op->_valueAssign));
  //  return true;
  //}

  //  Check for modifiers with constants at the end.

  if      (strncmp(_curParam, "#", 1)        == 0)   op->setValueConstant(_curParam+1, merylAssignValue::valueSet,      decodeInteger(_curParam, 1, 0, cs, _errors));
  else if (strncmp(_curParam, "min#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueMin,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "max#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueMax,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "add#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueAdd,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "sum#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueAdd,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "sub#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueSub,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "dif#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueSub,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "mul#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueMul,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "div#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueDiv,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "divzero#", 8) == 0)   op->setValueConstant(_curParam+8, merylAssignValue::valueDivZ,     decodeInteger(_curParam, 8, 0, cs, _errors));
  else if (strncmp(_curParam, "mod#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueMod,      decodeInteger(_curParam, 4, 0, cs, _errors));
  else if (strncmp(_curParam, "rem#", 4)     == 0)   op->setValueConstant(_curParam+4, merylAssignValue::valueMod,      decodeInteger(_curParam, 4, 0, cs, _errors));

  //  Check for modifiers without constants.  Set the constant to whatever
  //  the identity is for the given modifier.

  else if (strncmp(_curParam, "first", 6)    == 0)   op->setValueConstant(_curParam, merylAssignValue::valueFirst,    0);
  else if (strncmp(_curParam, "selected", 9) == 0)   op->setValueConstant(_curParam, merylAssignValue::valueSelected, kmvalumax);
  else if (strncmp(_curParam, "min", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueMin,      kmvalumax);
  else if (strncmp(_curParam, "max", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueMax,      0);
  else if (strncmp(_curParam, "add", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueAdd,      0);
  else if (strncmp(_curParam, "sum", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueAdd,      0);
  else if (strncmp(_curParam, "sub", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueSub,      0);
  else if (strncmp(_curParam, "dif", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueSub,      0);
  else if (strncmp(_curParam, "mul", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueMul,      1);
  else if (strncmp(_curParam, "div", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueDiv,      1);
  else if (strncmp(_curParam, "divzero", 8)  == 0)   op->setValueConstant(_curParam, merylAssignValue::valueDivZ,     1);
  else if (strncmp(_curParam, "mod", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueMod,      0);
  else if (strncmp(_curParam, "rem", 4)      == 0)   op->setValueConstant(_curParam, merylAssignValue::valueMod,      0);
  else if (strncmp(_curParam, "count", 6)    == 0)   op->setValueConstant(_curParam, merylAssignValue::valueCount,    0);

  //  Nope, don't know what this is.

  else {
    sprintf(_errors, "Unknown assign:value=<parameter> in '%s'.", _optString);
    return false;
  }

  resetClass();

  return true;
}


bool
merylCommandBuilder::isAssignLabel(void) {
  merylOpTemplate  *op = getCurrent();
  kmlabl            cs;

  if (_curPname != opPname::pnLabel)
    return false;

  fprintf(stderr, "isAssignLabel()-- '%s' with param '%s'\n", _optString, _curParam);

  //  Doesn't work; default value is not labelNOP!
  //if (op->_labelAssign != merylAssignLabel::labelNOP) {
  //  sprintf(_errors, "Assignment '%s' would replace existing %s assignment.", _optString, toString(op->_labelAssign));
  //  return true;
  //}

  //  Check for modifiers with constants at the end.

  if      (strncmp(_curParam, "#", 1)              == 0)   op->setLabelConstant(_curParam+1,  merylAssignLabel::labelSet,         decodeInteger(_curParam, 1,  0, cs, _errors));
  else if (strncmp(_curParam, "min#", 4)           == 0)   op->setLabelConstant(_curParam+4,  merylAssignLabel::labelMin,         decodeInteger(_curParam, 4,  0, cs, _errors));
  else if (strncmp(_curParam, "max#", 4)           == 0)   op->setLabelConstant(_curParam+4,  merylAssignLabel::labelMax,         decodeInteger(_curParam, 4,  0, cs, _errors));
  else if (strncmp(_curParam, "and#", 4)           == 0)   op->setLabelConstant(_curParam+4,  merylAssignLabel::labelAnd,         decodeInteger(_curParam, 4,  0, cs, _errors));
  else if (strncmp(_curParam, "or#", 3)            == 0)   op->setLabelConstant(_curParam+3,  merylAssignLabel::labelOr,          decodeInteger(_curParam, 3,  0, cs, _errors));
  else if (strncmp(_curParam, "xor#", 4)           == 0)   op->setLabelConstant(_curParam+4,  merylAssignLabel::labelXor,         decodeInteger(_curParam, 4,  0, cs, _errors));
  else if (strncmp(_curParam, "difference#", 11)   == 0)   op->setLabelConstant(_curParam+11, merylAssignLabel::labelDifference,  decodeInteger(_curParam, 11, 0, cs, _errors));
  else if (strncmp(_curParam, "lightest#", 9)      == 0)   op->setLabelConstant(_curParam+9,  merylAssignLabel::labelLightest,    decodeInteger(_curParam, 9,  0, cs, _errors));
  else if (strncmp(_curParam, "heaviest#", 9)      == 0)   op->setLabelConstant(_curParam+9,  merylAssignLabel::labelHeaviest,    decodeInteger(_curParam, 9,  0, cs, _errors));
  else if (strncmp(_curParam, "invert#", 7)        == 0)   op->setLabelConstant(_curParam+7,  merylAssignLabel::labelInvert,      decodeInteger(_curParam, 7,  0, cs, _errors));
  else if (strncmp(_curParam, "shift-left#", 11)   == 0)   op->setLabelConstant(_curParam+11, merylAssignLabel::labelShiftLeft,   decodeInteger(_curParam, 11, 0, cs, _errors));
  else if (strncmp(_curParam, "shift-right#", 12)  == 0)   op->setLabelConstant(_curParam+12, merylAssignLabel::labelShiftRight,  decodeInteger(_curParam, 12, 0, cs, _errors));
  else if (strncmp(_curParam, "rotate-left#", 12)  == 0)   op->setLabelConstant(_curParam+12, merylAssignLabel::labelRotateLeft,  decodeInteger(_curParam, 12, 0, cs, _errors));
  else if (strncmp(_curParam, "rotate-right#", 13) == 0)   op->setLabelConstant(_curParam+13, merylAssignLabel::labelRotateRight, decodeInteger(_curParam, 13, 0, cs, _errors));

  //  Check for modifiers without constants.  Set the constant to whatever
  //  the identity is for the given modifier.

  else if (strncmp(_curParam, "first", 6)         == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelFirst,       0);
  else if (strncmp(_curParam, "selected", 9)      == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelSelected,    0);
  else if (strncmp(_curParam, "min", 4)           == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelMin,         0);
  else if (strncmp(_curParam, "max", 4)           == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelMax,         0);
  else if (strncmp(_curParam, "and", 4)           == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelAnd,         kmlablmax);
  else if (strncmp(_curParam, "or", 3)            == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelOr,          0);
  else if (strncmp(_curParam, "xor", 4)           == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelXor,         kmlablmax);
  else if (strncmp(_curParam, "difference", 11)   == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelDifference,  0);
  else if (strncmp(_curParam, "lightest", 9)      == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelLightest,    kmlablmax);
  else if (strncmp(_curParam, "heaviest", 9)      == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelHeaviest,    0);
  else if (strncmp(_curParam, "invert", 7)        == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelInvert,      0);
  else if (strncmp(_curParam, "shift-left", 11)   == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelShiftLeft,   1);
  else if (strncmp(_curParam, "shift-right", 12)  == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelShiftRight,  1);
  else if (strncmp(_curParam, "rotate-left", 12)  == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelRotateLeft,  1);
  else if (strncmp(_curParam, "rotate-right", 13) == 0)   op->setLabelConstant(_curParam, merylAssignLabel::labelRotateRight, 1);

  //  Nope, don't know what this is.

  else {
    sprintf(_errors, "Unknown assign:label=<parameter> in '%s'.", _optString);
    return false;
  }

  resetClass();

  return true;
}



bool
merylCommandBuilder::isAssign(void) {

  if (_curClass != opClass::clAssign)      //  Silently ignore non-assign classes.
    return false;

  if ((_curPname != opPname::pnValue) &&   //  Noisly complain about unknown Pnames.
      (_curPname != opPname::pnLabel)) {
    sprintf(_errors, "expecting value or label parameter name in '%s'\n", _optString);
    return false;
  }

  if ((_curParam == nullptr) ||            //  Noisly complain about missing parameters.
      (_curParam[0] == 0)) {
    sprintf(_errors, "expecting parameters in '%s'\n", _optString);
    return false;
  }

  return ((isAssignValue() == true) ||     //  Process the parameter.
          (isAssignLabel() == true));
}
