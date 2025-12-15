
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

char const *begin   = "({c}";     //  Words start with a capture group.
char const *end     = ")[]\v]";   //  And end with a word-separator or an action close operator.

char const *col     = "\\s*:\\s*";
char const *equ     = "\\s*=\\s*";

char const *base2   = "({c}[+-]?[0-1]+[bB])";
char const *base8   = "({c}[+-]?[0-7]+[oO])";
char const *baseDd  = "({c}[+-]?[0-9]+[dD]?)";
char const *baseDs  = "({c}[+-]?[0-9]+[kKmMgGtTpPeE][iI]?)";
char const *baseH   = "({c}[+-]?[0-9a-fA-F]+[hH])";
char const *baseF   = "({c}[+-]?[0-9]*[.]?[0-9]+([eE][+-]?[0-9]+)?)";

//  baseF doesn't allow '0.'

char const *integer = "({c}([+-]?[01]+[bB])|"
                          "([+-]?[0-7]+[oO])|"
                          "([+-]?[0-9]+[dD]?)|"
                          "([+-]?[0-9]+[kKmMgGtTpPeE][iI]?)|"
                          "([+-]?[0-9a-fA-F]+[hH])"
                      ")";
char const *number  = "({c}([+-]?[01]+[bB])|"
                          "([+-]?[0-7]+[oO])|"
                          "([+-]?[0-9]+[dD]?)|"
                          "([+-]?[0-9]+[kKmMgGtTpPeE][iI]?)|"
                          "([+-]?[0-9a-fA-F]+[hH])|"
                          "([+-]?[0-9]*[.]?[0-9]+([eE][+-]?[0-9]+)?)"
                      ")";

//  Pathnames allow <SPACE>, <TAB> and all the printable characters, except
//  the last character cannot be a ']' because it interferes with our
//  detection of action-close operators.
//
char const *path    = "({c}[ \t[:print:]]*[ \t[:print:]~]])";

void
merylCommandBuilder::compileRegExes(void) {

  _regex[0x00].compile("({c}\v+)", nullptr);  //  Match word-seps, MUST be an extra capture group to fit with the rest.
  _regex[0x01].compile("({c}\\[)", nullptr);  //  Match a single '[' to start a new operation.
  _regex[0x02].compile("({c}\\])", nullptr);  //  Match a single ']' to complete an operation.

  //egex[0x03];
  //egex[0x04];

  //  COUNTING and COUNT OPTIONS
  //   - count[-forward | -reverse]
  //   - expect=<int>
  //   - suffix=<ACGTacgt>
  //   - segment=<int>/<int>
  //  

  //_regex[0x05].enableVerbose();
  _regex[0x05].compile(begin, "count(-({c}({p}forward)|({p}reverse)))?",  end, nullptr);
  _regex[0x06].compile(begin, "(({c}expect)", equ, integer, ")|(({c}suffix)", equ, "({c}[ACGTacgt]+))|(({c}segment)", equ, integer, "/", integer, ")", end, nullptr);


  //egex[0x07];
  //egex[0x08];

  //  COMBINING ALIASES
  //   - subtract,   union,     union-min,     union-max,     union-sum
  //   - difference, intersect, intersect-min, intersect-max, intersect-sum
  //
  _regex[0x09].compile(begin, "union"                 "(-({c}min|max|sum))?", end, nullptr);
  _regex[0x0a].compile(begin, "(isect|in({p}tersect))""(-({c}min|max|sum))?", end, nullptr);

  _regex[0x0b].compile(begin, "sub({p}tract)?",          end, nullptr);
  _regex[0x0c].compile(begin, "diff({p}erence)?",        end, nullptr);

  //egex[0x0d];
  //egex[0x0e];
  //egex[0x0f];

  //  FILTERING ALIASES
  //   - less-than     <tdw>
  //   - greater-than  <tdw>
  //   - at-least      <tdw>
  //   - at-most       <tdw>
  //   - equal-to      <tdw>
  //   - not-equal-to  <tdw>
  //
  //  tdw matches:
  //   - <integer>
  //   - threshold=<integer>
  //   - distinct=<float>
  //   - word-frequency=<float>  (and wf=<float>)
  //
  char *tdw = new char [128 + 2*strlen(integer) + 3*strlen(equ) + 2*strlen(baseF)];
  sprintf(tdw, "(%s|({c,p}threshold)%s%s|({c,p}distinct)%s%s|({c,p}w[of]rd-frequency)%s%s)",
          integer, equ, integer, equ, baseF, equ, baseF);

  _regex[0x10].compile(begin, "less-than\\s+",    tdw, end, nullptr);
  _regex[0x11].compile(begin, "greater-than\\s+", tdw, end, nullptr);
  _regex[0x12].compile(begin, "at-least\\s+",     tdw, end, nullptr);
  _regex[0x13].compile(begin, "at-most\\s+",      tdw, end, nullptr);
  _regex[0x14].compile(begin, "equal-to\\s+",     tdw, end, nullptr);
  _regex[0x15].compile(begin, "not-equal-to\\s+", tdw, end, nullptr);

  delete [] tdw;

  //egex[0x16];
  //egex[0x17];

  //  MODIFYING ALIASES
  //   - increase      <integer>
  //   - decrease      <integer>
  //   - multiply      <integer> - Floats make a little bit of sense
  //   - divide        <integer> - here, but are not supported.
  //   - divide-round  <integer> - 
  //   - modulo        <integer>
  //
  _regex[0x18].compile(begin, "increase\\s+",     integer, end, nullptr);
  _regex[0x19].compile(begin, "decrease\\s+",     integer, end, nullptr);
  _regex[0x1a].compile(begin, "multiply\\s+",     integer, end, nullptr);
  _regex[0x1b].compile(begin, "divide\\s+",       integer, end, nullptr);
  _regex[0x1c].compile(begin, "divide-round\\s+", integer, end, nullptr);
  _regex[0x1d].compile(begin, "modulo\\s+",       integer, end, nullptr);

  //egex[0x1e];

  //  Output parameters.
  //
  //  Compatibility mode 'output <path>' must be last so it doesn't match to 'output : db = file'.

  _regex[0x1f].compile(begin, "({p}output)", col, "s", end, nullptr);  //  Ambiguous.

  _regex[0x20].compile(begin, "({p}output)", col, "({p}d[ab]tabase)",        equ, path,       end, nullptr);
  _regex[0x21].compile(begin, "({p}output)", col, "({p}list)",               equ, path,       end, nullptr);
  _regex[0x22].compile(begin, "({p}output)", col, "({p}show)",                                end, nullptr);
  _regex[0x23].compile(begin, "({p}output)", col, "({p}pipe)",               equ, path,       end, nullptr);
  //_regex[0x24].enableVerbose();
  _regex[0x24].compile(begin, "({p}output)", col, "({p}histogram)",     "(", equ, path, ")?", end, nullptr);
  _regex[0x25].compile(begin, "({p}output)", col, "({p}stat[is]stics)", "(", equ, path, ")?", end, nullptr);

  _regex[0x26].compile(begin,    "(output)\\s",                                 path, end, nullptr);
  //egex[0x27];

  //  Input parameters.

  _regex[0x28].compile(begin, "({p}input)", col, "({p}database)", equ, path, end, nullptr);
  _regex[0x29].compile(begin, "({p}input)", col, "({p}list)",     equ, path, end, nullptr);
  _regex[0x2a].compile(begin, "({p}input)", col, "({p}pipe)",     equ, path, end, nullptr);
  _regex[0x2b].compile(begin, "({p}input)", col, "({p}action)",   equ, path, end, nullptr);

  //egex[0x2c];
  //egex[0x2d];
  //egex[0x2e];
  //egex[0x2f];

  //  Assign.

  _regex[0x30].compile(begin, "(({p}assign)|set)", col, "({p}value)" "", equ, "([^\v]*)", end, nullptr);
  _regex[0x40].compile(begin, "(({p}assign)|set)", col, "({p}label)" "", equ, "([^\v]*)", end, nullptr);

  //  Select.
  //   - 'and' 'or' 'not'
  //
  _regex[0x50].compile(begin, "(({p}select)|get)", col, "({p}value)" "", col, "([^\v]*)", end, nullptr);
  _regex[0x60].compile(begin, "(({p}select)|get)", col, "({p}label)" "", col, "([^\v]*)", end, nullptr);
  _regex[0x70].compile(begin, "(({p}select)|get)", col, "({p}bases)" "", col, "([^\v]*)", end, nullptr);
  _regex[0x80].compile(begin, "(({p}select)|get)", col, "({p}input)" "", col, "([^\v]*)", end, nullptr);


  //  segment n/m for canu
#ifndef CANU
  //sprintf(_errors, "option '%s' available only with Canu support.", _optString);
#endif

  _regex[0xff].compile(begin, path, end, nullptr);

  _regexLen = 0x100;
}


//  Search all regexes for a match, return the first one.
//  All state is saved in _regex; all we need to know is
//  which one fired.
//
uint32
merylCommandBuilder::matchProgramText(uint64 &pp) {
  uint32 rr = 0;

  if (globals.showDetails() == true) {
    fprintf(stderr, "Search for command word in text starting at position pp=%lu:\n", pp);
    fprintf(stderr, "  %s\n", displayString(_pTxt + pp));
  }

  for (rr=0; rr<_regexLen; rr++)
    if (_regex[rr].match(_pTxt + pp) == true)
      break;

  if (rr == _regexLen)  //  Whoops, no match!
    return _regexLen;

  if (globals.showDetails() == true) {
    for (uint32 ii=0; ii<_regex[rr].numCaptures(); ii++)
      fprintf(stderr, "  0x%02x/%02u: %s %03lu-%03lu '%s'\n", rr, ii,
              _regex[rr].isValid(ii) ? "valid" : "inval",
              _regex[rr].getBgn(ii),
              _regex[rr].getEnd(ii),
              displayString(_regex[rr].get(ii)));
    fprintf(stderr, "\n");
  }

  return rr;
}



void
merylCommandBuilder::processProgramText(void) {
  uint64      v64 = 0;
  double      vD  = 0;

  if (globals.showConstruction() == true)
    fprintf(stderr, "processProgramText()-\n");


  for (uint64 pp=0; pp<_pTxtLen; pp++) {
    uint32            rr = matchProgramText(pp);   //  Magic inside!

    //  For (potential) input matches, try to add the word as an input file.
    //  If it succeeds, we're done (we just need to do the rest of the loop
    //  to advance pp to the end of the word).  If it fais, reset rr
    //  to indicate a no-match word and fall through to reporting an error.
    //
    if ((rr == 0xff) && (isInput(_regex[rr].get(1)) == false))
      rr = _regexLen;

    //  If rr isn't valid, report the non-matching word.
    //
    if (rr >= _regexLen) {
      char *str = new char [_pTxtLen];

      for (uint32 oo=0; (pp < _pTxtLen) && (_pTxt[pp] != '\v'); pp++) {
        char const *sym = displayLetter(_pTxt[pp]);

        while (*sym)
          str[oo++] = *sym++;

        str[oo] = 0;
      }

      sprintf(_errors, "ERROR: word '%s' not recognized.", displayString(str));
      delete [] str;
    }

    //  Set pointers to the option strings of the matching regex.
    //    [0] - the whole match, including any word separators (not useful)
    //    [1] - the whole match, without word separators
    //    [2] - 1st option
    //    [3] - 2nd ...
    //
    //ar const *o0 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 0)) ? nullptr : _regex[rr][0];
    char const *o1 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 1)) ? nullptr : _regex[rr][1];
    char const *o2 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 2)) ? nullptr : _regex[rr][2];
    char const *o3 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 3)) ? nullptr : _regex[rr][3];
    char const *o4 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 4)) ? nullptr : _regex[rr][4];
    char const *o5 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 5)) ? nullptr : _regex[rr][5];
    char const *o6 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 6)) ? nullptr : _regex[rr][6];
    char const *o7 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 7)) ? nullptr : _regex[rr][7];
    char const *o8 = ((rr >= _regexLen) || (_regex[rr].numCaptures() <= 8)) ? nullptr : _regex[rr][8];

    //  If a valid operator, and one that doesn't mess with the operator
    //  stack (not '[' or ']') make sure there is a operator on the stack.

    if ((rr > 0x02) && (rr < _regexLen) && (_opStack.size() == 0))
      addNewOperation();

    //  Grab the current operation (it could be null) and process the token.

    merylOpTemplate  *op = getCurrent();

    //  Word-separator match, do nothing.
    if      (rr == 0x00) {
    }

    //  '['
    //
    //  If the existing operation is empty, do nothing,
    //  otherwise, close the existing operation and add a new one.
    else if (rr == 0x01) {
      if ((op == nullptr) || (op->isEmpty() == false)) {
        terminateOperations(1);
        addNewOperation();
      }
    }

    //  ']'
    //
    //  Close the existing operation, complaining if there isn't an
    //  operation to close.
    else if (rr == 0x02) {
      if (terminateOperations(1, true) == false)
        sprintf(_errors, "processWord()- Extra ']' encountered in command line.");
    }

    //  'count', 'count-forward', 'count-reverse', with options:
    //    'count-size   = <integer>'
    //    'count-suffix = [ACGTacgt]+'
    //
    else if (rr == 0x05) {
      if      (op->_isCounting)
        sprintf(_errors, "ERROR: operation is already a counting operation.\n");

      else if (op->_isSelector)
        sprintf(_errors, "ERROR: operation is a select operation, cannot also be a counting operation.\n");

      else {
        op->_isCounting = true;
        op->_counting   = new merylOpCounting(o2[0]);
      }
    }

    else if (rr == 0x06) {
      if (op->_isCounting == false) {
        sprintf(_errors, "processWord()- option '%s' requires a counting operation.", displayString(o1));
      } else {
        if (strcmp(o2, "expect")  == 0)  op->_counting->setExpectedNumberOfKmers(decodeInteger(o3, 0, 0, v64, _errors));
        if (strcmp(o2, "suffix")  == 0)  op->_counting->setCountSuffix(o3);
        if (strcmp(o2, "segment") == 0) {
#ifndef CANU
          sprintf(_errors, "option '%s' available only with Canu support.", displayString(o1));
#endif
          _segment    = decodeInteger(o3, 0, 0, v64, _errors);
          _segmentMax = decodeInteger(o4, 0, 0, v64, _errors);
        }
      }
    }


    //  OUTPUTS
    //
    else if (rr == 0x1f)  sprintf(_errors, "processWord()- '%s' is ambiguous; use output:show or output:stats.", displayString(o1));

    else if (rr == 0x20)  op->addOutputToDB  (o2,        _errors);  //  output:data=path
    else if (rr == 0x21)  op->addOutputToList(o2, false, _errors);  //  output:list=path
    else if (rr == 0x22)  op->addOutputToList(o2, false, _errors);  //  output:show
    else if (rr == 0x23)  op->addOutputToPipe(o2,        _errors);  //  output:pipe=path
    else if (rr == 0x24)  op->addHistoOutput (o2,        _errors);  //  output:histo(=path)
    else if (rr == 0x25)  op->addStatsOutput (o2,        _errors);  //  output:stats(=path)
    else if (rr == 0x26)  op->addOutputToDB  (o2,        _errors);  //  legacy output

    //  INPUTS
    //
    else if ((rr == 0x28) ||
             (rr == 0x29) ||
             (rr == 0x2a)) {
      merylInput *in = new merylInput;
      bool        sx = false;

      if (rr == 0x27)   sx = in->registerMerylDB   (o2, _errors);
      if (rr == 0x28)   sx = in->registerMerylList (o2, _errors);
      if (rr == 0x29)   sx = in->registerActionPipe(o2, _errors);
      // (rr == 0x2a)   sx = in->registerTemplate  (nullptr,   _errors);

      if (sx)
        op->addInput(in);
      else
        delete in;
    }

    //  Union and Intersect are very similar.
    //   - value assign is count (union) or first (intersect)
    //   - label assign is or    (union) or and   (intersect)
    //   - 
    else if ((rr == 0x09) ||
             (rr == 0x0a)) {            //  Union                      //  Intersect.
      op->_valueAssign = (rr == 0x09) ? merylAssignValue::valueCount : merylAssignValue::valueFirst;
      op->_labelAssign = (rr == 0x09) ? merylAssignLabel::labelOr    : merylAssignLabel::labelAnd;

      merylSelector  f(merylSelectorQuantity::isIndex,
                       merylSelectorRelation::isNOP, false, o1);

      f._input_num_any    = true;

      op->addSelectorToProduct(f);

      if       (o2[0] ==  0) {
        ;  //  Do nothing.
      }
      else if ((o2[0] == 'm') && (o2[1] == 'i') && (o2[2] == 'n')) {
        op->_valueAssign    = merylAssignValue::valueMin;
        op->_valueConstant  = kmvalumax;
        op->_labelAssign    = merylAssignLabel::labelSelected;
      }
      else if ((o2[0] == 'm') && (o2[1] == 'a') && (o2[2] == 'x')) {
        op->_valueAssign    = merylAssignValue::valueMax;
        op->_valueConstant  = kmvalumin;
        op->_labelAssign    = merylAssignLabel::labelSelected;
      }
      else if ((o2[0] == 's') && (o2[1] == 'u') && (o2[2] == 'm')) {
        op->_valueAssign    = merylAssignValue::valueAdd;
        //->_labelAssign    = { same as initially set: 'or' for union and 'and' for intersect }
      }
      else {
        sprintf(_errors, "unknown action '%s' encountered.", o1);   //  Can't actually happen unles regex is busted.
      }
    }

    //  SUBTRACT and DIFFERENCE
    //   - For subtract, the kmer must occur in the first input (input_idx == 1) but can
    //     occur in any number of other inputs (input_num_any).
    //   - For difference, the kmer must occur in the first input, and can only occur
    //     in that input.
    //
    else if ((rr == 0x0b) ||
             (rr == 0x0c)) {            //  Subtract                        //  Difference
      op->_valueAssign = (rr == 0x0b) ? merylAssignValue::valueSub        : merylAssignValue::valueFirst;
      op->_labelAssign = (rr == 0x0b) ? merylAssignLabel::labelDifference : merylAssignLabel::labelFirst;

      merylSelector  f(merylSelectorQuantity::isIndex,
                       merylSelectorRelation::isNOP, false, o1);

      if (rr == 0x0b)   f._input_num_any    = true;  //  Subtract can take any number of inputs
      if (rr == 0x0c)   f._input_num.push_back(1);   //  Difference 

      f._input_idx.push_back(1);

      op->addSelectorToProduct(f);
    }

    //
    else if ((rr == 0x10) ||   //  less-than
             (rr == 0x11) ||   //  greater-than
             (rr == 0x12) ||   //  at-least
             (rr == 0x13) ||   //  at-most
             (rr == 0x14) ||   //  equal-to
             (rr == 0x15)) {   //  not-equal-to
      op->_inputsMin      = 1;
      op->_inputsMax      = 1;

      op->_valueAssign    = merylAssignValue::valueFirst;
      op->_labelAssign    = merylAssignLabel::labelFirst;

      merylSelectorRelation   rel = merylSelectorRelation::isNOP;

      if (rr == 0x10)   rel = merylSelectorRelation::isLt;
      if (rr == 0x11)   rel = merylSelectorRelation::isGt;
      if (rr == 0x12)   rel = merylSelectorRelation::isGeq;
      if (rr == 0x13)   rel = merylSelectorRelation::isLeq;
      if (rr == 0x14)   rel = merylSelectorRelation::isEq;
      if (rr == 0x15)   rel = merylSelectorRelation::isNeq;

      merylSelector  f(merylSelectorQuantity::isValue, rel, false, o1);

      f._vIndex1 = 0;   //  Use value from output kmer.
      f._vValue2 = 0;   //  Reset below, or to distinct or word-freq.

      if      (o2[0] == 'd')   f._vValue2Distinct = strtodouble(o3);
      else if (o2[0] == 'w')   f._vValue2WordFreq = strtodouble(o3);
      else if (o2[0] == 't')   f._vValue2         = decodeInteger(o3, 0, 0, v64, _errors);
      else                     f._vValue2         = decodeInteger(o2, 0, 0, v64, _errors);

      //fprintf(stderr, "values %u %f %f\n", f._vValue2, f._vValue2Distinct, f._vValue2WordFreq);

      op->addSelectorToProduct(f);
    }

    //
    else if ((rr == 0x18) ||   //  increase
             (rr == 0x19) ||   //  decrease
             (rr == 0x1a) ||   //  multiply
             (rr == 0x1b) ||   //  divide
             (rr == 0x1c) ||   //  divide-round
             (rr == 0x1d)) {   //  modulo
      op->_inputsMin      = 1;
      op->_inputsMax      = 1;

      op->_valueAssign    = merylAssignValue::valueNOP;   //  Reset below.
      op->_valueConstant  = 0;                            //  This too.
      op->_labelAssign    = merylAssignLabel::labelFirst;

      if (rr == 0x18)    op->_valueAssign = merylAssignValue::valueAdd;
      if (rr == 0x19)    op->_valueAssign = merylAssignValue::valueSub;
      if (rr == 0x1a)    op->_valueAssign = merylAssignValue::valueMul;
      if (rr == 0x1b)    op->_valueAssign = merylAssignValue::valueDiv;
      if (rr == 0x1c)    op->_valueAssign = merylAssignValue::valueDivZ;
      if (rr == 0x1d)    op->_valueAssign = merylAssignValue::valueMod;

      //  Mul, Div and DivZ allow floats!
      op->_valueConstant = decodeInteger(o2, 0, 0, v64, _errors);
    }


    //  If a valid match, advance past the end of the consumed portion,
    //  then past any word separators, then backup one position
    //  so the loop can increment by one.
    //
    if ((rr < _regexLen) && (_regex[rr].isAccepted())) {
      pp += _regex[rr].getEnd(1);

      while ((pp < _pTxtLen) && (_pTxt[pp] == '\v'))
        pp++;

      pp--;
    }

    if (globals.showConstruction() == true)
      fprintf(stderr, "----------\n");
  }
}
