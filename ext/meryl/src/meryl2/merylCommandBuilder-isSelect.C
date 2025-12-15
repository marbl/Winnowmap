
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
//  isRelation() returns the length of any relation-operator found at
//  position p in string s.
//    0 - if none found
//    1 - = > <
//    2 - == eq != <> ne <= le >= ge lt gt
//
//  decodeRelation() returns the relation itself.  If no relation is found,
//  and error is added to the _errors list.
//
//  Note that invalid relations typically generate errors on the surrounding
//  numbers: "451equal123" will decode as "451", "eq", "ual123" and then
//  report that the latter isn't a valid number.
//

uint32
merylCommandBuilder::isRelation(char const *s, uint32 p) {

  if (((s[p] == '=') && (s[p+1] == '=')) ||
      ((s[p] == 'e') && (s[p+1] == 'q')) ||
      ((s[p] == '!') && (s[p+1] == '=')) ||
      ((s[p] == '<') && (s[p+1] == '>')) ||
      ((s[p] == 'n') && (s[p+1] == 'e')) ||
      ((s[p] == '<') && (s[p+1] == '=')) ||
      ((s[p] == 'l') && (s[p+1] == 'e')) ||
      ((s[p] == '>') && (s[p+1] == '=')) ||
      ((s[p] == 'g') && (s[p+1] == 'e')) ||
      ((s[p] == 'l') && (s[p+1] == 't')) ||
      ((s[p] == 'g') && (s[p+1] == 't')))
    return 2;

  if ((s[p] == '=') ||
      (s[p] == '<') ||
      (s[p] == '>'))
    return 1;

  return 0;
}


merylSelectorRelation
merylCommandBuilder::decodeRelation(char const *s, uint32 p) {
  merylSelectorRelation  relation = merylSelectorRelation::isNOP;

  if      (strncmp(s+p, "==", 2) == 0)   { relation = merylSelectorRelation::isEq;  }
  else if (strncmp(s+p, "=",  1) == 0)   { relation = merylSelectorRelation::isEq;  }
  else if (strncmp(s+p, "eq", 2) == 0)   { relation = merylSelectorRelation::isEq;  }

  else if (strncmp(s+p, "!=", 2) == 0)   { relation = merylSelectorRelation::isNeq; }
  else if (strncmp(s+p, "<>", 2) == 0)   { relation = merylSelectorRelation::isNeq; }
  else if (strncmp(s+p, "ne", 2) == 0)   { relation = merylSelectorRelation::isNeq; }

  else if (strncmp(s+p, "<=", 2) == 0)   { relation = merylSelectorRelation::isLeq; }
  else if (strncmp(s+p, "le", 2) == 0)   { relation = merylSelectorRelation::isLeq; }
  else if (strncmp(s+p, ">=", 2) == 0)   { relation = merylSelectorRelation::isGeq; }
  else if (strncmp(s+p, "ge", 2) == 0)   { relation = merylSelectorRelation::isGeq; }

  else if (strncmp(s+p, "<",  1) == 0)   { relation = merylSelectorRelation::isLt;  }
  else if (strncmp(s+p, "lt", 2) == 0)   { relation = merylSelectorRelation::isLt;  }
  else if (strncmp(s+p, ">",  1) == 0)   { relation = merylSelectorRelation::isGt;  }
  else if (strncmp(s+p, "gt", 2) == 0)   { relation = merylSelectorRelation::isGt;  }

  else {
    sprintf(_errors, "No comparison operator found in '%s',", _optString);
    sprintf(_errors, "  expecting one of '==', 'eq', '!=', 'ge', '<', etc.");
    sprintf(_errors, "");
  }

  if (globals.showConstruction() == true)
    fprintf(stderr, "decodeRelation()- Found relation '%s'\n", toString(relation));

  return relation;
}



void
merylCommandBuilder::decodeSelector(char const *s, uint32 p, merylSelector &f) {

  //  Output values.  Each arg is either a db-index or a constant.

  uint32   index1 = uint32max,  index2 = uint32max;
  uint64   const1 = 0,          const2 = 0;

  //  Find pointers to the various bits of the string.
  //    ab..ae - the first argument word
  //    rb..re - the first relation encountered in the string
  //    bb..be - the second argument word

  uint32 ab = p;                           //  start of the first arg is easy, we're already there.
  uint32 rb = ab;                          //  Search for the first relation.
  while ((s[rb] != 0) &&
         (isRelation(s, rb) == 0))
    rb++;

  uint32 ae = rb;                          //  The end of the first arg is the beginning of the relation,
  //if (s[ae-1] == ':')                      //  unless there is a ':' involved.
  //  ae--;                                  //  I THINK was to catch "value:>3"

  uint32 bb = rb + isRelation(s, rb);      //  second arg starts after the relation.
  //if (s[bb] == ':')                        //   and any ':'.
  //  bb++;

  uint32 be = bb;                          //  The end of the second arg is the end of the string.
  while (s[be] != 0)
    be++;

  if (s[bb] == 0) {                        //  Fail if there is no second argument.
    sprintf(_errors, "Invalid selector '%s': no second argument to comparison operator found.", s);
    sprintf(_errors, "");
  }

  if (globals.showConstruction() == true) {
    fprintf(stderr, "decodeSelector()- WORD          '%s'\n",               s);
    fprintf(stderr, "decodeSelector()- ARG1  %3d-%3d '%s'\n", ab, ae, s + ab);
    fprintf(stderr, "decodeSelector()- RELA  %3d     '%s'\n", rb,        s + rb);
    fprintf(stderr, "decodeSelector()- ARG2  %3d-%3d '%s'\n", bb, be, s + bb);
  }

  //  Decode the first argument.
  if      (ab == rb)      index1 = 0;                                    //  No first arg.
  else if (s[ab] == '@')  decodeInteger(s, ab+1, ae, index1, _errors);   //  First arg is file index.
  else if (s[ab] == '#')  decodeInteger(s, ab+1, ae, const1, _errors);   //  First arg is integer constant.
  else                    decodeInteger(s, ab,   ae, const1, _errors);   //  First arg is integer constant, by default.

  //  Decode the relation.
  f._r = decodeRelation(s, rb);

  //  Decode the decond argument, with some special cases for weird constructs.

  if      (strncmp(s+bb, "distinct=",        9) == 0)  f._vValue2Distinct = strtodouble(s+bb+9);
  else if (strncmp(s+bb, "word-freq=",      10) == 0)  f._vValue2WordFreq = strtodouble(s+bb+10);
  else if (strncmp(s+bb, "word-frequency=", 15) == 0)  f._vValue2WordFreq = strtodouble(s+bb+15);
  else if (strncmp(s+bb, "threshold=",      10) == 0)  decodeInteger(s, bb+10, be, const2, _errors);
  else if (s[bb] == '@')                               decodeInteger(s, bb+ 1, be, index2, _errors);
  else if (s[bb] == '#')                               decodeInteger(s, bb+ 1, be, const2, _errors);
  else                                                 decodeInteger(s, bb,    be, const2, _errors);

  //  Set the index and constants in the selector.
  //
  //  If an index was not specified, we'll set vIndex to uint32max (the
  //  default value for index1 and index2), which is code for 'use the
  //  constant', and the correct constant will be set to whatever was
  //  specified.
  //
  //  If an index is specified, index1 (index2) will not be uint32max, code
  //  for 'use the value from the kmer in database i; the output kmer if i=0`
  //  and the constants will be at their default value of zero.
  //
  //  Finally, if the two indexes are the same, the selector is constant true
  //  of false.  They're either both uint32max, and so the selector is
  //  comparing two constants, or both database indeces refering to the same
  //  kmer.

  f._vIndex1 = index1;
  f._vIndex2 = index2;

  if (f._vIndex1 == f._vIndex2) {
    sprintf(_errors, "Invalid selector '%s': always true (or false).", s);
    sprintf(_errors, "");
  }

  switch (f._q) {
    case merylSelectorQuantity::isValue:
      f._vValue1 = const1;
      f._vValue2 = const2;
      break;
    case merylSelectorQuantity::isLabel:
      f._vLabel1 = const1;
      f._vLabel2 = const2;
      break;
    case merylSelectorQuantity::isBases:
      f._vBases1 = const1;
      f._vBases2 = const2;
      break;
      //case merylSelectorQuantity::isIndex:
      //  f._vIndex1 = const1;
      //  f._vIndex2 = const2;
      //  break;
    default:
      assert(0);
      break;
  }
}





//  Decide if _curParam is a value selector.  If so, decode the
//  stuff and add it to the current selector product term.
//
//  Value selectors can look like (omitting the spaces):
//    value:          OP constant   - both of these use an implicit @1 on the
//    value:          OP @index     - left side; and are 'more natural'
//
//    value: @index   OP constant   - these allow crazy stuff like comparing
//    value: @index   OP @index     - @3<@2 then outputting the value from @1
//    value: constant OP @index     - 
//
//    value: constant OP            - technically will complete the set, but
//    value: @index   OP            - seem awkward to use
//  
bool
merylCommandBuilder::isValueSelector(void) {

  if ((_curClass != opClass::clSelect) &&
      (_curPname != opPname::pnValue))
    return false;

  merylSelector  f(merylSelectorQuantity::isValue,
                   merylSelectorRelation::isNOP,
                   _invertNextSelector,
                   _curParam);

  decodeSelector(_curParam, 0, f);

  getCurrent()->addSelectorToProduct(f);

  return true;
}


//  Decide if _curParam is a label selector.  If so, decode the
//  stuff and add it to the current selector product term.
//
//  Label selectors look like value selectors.
//
bool
merylCommandBuilder::isLabelSelector(void) {

  if ((_curClass != opClass::clSelect) &&
      (_curPname != opPname::pnLabel))
    return false;

  merylSelector  f(merylSelectorQuantity::isLabel,
                   merylSelectorRelation::isNOP,
                   _invertNextSelector,
                   _curParam);

  decodeSelector(_curParam, 0, f);

  getCurrent()->addSelectorToProduct(f);

  return true;
}



//  Decide if _curParam is a base content selector.  If so, decode the
//  stuff and add it to the current selector product term.
//
//  Base content selectors are slightly different than value and label selectors
//  as they also specify what bases to count, and it makes no sense to
//  compare kmers in different databases (they are all the same).
//    bases:acgt: OP constant
//
bool
merylCommandBuilder::isBasesSelector(void) {

  if ((_curClass != opClass::clSelect) &&
      (_curPname != opPname::pnBases))
    return false;

  merylSelector  f(merylSelectorQuantity::isBases,
                   merylSelectorRelation::isNOP,
                   _invertNextSelector,
                   _curParam);

  //  Decode the bases string itself.

  uint32  bpos = 0;

  while ((_curParam[bpos] != ':') && (_curParam[bpos] != 0)) {
    switch (_curParam[bpos]) {
      case 'a':
      case 'A':
        f._countA = true;
        break;
      case 'c':
      case 'C':
        f._countC = true;
        break;
      case 't':
      case 'T':
        f._countT = true;
        break;
      case 'g':
      case 'G':
        f._countG = true;
        break;
      default:
        sprintf(_errors, "Invalid 'bases' letter in selector '%s'.", _curParam);
        sprintf(_errors, "");
        break;
    }

    bpos++;
  }

  if (_curParam[bpos] != ':') {
    sprintf(_errors, "Failed to parse 'bases' selector '%s'.", _curParam);
    sprintf(_errors, "");
    return true;
  }

  //  Pass the rest of the string to the usual selector decoding to get the
  //  operation and constant.

  decodeSelector(_curParam, bpos, f);

  //  Make sure the user didn't specify a useless index.
  //
  //  vIndex must be either 0 (use the output kmer) or uint32max (use the constant).
  //  vIndex is 0 if nothing is supplied for this side: "bases:acgt:ge4" will set vIndex1 to 0.

  if      ((f._vIndex1 != uint32max) &&
           (f._vIndex1  > 0)) {
    sprintf(_errors, "selector '%s' right hand side cannot specify a database input.", _curParam);
    sprintf(_errors, "");
  }

  if      ((f._vIndex2 != uint32max) &&
           (f._vIndex2  > 0)) {
    sprintf(_errors, "selector '%s' left hand side cannot specify a database input.", _curParam);
    sprintf(_errors, "");
  }
  
  getCurrent()->addSelectorToProduct(f);

  return true;
}



bool
merylCommandBuilder::isInputSelector(void) {

  if ((_curClass != opClass::clSelect) &&
      (_curPname != opPname::pnInput))
    return false;

  merylSelector  f(merylSelectorQuantity::isIndex,
                   merylSelectorRelation::isNOP,
                   _invertNextSelector,
                   _curParam);

  //  The 'input' selector is a ':' or ',' separated list specifying how
  //  many and which input databases a kmer must be present in.
  //
  //  How many input databases the kmer must be present in:
  //    'n'      = in exactly n files
  //    'n-m'    = in between n and m inclusive
  //    'all'    - in all
  //    'any'    - in any number (== '1-all', the default)
  //    'n-all'  = in at least n
  //
  //  Which input files the kmer must be present in:
  //    'first'  - in the first input file (== '@1')
  //    '@n'     = in the nth input file
  //    '@n-@m'  = in input files n-m inclusive
  //

  //for (uint32 ii=0; ii<_curParamLen; ii++)
  //  if (_curParam[ii] == ',')
  //    _curParam[ii] = ':';

  splitToWords  W(_curParam, ":,");

  for (uint32 ww=1; ww<W.numWords(); ww++) {     //  Skip the first word; it is 'input'.
    splitToWords  P(W[ww], '-');

    //  Kmer must be present in all input databases.
    if      (strcmp(W[ww], "all") == 0) {
      //fprintf(stderr, "PUSH input_num_all\n");
      f._input_num_all = true;
    }

    //  Kmer must be present in any number of input databases.
    //  This is the default if nothing else is specified.
    else if (strcmp(W[ww], "any") == 0) {
      //fprintf(stderr, "PUSH input_num_any\n");
      f._input_num_any = true;
    }

    //  Kmer must be present in the first database.
    //  Equivalent to @1, and implemented as such.
    else if (strcmp(W[ww], "first") == 0) {
      //fprintf(stderr, "PUSH input_idx 1\n");
      f._input_idx.push_back(1);
    }

    //  @a: Kmer must be present in a specific input file.
    else if ((P.numWords() == 1) &&
             (P[0][0] == '@') && (isDecInteger(P[0]+1) == true)) {
      uint32 a = strtouint32(P[0]+1);

      //fprintf(stderr, "PUSH input_idx a   <- %u\n", a);
      f._input_idx.push_back(a);
    }

    //  @a-@b:  Kmer must be present in input files a through b, inclusive.
    else if ((P.numWords() == 2) &&
             (P[0][0] == '@') && (isDecInteger(P[0]+1) == true) &&
             (P[1][0] == '@') && (isDecInteger(P[1]+1) == true)) {
      uint32 a = strtouint32(P[0]+1);
      uint32 b = strtouint32(P[1]+1);

      for (uint32 x=a; x<=b; x++) {
        //fprintf(stderr, "PUSH input_idx a-b <- %u\n", x);
        f._input_idx.push_back(x);
      }
    }

    //  a:  Kmer must occur in a input files.
    else if ((P.numWords() == 1) &&
             (isDecInteger(P[0]) == true)) {
      uint32 a = strtouint32(P[0]);

      //fprintf(stderr, "PUSH a   input_num <- %u\n", a);
      f._input_num.push_back(a);
    }

    //  a-b:  Kmer must occur in between a and b input files, inclusive.
    else if ((P.numWords() == 2) &&
             (isDecInteger(P[0]) == true) &&
             (isDecInteger(P[1]) == true)) {
      uint32 a = strtouint32(P[0]);
      uint32 b = strtouint32(P[1]);

      for (uint32 x=a; x<=b; x++) {
        //fprintf(stderr, "PUSH a-b input_num <- %u\n", x);
        f._input_num.push_back(x);
      }
    }

    //  a-all:  Kmer must occur in at least a input files.
    else if ((P.numWords() == 2) &&
             (isDecInteger(P[0]) == true) &&
             (strcmp(P[1], "all") == 0)) {
      uint32 a = strtouint32(P[0]);

      //fprintf(stderr, "PUSH input_num_at_least <- u\n", a);
      f._input_num_at_least = std::min(a, f._input_num_at_least);
    }

    else {
      sprintf(_errors, "selector '%s' cannot be decoded: unknown word '%s'.", _curParam, W[ww]);
      sprintf(_errors, "");
    }
  }

  getCurrent()->addSelectorToProduct(f);

  return true;
}


bool
merylCommandBuilder::isSelect(void) {

  if ((isValueSelector() == false) &&
      (isLabelSelector() == false) &&
      (isBasesSelector() == false) &&
      (isInputSelector() == false))
    return false;

  _invertNextSelector = false;           //  Successfully decoded, reset state.

  _curClass = opClass::clNone;
  _curPname = opPname::pnNone;
  _curParam = nullptr;

  return true;
}


bool
merylCommandBuilder::isSelectConnective(void) {
  merylOpTemplate *t = getCurrent();
  bool             a = false;

  if      (strcmp(_optString, "not") == 0)         //  'not' inverts the sense of the next term.
    _invertNextSelector = !_invertNextSelector;    //
  else if (strcmp(_optString, "and") == 0)         //  'and' is syntactic sugar; don't need or use it.
    ;                                              //
  else if (strcmp(_optString, "or") == 0)          //  'or' adds a new product term (later).
    a = true;                                      //
  else                                             //  Whatever it is, it isn't a connective,
    return false;                                  //  and we're all done with selectors.

  if (a) {                                         //  Add a new product term, here, outside the if block above
    if (t->isSelectorProductEmpty() == true) {     //  to keep that block clean and tidy.
      sprintf(_errors, "attempt to add new selector product when existing selector product is empty.");
      sprintf(_errors, "");
    }
    else
      t->addNewSelectorProduct();
  }

  return true;
}
