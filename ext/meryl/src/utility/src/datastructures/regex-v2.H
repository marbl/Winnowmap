/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
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

#ifndef MERYLUTIL_REGEX_H
#define MERYLUTIL_REGEX_H

//
//  A simple regex parser.
//    BRE - Basic Regular Expressions (POSIX)
//    ERE - Extended Regular Expressions (POSIX)
//    SRE - Simple Regular Expressions (deprecated in favor of BRE)
//    PRE - Perl Regular Expressions
//
//  ERE adds ?, + and |, and removes the requirement that () and {} are escaped.
//  
//  ^      -- matches start of line
//  $       -- matches end of line (either end of string or newline)
//  .      -- matches any character (except newline)
//  []     -- matches any of a set of characters
//  [^]    -- matches any except a set of characters
//  \n     -- BRE only, matches a previously matched group
//  *      -- matches previous element zero or more times (append ? to match minimally)
//  {m,n}  -- matches previous element at least m and at most n times (append ? to match minimally)
//
//  ?      -- matches previous element zero or one times (append ? to match minimally)
//  +      -- matches previous element at lest once (append ? to match minimally)
//  |      -- matches either the expression before or the expression after
//
//  (..)   -- grouping and capturing
//  <..>   -- grouping only - (?:..) is the standard notation
//
//  NOTE:  Edge cases in nested captures exist.
//              (ab(c?)) - will return 'ab' and 'c' for captures, instead of 'abc' and 'c'.
//
//  Still need something for prefix match in a group.
//    ({g}pattern)  - disables capturing
//    ({c}pattern)  - explicitly enables capturing
//    ({p}pattern)  - captures any prefix of the pattern
//
//    ({prefix}example) --> (e(x(a(m(p(le?)?)?)?)?)?)
//
//  Character classes, used within bracket expressions:
//    [:alnum:]  - [A-Za-z0-9]
//    [:alpha:]  - [A-Za-z]
//    [:blank:]  - [ \t]            - space and tab
//    [:cntrl:]  - [\x00-\x1F\x7F]
//    [:digit:]  - [0-9]
//    [:graph:]  - [\x21-\x7E]      - visible characters
//    [:lower:]  - [a-z]
//    [:print:]  - [\x20-\x7E]
//    [:punct:]  - [][!"#$%&'()*+,./:;<=>?@\^_`{|}~-]
//    [:space:]  - [ \t\r\n\v\f]
//    [:upper:]  - [A-Z]
//    [:xdigit:] - [A-Fa-f0-9]
//

#include "types.H"
#include "arrays.H"

namespace merylutil::inline regex::inline v2 {

//  Precedence from high to low
enum class regExTokenType {
  rtGroupBegin,
  rtGroupEnd,
  rtClosure,
  rtConcat,
  rtAlternation,
  rtLineStart,
  rtLineEnd,
  rtCharClass,
  rtEpsilon,
  rtNone,
};



class regExToken {
public:
  char const *display(char *str=nullptr);

  void  unsetMatch(uint8 c)   {             (_syms[c>>6] &= ~(uint64one << (c & 0x3f)));  }
  void    setMatch(uint8 c)   {             (_syms[c>>6] |=  (uint64one << (c & 0x3f)));  }
  bool    isMatch(uint8 c)    { return 0 != (_syms[c>>6] &   (uint64one << (c & 0x3f)));  }

public:
  void  unmatchLetter(uint8 c)                  {                               unsetMatch(c);              }
  void    matchLetter(uint8 c)                  {                                 setMatch(c);              }
  void    matchLetters(char const *s)           { for (uint16 ii=0; s[ii]; ii++)  setMatch(s[ii]);          }
  void    matchRange(uint8 a, uint8 b)          { for (uint16  c=a;  c<=b;  c++)  setMatch(c);              } 

  void    invertMatch(void)                     { for (uint8  ii=0; ii< 4; ii++)  _syms[ii] = ~_syms[ii];   }

  void    mergeMatches(regExToken const &d)     { for (uint8  ii=0; ii< 4; ii++)  _syms[ii] |= d._syms[ii];  }

public:
  void    matchAllSymbols(void)                { _type = regExTokenType::rtCharClass;  matchRange(0x00, 0xff); }

public:
  bool    matchCharacterToken(char const *str, uint64 &nn, bool allowCharacters=false);
public:
  void    matchCharacterClass     (char const *str, uint64 &nn);

  void    makeClosure(char const *str, uint64 ss, uint64 &nn);
  void    makeClosure(uint64 min, uint64 max) { _type = regExTokenType::rtClosure; _min=min; _max=max; };
  void    makeAlternation(void)               { _type = regExTokenType::rtAlternation;  };
  void    makeConcatenation(void)             { _type = regExTokenType::rtConcat;       };

  void    makeGroupBegin(bool pfx, bool cap,                     uint64 &grpIdent, uint64 capActiv, uint64 &capIdent);
  void    makeGroupBegin(char const *str, uint64 ss, uint64 &nn, uint64 &grpIdent, uint64 capActiv, uint64 &capIdent);
  void    makeGroupEnd  (                                        uint64  grpIdent,                  uint64  capIdent);

  void    matchLineStart(void)                { _type = regExTokenType::rtLineStart; };
  void    matchLineEnd(void)                  { _type = regExTokenType::rtLineEnd;   };

public:
  uint64             _id        = 0;
  regExTokenType     _type      = regExTokenType::rtNone;

  uint64             _syms[4]   = { 0, 0, 0, 0 };

  bool               _pfx       = false;
  bool               _cap       = false;
  uint64             _capIdent  = 0;        //  ID of the capture group this token belongs to.
  uint64             _grpIdent  = 0;        //  If a group begin/end, the id of that group - independent from capIdent!

  //std::set<uint64> _capGroups = std::set<uint64>();
  //  replace this with an index to a list of sets

  uint32             _capGrpIdx = uint32max;
  linearset<uint32> *_capGrp    = nullptr;

  uint64             _min       = 0;
  uint64             _max       = uint64max;
};



class regExState {
public:
  void         addMatch(uint64 tid, regExToken tok, bool verbose=false);
  void         addEpsilon(uint64 tid, bool verbose=false);

  uint64      _id         = 0;
  bool        _accepting  = false;

  regExToken  _tok;

  uint64      _mid    =   uint64max;               //  State for tok match.
  uint64      _lid[2] = { uint64max, uint64max };  //  States for lambda transitions.
};


class regExExpr {
public:
  regExState *bgn(regExState *rs) { return rs + _bid; }
  regExState *end(regExState *rs) { return rs + _eid; }

  uint64  _bid = uint64max;
  uint64  _eid = uint64max;
};


class regExMatch {
public:
  regExState *state(regExState *rs) { return rs + _sid; }

  uint64  _sid       = uint64max;  //  Pointer to the state node that we matched.
  uint64  _prevMatch = uint64max;  //  Pointer to the previous regExMatch that led to here.
  uint64  _stringIdx = uint64max;  //  Pointer to the character in the text we matched.
};



class regEx {
public:
  //regEx()  {
  //};
  ~regEx() {
    delete [] tl;
    delete [] rs;
    delete [] capstor;
    delete [] bgnP;
    delete [] endP;
    delete [] lenP;
    delete [] caps;
  };

public:
  void        enableVerbose(char which='a') {
    if ((which == 'a') || (which == 'p'))  vParse   = true;
    if ((which == 'a') || (which == 'c'))  vConvert = true;
    if ((which == 'a') || (which == 'b'))  vBuild   = true;
    if ((which == 'a') || (which == 'm'))  vMatch   = true;
  }
  void        disableVerbose(char which='a') {
    if ((which == 'a') || (which == 'p'))  vParse   = false;
    if ((which == 'a') || (which == 'c'))  vConvert = false;
    if ((which == 'a') || (which == 'b'))  vBuild   = false;
    if ((which == 'a') || (which == 'm'))  vMatch   = false;
  }

public:
  bool        compile(char const *p1, char const *p2=nullptr, ...) {
    char const *pattern = p1;
    char const *pconcat = nullptr;

    if (p2) {
      uint64      pl = 1 + strlen(p1) + strlen(p2);
      char       *p = nullptr;
      va_list     ap;

      va_start(ap, p2);
      for (char const *a = va_arg(ap, char const *); a; a = va_arg(ap, char const *))
        pl += strlen(a);
      va_end(ap);

      pattern = pconcat = p = new char [pl];

      for (char const *a=p1; *a; )  *p++ = *a++;    //  Append p1.
      for (char const *a=p2; *a; )  *p++ = *a++;    //  Append p2.

      va_start(ap, p2);                             //  Append additional patterns.
      for (char const *a = va_arg(ap, char const *); a; a = va_arg(ap, char const *))
        while (*a)
          *p++ = *a++;
      va_end(ap);

      *p = 0;
    }

    parse(pattern);
    convert();
    build();

    delete [] pconcat;

    return true;
  }
  bool        match(char const *query);

  uint64      numCaptures(void)   { return capsLen;   }

  bool        isAccepted(void)    { return accepting; }

  bool        isValid(uint64 cap) { return accepting && (bgnP[cap] < endP[cap]); }
  uint64      getBgn (uint64 cap) { return bgnP[cap]; }
  uint64      getEnd (uint64 cap) { return endP[cap]; }
  char const *get    (uint64 cap) { return caps[cap]; }

  char const *operator[]     (uint64 cap) { return caps[cap]; }

private:
  regExExpr   concat   (regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExExpr a, regExExpr b);
  regExExpr   alternate(regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExExpr a, regExExpr b);
  regExExpr   closure  (regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExExpr a, uint64 min, uint64 max);
  regExExpr   symbol   (regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExToken tok);
  regExExpr   epsilon  (regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExToken tok);

  bool        isAccepting(uint64 sid);
  void        findNodes(uint64      previous,
                        uint64      sid,
                        uint64      currentPos,
                        regExMatch *&mM, uint64 &mMlen, uint64 &mMmax,
                        uint64     *&sB, uint64 &sBlen, uint64 &sBmax);

  bool        parse(char const *pattern);
  bool        convert(void);
  bool        build(void);

private:
  bool         vParse   = false;
  bool         vConvert = false;
  bool         vBuild   = false;
  bool         vMatch   = false;

  uint64       tlMax = 0;         //  A list of tokens, the output of parse().
  uint64       tlLen = 0;
  regExToken  *tl    = nullptr;

  uint64       rsMax = 0;         //  A set of allocated NFA states.
  uint64       rsLen = 0;
  regExState  *rs    = nullptr;

  regExExpr    re;

  uint64       capstorLen = 0;
  uint64       capstorMax = 0;
  char        *capstor    = nullptr;

  bool         accepting  = false;

  linearsets<uint32>   _cgs;

  uint64       capsLen    = 0;
  uint64       capsMax    = 0;
  uint64      *bgnP       = nullptr;   // Begin position of any matches 
  uint64      *endP       = nullptr;   // End position of matches
  uint64      *lenP       = nullptr;   // Length of text in match, internal only!
  char       **caps       = nullptr;   // Match string
};


}  //  merylutil::regex::v2

#endif  //  MERYLUTIL_REGEX_H
