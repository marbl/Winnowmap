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

#include "regex-v2.H"
#include "arrays.H"
#include "strings.H"

//  Convert the input string to a list of tokens.
//  Each token represents either a letter to match (or a character class)
//  or an operation.
//    character classes ("[abcde]") are smashed to one token, as are
//    ranges ("{4,19}").  Ranges, '*', '?' and '+' are converted to
//    type Closure.

namespace merylutil::inline regex::inline v2 {




//  Recognizes character classes:
//     '\d'        -> decimal digits
//     '[:digit:]' -> decimal digits
//
//  Recognized non-printable letters and uselessly escaped letters (both only when allowCharacters == true):
//     '\xHH'      -> the 8-bit character represented by hex HH
//     '\OOO'      -> the 8-bit character represented by octal OOO
//
//  Handles escaped letters:
//     '\t', '\n', '\v', '\f', '\r'    -> the usual whitespace letters
//     '\N' (where N is anything else) -> the letter N
//
bool
regExToken::matchCharacterToken(char const *str, uint64 &nn, bool allowCharacters) {
  uint64  zz = 0;     //  A dummy.
  uint64  n1 = nn+1;  //  A handy alias.

  _type = regExTokenType::rtCharClass;

  //fprintf(stderr, "matchCharacterToken()- nn=%lu '%s'\n", nn, str+nn);

  if       (str[nn] == '\\') {
    if      (str[n1] == 'a')  matchCharacterToken("[:alpha:]", zz);
    else if (str[n1] == 'w')  matchCharacterToken("[:alpha:]", zz), matchLetter('_');
    else if (str[n1] == 's')  matchCharacterToken("[:space:]", zz);
    else if (str[n1] == 'd')  matchCharacterToken("[:digit:]", zz);
    else if (str[n1] == 'l')  matchCharacterToken("[:lower:]", zz);
    else if (str[n1] == 'u')  matchCharacterToken("[:upper:]", zz);
    else if (str[n1] == 'p')  matchCharacterToken("[:print:]", zz);

    else if (str[n1] == 'W')  { regExToken d;  d.matchCharacterToken("[:alnum:]", zz); d.matchLetter('_'); d.invertMatch(); mergeMatches(d); }
    else if (str[n1] == 'D')  { regExToken d;  d.matchCharacterToken("[:digit:]", zz);                     d.invertMatch(); mergeMatches(d); }
    else if (str[n1] == 'S')  { regExToken d;  d.matchCharacterToken("[:space:]", zz);                     d.invertMatch(); mergeMatches(d); }

    else if (str[n1] == 't')  matchLetter('\t');  //  0x09 - horizontal tab
    else if (str[n1] == 'n')  matchLetter('\n');  //  0x0a - line feed
    else if (str[n1] == 'v')  matchLetter('\v');  //  0x0d - vertical tab
    else if (str[n1] == 'f')  matchLetter('\f');  //  0x0c - form feed
    else if (str[n1] == 'r')  matchLetter('\r');  //  0x0d - carriage return

    else if (allowCharacters == false)   //  Characters not allowed, so
      return false;                      //  no match found.

    else if  ((str[n1] == 'x')                      && isHexDigit(str[nn+2]) && isHexDigit(str[nn+3])) {
      matchLetter(asciiHexToInteger(str[nn+2]) << 4 |
                  asciiHexToInteger(str[nn+3]));
      nn += 2;  //  We match 2 more than the default 2 letters.
    }

    else if (((str[n1] == '0') || (str[n1] == '1')) && isOctDigit(str[nn+2]) && isOctDigit(str[nn+3])) {
      matchLetter(asciiOctToInteger(str[nn+1]) << 6 |
                  asciiOctToInteger(str[nn+2]) << 3 |
                  asciiOctToInteger(str[nn+3]));
      nn += 2;  //  We match 2 more than the default 2 letters.
    }

    else {
      matchLetter(str[n1]);
    }

    nn += 2;
  }

  else if ((str[nn] == '[') &&
           (str[n1] == ':')) {
    if      (strncmp(str+nn, "[:alnum:]",   9) == 0)  { matchCharacterClass("[A-Za-z0-9]", zz);              nn +=  9; }
    else if (strncmp(str+nn, "[:alpha:]",   9) == 0)  { matchCharacterClass("[A-Za-z]", zz);                 nn +=  9; }
    else if (strncmp(str+nn, "[:blank:]",   9) == 0)  { matchLetters("\t ");                                 nn +=  9; }
    else if (strncmp(str+nn, "[:cntrl:]",   9) == 0)  { matchCharacterClass("[\\x00-\\x1F\\x7F]", zz);       nn +=  9; }
    else if (strncmp(str+nn, "[:digit:]",   9) == 0)  { matchCharacterClass("[0-9]", zz);                    nn +=  9; }
    else if (strncmp(str+nn, "[:graph:]",   9) == 0)  { matchCharacterClass("[\\x21-\\x7E]", zz);            nn +=  9; }
    else if (strncmp(str+nn, "[:lower:]",   9) == 0)  { matchCharacterClass("[a-z]", zz);                    nn +=  9; }
    else if (strncmp(str+nn, "[:print:]",   9) == 0)  { matchCharacterClass("[\\x20-\\x7E]", zz);            nn +=  9; }
    else if (strncmp(str+nn, "[:punct:]",   9) == 0)  { matchLetters("][!\"#$%&'()*+,./:;<=>?@\\^_`{|}~-");  nn +=  9; }
    else if (strncmp(str+nn, "[:space:]",   9) == 0)  { matchLetters("\t\n\v\f\r ");                         nn +=  9; }
    else if (strncmp(str+nn, "[:upper:]",   9) == 0)  { matchCharacterClass("[A-Z]", zz);                    nn +=  9; }
    else if (strncmp(str+nn, "[:xdigit:]", 10) == 0)  { matchCharacterClass("[0-9A-Za-z]", zz);              nn += 10; }
    else
      return false;
  }

  else if (allowCharacters == false) {  //  Not starting with either '\' or '[:', and not
    return false;                       //  allowing letters, so not something we decode.
  }

  else {                                //  Otherwise, we do allow matches to single letters
    matchLetter(str[nn++]);             //  so make the match and advance.
  }

  return true;
}


          
//  Extracts a full character class from within brackets:  '[ (stuff) ]'
//  Iterates over every letter in the brackets, checking:
//    Is it a symbolic token - perl style '\d' or ISO style '[:alpha:]'
//    If not, decode it as a plain letter ('a'), an escaped letter ('\]') or a hex byte ('\xDD').
//    If after the letter there is a dash ('-') decode a second letter and set a range.
//
//  A single '^' at the start of the class will invert the sense.
//
//  But see https://perldoc.perl.org/perlrecharclass
//
void
regExToken::matchCharacterClass(char const *str, uint64 &nn) {
  uint64 iv = 0;   //  If non-zero, invert the class at the end.
  uint8  c1 = 0;   //  The begin character in the range.
  uint8  c2 = 0;   //  The ending character in the range (possible the same as c1).

  //  Decode the letter ('\xHH' or '\OOO' or 'C') at the current position and advance over it.
  auto decode = [](char const *str, uint64 &nn) -> uint8 {
    bool  ishex = (str[nn+0] == '\\') &&  (str[nn+1] == 'x')                        && isHexDigit(str[nn+2]) && isHexDigit(str[nn+3]);
    bool  isoct = (str[nn+0] == '\\') && ((str[nn+1] == '0') || (str[nn+1] == '1')) && isOctDigit(str[nn+2]) && isOctDigit(str[nn+3]);

    if      (ishex) { nn += 4;  return                                     asciiHexToInteger(str[nn-2]) << 4 | asciiHexToInteger(str[nn-1]); }
    else if (isoct) { nn += 4;  return asciiOctToInteger(str[nn-3]) << 6 | asciiOctToInteger(str[nn-2]) << 3 | asciiOctToInteger(str[nn-1]); }
    else            { nn += 1;  return str[nn-1];                                                                                            }
  };

  _type = regExTokenType::rtCharClass;

  assert(str[nn] == '[');
  //fprintf(stderr, "matchCharacterClass()- nn=%lu '%s' (on ENTRY)\n", nn, str+nn);

  if (str[++nn] == '^')      //  Skip over the opening '[' then test for the invert symbol;
    iv = ++nn;               //  if found, set the flag and move past the invert symbol.

  if (str[nn] == ']')        //  If the first letter in the class is a ']', add it to the
    matchLetter(str[nn++]);  //  the class and skip over it; do not just make an empty class.

  while ((str[nn] != 0) && (str[nn] != ']')) {
    //fprintf(stderr, "matchCharacterClass()- nn=%lu '%s' (in LOOP)\n", nn, str+nn);

    //  Handle character classes and escaped letters.
    if (matchCharacterToken(str, nn) == true)
      continue;

    //  Handle inverted letters.
    if (str[nn] == '~') {
      unmatchLetter(decode(str, ++nn));
      continue;
    }

    //  Handle single letters, then ranges.
    c1 = decode(str, nn);                         //  Decode the first letter.

    if (str[nn] != '-') {                         //  If this isn't the start of a
      matchLetter(c1);                            //  range, add the letter and
      continue;                                   //  skip the rest.
    }

    if (str[nn+1] == ']') {                       //  [nn] == '-'.  If the next letter
      matchLetter(c1);                            //  is ']' (to close the class), add
      continue;                                   //  the dash and skip the rest.
    } 

    c2 = decode(str, ++nn);                       //  A range!  Skip the dash (++n),
                                                  //  decode the letter, then move over it.
    matchRange(c1, c2);
    //fprintf(stderr, "matchCharacterClass()- range c1=0x%02x c2=0x%02x nn=%lu ch=%c\n", c1, c2, nn, str[nn]);
  }

  assert(str[nn] != 0);   //  We failed if we're at the end of the string.

  if (iv)
    invertMatch();
}





//  Makes a group begin with specified features.
//
void
regExToken::makeGroupBegin(bool pfx, bool cap,
                           uint64 &grpIdent,      //  Global next new group ident, always incremented
                           uint64  capActiv,      //  Currently active capture group
                           uint64 &capIdent) {    //  Global next new capture ident
  _type     = regExTokenType::rtGroupBegin;

  _pfx      =  pfx;
  _cap      =  cap;
  _capIdent = (cap) ? ++capIdent : capActiv;
  _grpIdent =         ++grpIdent;
}


//  Makes a group begin, decoding features "({group}", "({capture}",
//  "({prefix}" (and abbreviations).
//
void
regExToken::makeGroupBegin(char const *str, uint64 ss, uint64 &nn,
                           uint64 &grpIdent,
                           uint64  capActiv,
                           uint64 &capIdent) {
  bool    pfx = false;
  bool    cap = false;
  uint64  err = 0;

  if ((str != nullptr) &&     //  If a string supplied,
      (str[nn+1] == '{')) {   //  Skip over the opening '(' then test for the modifier symbol;
    nn += 2;                  //  if found, skip both '(' and '{' and decode the modifiers.

    while ((str[nn] != 0) && (str[nn] != '}')) {
      if      (strncmp(str+nn, "group", 5) == 0)    { cap = false;  nn += 5;  }
      else if (strncmp(str+nn, "capture", 7) == 0)  { cap = true;   nn += 7;  }
      else if (strncmp(str+nn, "prefix", 6) == 0)   { pfx = true;   nn += 6;  }
      else if (str[nn] == 'g')                      { cap = false;  nn += 1;  }
      else if (str[nn] == 'c')                      { cap = true;   nn += 1;  }
      else if (str[nn] == 'p')                      { pfx = true;   nn += 1;  }
      else if (str[nn] == ',')                      {               nn += 1;  }
      else
        err = ++nn;
    }
  }

  if (err > 0) {
    fprintf(stderr, "ERROR: expecting 'group', 'capture', 'prefix' in '%s'.\n", displayString(str+ss));
    exit(1);
  }

  makeGroupBegin(pfx, cap, grpIdent, capActiv, capIdent);
}


void
regExToken::makeGroupEnd(uint64 grpIdent, uint64 capIdent) {
  _type     = regExTokenType::rtGroupEnd;
  _grpIdent = grpIdent;
  _capIdent = capIdent;
}



//  Decodes a closure range:
//    #1 - {a,b}
//    #2 - {a}
//    #3 - {a,}
//    #4 - {,b}
//
//  It's pretty tricky code.
//    1)  Skip any open brace.  We'll either be on 'a' or a comma now.
//    2)  If on 'a' (#1, #2, #3), decode it into _min.
//        We'll be on the comma or closing brace now.
//    3)  If on the closing brace (#2), set _max to _min.  We're done.
//    4)  If on the comma (#1, #3, #4), skip over it.
//    5)  If on 'b' (#1 and #4), decode it into _max.
//        We'll be on the closing brace now.
//        For case #3, _max will remain at uint64max.
//    6)  Update string pointer and ensure we're at the closing brace.
void
regExToken::makeClosure(char const *str, uint64 ss, uint64 &nn) {
  char const *dec = str + nn;

  _type = regExTokenType::rtClosure;
  _min  = 0;
  _max  = uint64max;

  //fprintf(stderr, "makeClosure()-- '%s'\n", str+nn);

  if (dec[0] == '{')  dec++;
  if (dec[0] != ',')  dec = strtonumber(dec, _min);
  if (dec[0] == '}')  _max = _min;
  if (dec[0] == ',')  dec++;
  if (dec[0] != '}')  dec = strtonumber(dec, _max);
  if (dec[0] != '}')  fprintf(stderr, "Failed to decode closure range '%s'\n", str + nn);

  if (dec[0] != '}')
    fprintf(stderr, "ERROR: expecting closure range in '%s'\n", str + ss);
  assert(dec[0] == '}');

  //fprintf(stderr, "makeClosure()-- '%s' -> %lu %lu\n", str+nn, _min, _max);

  nn = dec - str;
}






struct groupState {
  bool    pfx   = false;
  bool    cap   = false;
  uint64  depth = 0;
  uint64  cid   = 0;
  uint64  gid   = 0;
};



bool
regEx::parse(char const *str) {
  uint64                        tokIdent = uint64max;
  uint64                        grpIdent = uint64max;
  uint64                        capActiv = 0;
  uint64                        capIdent = uint64max;
  regExToken                    toka     = { ._id=++tokIdent };

  merylutil::stack<groupState>  groupStates;     //  Stack of group information
  merylutil::stack<groupState>  prefxStates;     //  Stack of group information

  uint32                        groupIndex = 0;  //  Index into current active capGroups


  if (vParse) {
    fprintf(stderr, "\n");
    fprintf(stderr, "PARSING\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s\n", displayString(str));
    fprintf(stderr, "\n");
  }

  resizeArray(tl, 0, tlMax, 1024);   //  Allocate an initial 1024 nodes.
  tlLen = 0;                         //  Recycle any existing token list.

  //  Make an initial capture group to get the entirety of the match.

  toka.makeGroupBegin(false, true, grpIdent, capActiv, capIdent);

  tl[tlLen++] = toka;

  groupStates.push(groupState{ .pfx=toka._pfx, .cap=toka._cap, .depth=0, .cid=toka._capIdent, .gid=toka._grpIdent });

  {
    linearset<uint32> grp(groupStates.depth());

    for (uint32 ii=0; ii<groupStates.depth(); ii++)
      grp.insert(groupStates[ii].cid);

    groupIndex = _cgs.insert(grp);
  }


                            
  for (uint64 ss=0, nn=0; str[ss]; ss = nn) {    //  ss is the (constant) start of this token
    toka = { ._id=++tokIdent };                  //  nn is the (computed) start of the next token

    //fprintf(stderr, "parse()-- ss %lu ch %c\n", ss, str[ss]);

    if      (str[ss] == '\\')   toka.matchCharacterToken(str, nn, true), nn--;
    else if (str[ss] == '.')    toka.matchAllSymbols();
    else if (str[ss] == '[')    toka.matchCharacterClass(str, nn);

    else if (str[ss] == '(')  toka.makeGroupBegin(str, ss, nn, grpIdent, capActiv, capIdent);
    else if (str[ss] == ')')  toka.makeGroupEnd(groupStates.top().gid, groupStates.top().cid);
#warning mismatched parens will crash the above

    else if (str[ss] == '*')  toka.makeClosure(0, uint64max);
    else if (str[ss] == '?')  toka.makeClosure(0, 1);
    else if (str[ss] == '+')  toka.makeClosure(1, uint64max);
    else if (str[ss] == '{')  toka.makeClosure(str, ss, nn);

    else if (str[ss] == '|')  toka.makeAlternation();

    else                      toka.matchCharacterToken(str, nn, true), nn--;

    nn++;  //  Move to the next input character.


    //  Hooray!  Parsing of the token is finished.

    //  Ensure space exists for any new tokens we add:
    //    Possibly one extra for a prefix match.
    //    Always one for the symbol we just parsed.
    //    Up to 2 * prefixDepth if we're at the end of a prefix group.
    //    Possible one concatenation operator,
    //    
    while (tlLen + 1 + 1 + 2 * groupStates.top().depth + 1 >= tlMax)
      resizeArray(tl, tlLen, tlMax, tlMax + 1024);

    //  If the last symbol is a group-begin or an alternation and this is an
    //  alternation, add an epsilon match (use the id from the token to add,
    //  and bump that id).  We always have a previous symbol, from the
    //  initial capture group.
    //
    if ((toka._type == regExTokenType::rtAlternation) &&
        ((tl[tlLen-1]._type == regExTokenType::rtGroupBegin) ||
         (tl[tlLen-1]._type == regExTokenType::rtAlternation))) {
      tl[tlLen++] = { ._id=toka._id++, ._type=regExTokenType::rtEpsilon, ._capIdent=toka._capIdent, ._grpIdent=toka._grpIdent };
    }

    //  If we're a new closure, the last symbol must be a group end or match.
    if ((toka._type == regExTokenType::rtClosure) &&
        ((tl[tlLen-1]._type != regExTokenType::rtGroupEnd) &&
         (tl[tlLen-1]._type != regExTokenType::rtCharClass))) {
      fprintf(stderr, "Closure must follow only group end or matches.\n");
      assert(0);
    }

    //  If we're in a prefix group, append a non-capture group begin before
    //  every symbol except the first, and count how many were appended:
    //
    //    ({p}1234) -> (1(2(3(4?)?0?)?)?)
    //                   ^ ^ ^ - appended group begin

    if ((groupStates.top().pfx == true) &&                //  The last one is a prefix group
        ((toka._type == regExTokenType::rtCharClass) ||   //  The current token is a matching operator,
         (toka._type == regExTokenType::rtGroupBegin)) && //   or a group begin
        (groupStates.top().depth++ > 0)) {                //  And not the first   (Also count how many)
      tl[tlLen++] = { ._id=++tokIdent, ._type=regExTokenType::rtGroupBegin, ._grpIdent=++grpIdent };
      prefxStates.push(groupState{ .gid=grpIdent });
    }

    if ((groupStates.top().pfx == true) &&                //  Blow up on alternation in prefix groups,
        (toka._type == regExTokenType::rtAlternation)) {  //  because it makes no sense.
      fprintf(stderr, "Alteration not supported in prefix groups.\n");
      assert(0);
    }

    //  If this is a new group, remember group state.

    if (toka._type == regExTokenType::rtGroupBegin) {
      groupStates.push(groupState{ .pfx=toka._pfx, .cap=toka._cap, .depth=0, .cid=toka._capIdent, .gid=toka._grpIdent });

      linearset<uint32> grp(groupStates.depth());

      for (uint32 ii=0; ii<groupStates.depth(); ii++)
        grp.insert(groupStates[ii].cid);

      groupIndex = _cgs.insert(grp);
    }

    //  If this is a match operator, assign the match
    //  to the current capture group.

    if (toka._type == regExTokenType::rtCharClass) {
      toka._capGrpIdx =      groupIndex;
      toka._capGrp    = _cgs[groupIndex];

      //for (uint32 ii=0; ii<groupStates.depth(); ii++)
      //  toka._capGroups.insert(groupStates[ii].cid);
    }

    //  If we've just seen a group-end symbol, pop off the group-info for this group.
    //  If this is closing a prefix-group, terminate all the interal groups we made.

    if (toka._type == regExTokenType::rtGroupEnd) {
      groupState gs = groupStates.pop();

      if (gs.pfx == true) {
        while (gs.depth-- > 1) {
          groupState ps = prefxStates.pop();

          tl[tlLen++] = { ._id=++tokIdent, ._type=regExTokenType::rtGroupEnd, ._grpIdent=ps.gid };
          tl[tlLen++] = { ._id=++tokIdent, ._type=regExTokenType::rtClosure,  ._min=0, ._max=1 };
        }
      }

      {
        linearset<uint32> grp(groupStates.depth());

        for (uint32 ii=0; ii<groupStates.depth(); ii++)
          grp.insert(groupStates[ii].cid);

        groupIndex = _cgs.insert(grp);
      }

      //  Find the next active capture group
      //capIdent = groupStates.top().cid;
      //capIdent = gs.cid;
    }

    //  Append whatever token we made.

    tl[tlLen++] = toka;

    if (vParse) {
      fprintf(stderr, "tl[%03lu] -- %s\n", tlLen-1, toka.display());
      for (uint32 xx=groupStates.depth(); xx--; )
        fprintf(stderr, "        -- groupState[%02d] pfx:%d cap:%d depth:%lu cid:%lu gid:%lu\n", xx,
                groupStates[xx].pfx,
                groupStates[xx].cap,
                groupStates[xx].depth,
                groupStates[xx].cid,
                groupStates[xx].gid);
      fprintf(stderr, "\n");
    }


    //  Append a concatentaion operator if there is a next symbol and:
    //   - this was NOT a start group or alternation operator
    //        SYM ( concat     - makes no sense
    //        SYM | concat     - makes no sense
    //
    //   - the next symbol is NOT a closure operator (*, {, ?, +)
    //        SYM concat *     - makes no sense
    //
    //   - the next symbol is NOT an alternate operator or a group end
    //        SYM concat |     - makes no sense
    //        SYM concat )     - makes no sense
    //

    if ((str[nn] != 0) &&
        (str[ss] != '(') && (str[ss] != '|') &&
        (str[nn] != '*') && (str[nn] != '{') && (str[nn] != '?') && (str[nn] != '+') &&
        (str[nn] != '|') && (str[nn] != ')'))
      tl[tlLen++] = { ._id=++tokIdent, ._type=regExTokenType::rtConcat };

    assert(tlLen <= tlMax);
  }

  //  Close the overall capture group.
  tl[tlLen++] = { ._id=++tokIdent, ._type=regExTokenType::rtGroupEnd, ._capIdent=0, ._grpIdent=0 };

  if (vParse) {
    for (uint64 ii=0; ii<tlLen; ii++)
      fprintf(stderr, "tl[%03lu] -- %s\n", ii, tl[ii].display());
    fprintf(stderr, " -- %lu capture groups\n", capsLen+1);
  }

  assert(groupStates.top().gid == 0);
  assert(groupStates.top().cid == 0);
  assert(groupStates.depth()   == 1);

  //  Allocate space for resuts.  By construction, there is always at least
  //  one capture group - the one enclosing the whole match.

  capsLen = capIdent + 1;

  resizeArray(bgnP, endP, lenP, caps, 0, capsMax, capsLen);


  return true;
}

}  //  merylutil::regex::v2
