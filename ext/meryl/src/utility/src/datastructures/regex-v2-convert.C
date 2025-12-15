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

//  Convert the tokenized and concat-inserted regExToken list from parse()
//  and converts it to reverse-polish form.
//
//  Iterates over each symbol:
//
//    1) If an alternation, concat or closure operation, pop off any higher
//       precedence operations from the stack and apply them.
//
//    2) If the start of a group, push on the group begin operator.
//
//    3) If the end of a group, pop the stack until we get to the first group
//       open, then discard the open operation.
//
//    4) Otherwise, just move the token - rtNone, rtLineStart, rtLineEnd,
//       rtCharClass - to the output string.

namespace merylutil::inline regex::inline v2 {

bool
regEx::convert(void) {
  uint64       olLen = 0;

  if (vConvert) {
    fprintf(stderr, "\n");
    fprintf(stderr, "CONVERTING\n");
    fprintf(stderr, "\n");
  }

  merylutil::stack<regExToken>  st;   //  Operation stack.

  for (uint64 tt=0; tt<tlLen; tt++) {
    if (vConvert) {
      fprintf(stderr, "tl[%03lu] -- process %s\n", tt, tl[tt].display());
    }

    if ((tl[tt]._type == regExTokenType::rtAlternation) ||
        (tl[tt]._type == regExTokenType::rtConcat) ||
        (tl[tt]._type == regExTokenType::rtClosure)) {
      while ((st.depth() > 0) && (st.top()._type != regExTokenType::rtGroupBegin) && (st.top()._type <= tl[tt]._type)) {
        tl[olLen++] = st.pop();
        if (vConvert)
          fprintf(stderr, "tl[%03lu] <- pop-op  %s\n", olLen-1, tl[olLen-1].display());
      }
      if (vConvert)
        fprintf(stderr, "        -- push\n");
      st.push(tl[tt]);
    }

    else if (tl[tt]._type == regExTokenType::rtGroupBegin) {
      if (vConvert)
        fprintf(stderr, "        -- push\n");
      st.push(tl[tt]);
    }

    else if (tl[tt]._type == regExTokenType::rtGroupEnd) {
      while ((st.depth() > 0) && (st.top()._type != regExTokenType::rtGroupBegin)) {
        tl[olLen++] = st.pop();
        if (vConvert)
          fprintf(stderr, "tl[%03lu] <- pop-end %s\n", olLen-1, tl[olLen-1].display());
      }

      //  Should always have an rtGroupBegin on the stack, if not,
      //  parens are unbalanced.
      if ((st.depth() == 0) || (st.top()._type != regExTokenType::rtGroupBegin))
        fprintf(stderr, "Unbalanced parentheses.\n");

      assert(st.depth() > 0);
      assert(st.top()._type == regExTokenType::rtGroupBegin);

      st.pop();
    }

    else {
      if (olLen != tt)
        tl[olLen] = tl[tt];

      olLen++;
      if (vConvert)
        fprintf(stderr, "tl[%03lu] <- copy    %s\n", olLen-1, tl[olLen-1].display());
    }

    if (vConvert) {
      for (uint64 ii=st.depth(); ii-- > 0; )
        fprintf(stderr, "        -- st[%03lu] %s\n", ii, st[ii].display());
      fprintf(stderr, "\n");
    }
  }

  //  Clear the stack.
  while (st.depth() > 0)
    tl[olLen++] = st.pop();

  tlLen = olLen;

  return true;
}

}  //  merylutil::regex::v2
