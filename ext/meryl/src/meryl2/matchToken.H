
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

//
//  Expects inputs of
//    opt='<option>[:=]<value>' - the input word to decode
//    pat='<pattrn>[:=]'        - the pattern we're expecting to find
//
//  Returns true if <option> is a prefix of <pattrn> and the separator
//  matches.  optout is set to the letter after the separator, e.g., <value>.
//      opt='abc:xyz'  pat='abcde:'  ->  true, optout='xyz'
//      opt='abc:'     pat='abcde:'  ->  true, optout=NUL
//
//  If opt does not end with a separator but is still a prefix of <pattrn>,
//  true is returned and optout is set to nullptr.  pat is not required to
//  end with a separator.
//      opt='abc'      pat='abcde:'  ->  true, optout=nullptr
//      opt='abc'      pat='abcde'   ->  true, optout=nullptr
//
//  If pat does not end with a separator, opt must also not end with a
//  separator.
//      opt='abc:xyz'  pat='abcde'   ->  false
//      opt='abc'      pat='abcde'   ->  true, optout=nullptr
//
//  If 'exact' is set, opt must equal pat; if opt has no separator, the
//  separator on pat is ignored.
//
inline
bool
matchToken(char const *opt, char const *&optout, char const *pat, bool exact=false) {
  optout = nullptr;

  while ((*opt !=  0) && (*opt != ':') && (*opt != '=') &&   //  Advance opt and pat to
         (*pat !=  0) &&                                     //  the separator in opt or
         (*opt == *pat))                                     //  first disagreement.
    opt++, pat++;

  if ((*opt != 0) && (*opt != ':') && (*opt != '='))         //  If opt isn't at its end,
    return false;                                            //  there is no match.

  if (exact == false)                                        //  Advance pat to its end
    while ((*pat != 0) && (*pat != ':') && (*pat != '='))    //  unless we're looking
      pat++;                                                 //  for an exact match.

  if ((*opt == 0) && ((*pat == ':') || (*pat == '=')))       //  If opt is at its end,
    pat++;                                                   //  skip over any pat sep.

  if (*opt != *pat)                                          //  A match will be:
    return false;                                            //    0, ':' or '='.

  if (*opt)                                                  //  If a sep, return
    optout = opt + 1;                                        //  the next letter.

  return true;
}
