
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
//  Print a full command tree with indent 'indent'.
//
void
merylCommandBuilder::printTree(uint32 tID, uint32 indent) {
  printTree(getTree(tID), tID+1, 0, "");
}



//
//  Print a single tree node, recursing into child nodes.
//    op     - the tree node to print
//    tID    - integer index representing which tree is being printed
//    iID    - integer index representing which input this node is to some other node
//    indent - amount to indent the display
//
//  tID and iID are mutually exclusive.
//    if tID > 0 - this node is the root of a tree, and its output goes nowhere.
//    if iID > 0 - this node is the input to some other node.
//

void
merylCommandBuilder::printTree(merylOpTemplate *op, uint32 tID, uint32 iID, char const *istr) {
  char   sA[1024];

  //
  //  Report where our kmers go and make an indentation string.
  //
  //    istr - the previous indentation string
  //    jstr - the indentation string for this level
  //         - formed from "{istr-last-character}|  "
  //
  char *jstr = new char [strlen(istr) + 8];
  char *kstr = new char [strlen(istr) + 8];

  if (tID > 0) {
    fprintf(stderr, "%s+-TREE #%u\n",  istr, tID);
    jstr[0] = ' ';
    jstr[1] = ' ';
    jstr[2] = ' ';
    jstr[3] =  0;
  }
  else {
    uint32 jj=0;
    for (; istr[jj]; jj++)
      jstr[jj] = kstr[jj] = istr[jj];

    jstr[jj] = '|';  kstr[jj] = ' ';  jj++;
    jstr[jj] = ' ';  kstr[jj] = ' ';  jj++;
    jstr[jj] = ' ';  kstr[jj] = ' ';  jj++;
    jstr[jj] =  0;   kstr[jj] =  0;
  }

  //
  //  Report all OUTPUTS.
  //

  {
    bool hasOutput = ((op->_outDbseName  != nullptr) ||
                      (op->_outListName  != nullptr) ||
                      (op->_outShowName  != nullptr) ||
                      (op->_outPipeName  != nullptr) ||
                      (op->_outStatsName != nullptr) ||
                      (op->_outHistoName != nullptr));

    if (hasOutput         == true)      fprintf(stderr, "%s|\n",                 jstr);
    if (op->_outDbseName  != nullptr)   fprintf(stderr, "%s|-OUTPUT meryl database '%s'\n",    jstr, op->_outDbseName);
    if (op->_outListName  != nullptr)   fprintf(stderr, "%s|-OUTPUT text list '%s'%s\n",       jstr, op->_outListName, (op->_outList == nullptr) ? "" : " (parallel mode)");
    if (op->_outShowName  != nullptr)   fprintf(stderr, "%s|-OUTPUT text list (STDOUT)\n",     jstr);
    if (op->_outPipeName  != nullptr)   fprintf(stderr, "%s|-OUTPUT meryl pipe '%s'\n",        jstr, op->_outPipeName);
    if (op->_outStatsName != nullptr)   fprintf(stderr, "%s|-OUTPUT statistics report '%s'\n", jstr, op->_outStatsName);
    if (op->_outHistoName != nullptr)   fprintf(stderr, "%s|-OUTPUT histogram report '%s'\n",  jstr, op->_outHistoName);
  }

  //
  //  Report kmer/value/label assignment.
  //

  {
    char ts[1024];

    fprintf(stderr, "%s|\n",                 jstr);
    fprintf(stderr, "%s|-SET value to %s\n", jstr, op->displayValueAssignment(ts));
    fprintf(stderr, "%s|-SET label to %s\n", jstr, op->displayLabelAssignment(ts));
  }

  //
  //  Report selectors
  //

  if (op->_select.size() > 0) {
    uint32  ss = op->_select.size();

    fprintf(stderr, "%s|\n",                 jstr);
    fprintf(stderr, "%s|-SELECT kmer if:\n", jstr);

    for (uint32 ii=0; ii<op->_select.size(); ii++) {
      uint32  ses = op->_select[ii].size();
      auto   &SES = op->_select[ii];
      char    lin = (ii < op->_select.size() - 1) ? '|' : ' ';

      if (ii == 0)
        ;

      for (uint32 jj=0; jj<ses; jj++) {
        if        ((ii == 0) && (jj == 0)) {
          //sparse
          fprintf(stderr, "%s|  |\n", jstr);
          fprintf(stderr, "%s|  +-----+---- %s", jstr, SES[0].describe(sA));
        } else if ((ii > 0) && (jj == 0)) {
          fprintf(stderr, "%s|  +-or--+---- %s", jstr, SES[0].describe(sA));
        } else if (jj < ses-1) {
          fprintf(stderr, "%s|  %c     | and %s", jstr, lin, SES[jj].describe(sA));
        } else {
          fprintf(stderr, "%s|  %c     \\-and %s", jstr, lin, SES[ses-1].describe(sA));
        }
      }

      //sparse
      fprintf(stderr, "%s|  %c\n", jstr, lin);
    }
  }

  //
  //  Report inputs
  //

  for (uint32 ii=0; ii<op->_inputs.size(); ii++) {
    merylInput  *in = op->_inputs[ii];
    char         sp = (ii < op->_inputs.size()-1) ? '+' : '\\';

    if (in->isFromTemplate() == true)   fprintf(stderr, "%s%c-INPUT @%u: action #%d\n",                      jstr, sp, ii+1, in->_template->_ident);
    if (in->isFromDatabase() == true)   fprintf(stderr, "%s%c-INPUT @%u: meryl database '%s'\n",             jstr, sp, ii+1, in->inputName());
    if (in->isFromSequence() == true)   fprintf(stderr, "%s%c-INPUT @%u: sequence file '%s'%s\n",            jstr, sp, ii+1, in->inputName(), in->_squish ? " (homopoly compressed)" : "");
    if (in->isFromStore()    == true)   fprintf(stderr, "%s%c-INPUT @%u: Canu sqStore '%s' (reads %u-%u)\n", jstr, sp, ii+1, in->inputName(), in->_sqBgn, in->_sqEnd);

    if (in->isFromTemplate() == true)   printTree(in->_template, 0, ii+1, jstr);
  }
  fprintf(stderr, "%s\n", jstr);

#if 0
  if (op->_presentInFirst == true)
    fprintf(stderr, "%s|--kmer must be present in first input\n", jstr);

  if (op->_presentInAll == true)
    fprintf(stderr, "%s|--kmer has no database presence requirements\n", jstr);

  if (op->_presentInAll == true)
    fprintf(stderr, "%s|--kmer must be present in all inputs\n", jstr);

  for (uint32 pp=0; pp<op->_inputs.size(); pp++)
    if (op->_presentIn[pp] == true)
      fprintf(stderr, "%s|--kmer must be present in %u inputs\n", jstr, pp);
#endif

#if 0
  if (op->_mathConstant > 0)
    fprintf(stderr, "%sconstant=%lu\n", jstr, op->_mathConstant);

  if (op->_threshold != UINT64_MAX)
    fprintf(stderr, "%sthreshold=%lu\n", jstr, op->_threshold);

  if (op->_fracDist != DBL_MAX)
    fprintf(stderr, "%sfraction-distinct=%f\n", jstr, op->_fracDist);

  if (op->_wordFreq != DBL_MAX)
    fprintf(stderr, "%sword-frequenct=%f\n", jstr, op->_wordFreq);
#endif

  //
  //  If a tree, flag the end.
  //

  if (tID > 0) {
    //fprintf(stderr, "%s|----TREE %u ends.\n", indent, "", tID);
    fprintf(stderr, "\n");
  }

  delete [] kstr;
  delete [] jstr;
}
