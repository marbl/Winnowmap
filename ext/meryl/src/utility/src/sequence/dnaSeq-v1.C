
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

#include "dnaSeq-v1.H"
#include "arrays.H"

namespace merylutil::inline sequence::inline v1 {

dnaSeq::dnaSeq() {
}


dnaSeq::~dnaSeq() {
  releaseIdent();
  releaseBases();
}


void
dnaSeq::releaseIdent(void) {
  delete [] _headr;   _headr = nullptr;  _headrMax = 0;
  delete [] _ident;   _ident = nullptr;  _identMax = 0;
  delete [] _flags;   _flags = nullptr;  _flagsMax = 0;
}


void
dnaSeq::releaseBases(void) {
  delete [] _seq;     _seq   = nullptr;  _seqMax   = 0;
  delete [] _qlt;     _qlt   = nullptr;  _seqLen   = 0;
}


bool
dnaSeq::copy(char  *bout,
             uint32 bgn, uint32 end, bool terminate) {

  if ((end < bgn) || (_seqLen < end))
    return(false);

  for (uint32 ii=bgn; ii<end; ii++)
    bout[ii-bgn] = _seq[ii];

  if (terminate)
    bout[end-bgn] = 0;

  return(true);
}


bool
dnaSeq::copy(char  *bout,
             uint8 *qout,
             uint32 bgn, uint32 end, bool terminate) {

  if ((end < bgn) || (_seqLen < end))
    return(false);

  for (uint32 ii=bgn; ii<end; ii++) {
    bout[ii-bgn] = _seq[ii];
    qout[ii-bgn] = _qlt[ii];
  }

  if (terminate) {
    bout[end-bgn] = 0;
    qout[end-bgn] = 0;
  }

  return(true);
}


void
dnaSeq::findNameAndFlags(void) {
  uint32 hh=0, ii=0, ff=0, td=0;

  if (_identMax < _headrMax)   allocateArray(_ident, _identMax, _headrMax);
  if (_flagsMax < _headrMax)   allocateArray(_flags, _flagsMax, _headrMax);

  while (isWhiteSpace(_headr[hh]) == true)       //  Skip white space before the name.
    hh++;                                        //  Why do you torture us?

  while (isVisible(_headr[hh]) == true)          //  At the start of a name;
    _ident[ii++] = _headr[hh++];                 //  copy it to the ident.
  _ident[ii] = 0;                                //

  if (isTab(_headr[hh]) == true)                 //  If a tab after the ident, assume
    td = hh++;                                   //  flags are tab-delimited,
  else                                           //  else skip whitespace until the
    while (isWhiteSpace(_headr[hh]) == true)     //  first flag.
      hh++;                                      //

  while (isNUL(_headr[hh]) == false)             //  Copy flags until
    _flags[ff++] = _headr[hh++];                 //  the end of string.
  _flags[ff] = 0;

  if (td == 0)                                   //  If no tabs found, trim back
    while ((ff > 0) &&                           //  any whitespace at the end.
           (isWhiteSpace(_flags[ff-1])))
      _flags[--ff] = 0;
}

}  //  namespace merylutil::sequence::v1
