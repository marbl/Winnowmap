
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

#include "dnaSeqFile-v1.H"

#include "arrays.H"
#include "strings.H"

namespace merylutil::inline sequence::inline v1 {

//
//  File management.
//

bool
htsSeqFile::open(char const *fn, bool indexed) {
  dnaSeqFile::filename(fn);

  _hts        = hts_open(filename(), "r");
  _htshdr     = sam_hdr_read(_hts);
  _indexable  = false;
  _reopenable = false;
  _compressed = false;

  if (_hts == nullptr)
    return false;

  else if (_hts->format.category != htsFormatCategory::sequence_data) {
    fprintf(stderr, "Unable to open '%s': HTS type '%s' contains no sequence data.\n",
            filename(), hts_format_description(&_hts->format));
    return close();
  }

  //fprintf(stderr, "Opened '%s' file '%s'.\n",
  //        hts_format_description(&_hts->format), filename());

  _htsbam = bam_init1();

  return true;
}


bool
htsSeqFile::close(void) {
  if (_hts) {
    bam_destroy1(_htsbam);     _htsbam = nullptr;
    sam_hdr_destroy(_htshdr);  _htshdr = nullptr;
    sam_close(_hts);           _hts    = nullptr;
    unloadIndex();
  }
  return false;
}


bool
htsSeqFile::reopen(void) {
  close();
  if (_reopenable == false)
    return false;
  return open(filename());
}

//
//  Index loading/saving/creation - NOT SUPPORTED.
//

bool
htsSeqFile::loadIndex(bool create) {
  return false;
}

void
htsSeqFile::unloadIndex(void) {
}

void
htsSeqFile::createIndex(void) {
}

void
htsSeqFile::destroyIndex(void) {
}

//
//  Sequence accessors
//

bool
htsSeqFile::findSequence(uint64 i) {
  return false;
}

uint64
htsSeqFile::sequenceLength(uint64 i) {
  assert(0);
  return 0;
}

//
//
//

bool
htsSeqFile::loadSequence(char  *&name, uint32 &nameMax,
                         char  *&seq,
                         uint8 *&qlt,  uint64 &seqMax, uint64 &seqLen, uint32 &error) {

  do                                              //  Load the next sequence, skipping
    if (sam_read1(_hts, _htshdr, _htsbam) < 0)    //  non-primary records and returning
      return false;                               //  end-of-file on EOF AND errors.
  while ((_htsbam->core.flag & BAM_FSECONDARY) ||
         (_htsbam->core.flag & BAM_FSUPPLEMENTARY));

  //  Resize name, seq and qlt arrays for the new sequence.

  resizeArray    (name,      0, nameMax, (uint32)strlen(bam_get_qname(_htsbam)) + 1);
  resizeArrayPair(seq,  qlt, 0, seqMax,  (uint64)_htsbam->core.l_qseq + 1);

  strcpy(name, bam_get_qname(_htsbam));          //  Copy name into output.

  //fprintf(stderr, "%s flags 0x%x %s\n", name, _htsbam->core.flag, bam_flag2str(_htsbam->core.flag));

  uint8_t *bamseq = bam_get_seq(_htsbam);        //  Copy bases into output.
  uint8_t *bamqlt = bam_get_qual(_htsbam);

  if (bam_is_rev(_htsbam) == false) {
    for (seqLen=0, _htspos=0; (_htspos < _htsbam->core.l_qseq); seqLen++, _htspos++) {
      seq[seqLen] = ".AC.G...T......N"[bam_seqi(bamseq,_htspos)];   //  Forward from BGN.
      qlt[seqLen] =                             bamqlt[_htspos];
    }
  }
  else {
    for (seqLen=0, _htspos=_htsbam->core.l_qseq; (_htspos-- > 0); seqLen++) {
      seq[seqLen] = ".TG.C...A......N"[bam_seqi(bamseq,_htspos)];   //  Reverse from END.
      qlt[seqLen] =                             bamqlt[_htspos];
    }
  }

  seq[seqLen] = 0;
  qlt[seqLen] = 0;

  return true;
}


bool
htsSeqFile::loadBases(char    *seq,
                      uint64   maxLength,
                      uint64  &seqLength,
                      bool    &endOfSequence) {

  if (_htspos >= _htsbam->core.l_qseq) {             //  Load a new sequence if needed,
    do                                               //  skipping non-primary records,
      if (sam_read1(_hts, _htshdr, _htsbam) < 0)     //  and returning no-more-bases on EOF
        return false;                                //  AND errors.
    while ((_htsbam->core.flag & BAM_FSECONDARY) ||
           (_htsbam->core.flag & BAM_FSUPPLEMENTARY));

    _htspos = (bam_is_rev(_htsbam) == false) ? 0 : _htsbam->core.l_qseq;
  }

  if (bam_is_rev(_htsbam) == false) {
    for (uint8_t *bamseq = bam_get_seq(_htsbam);  //  Copy bases (starting at wherever
         ((seqLength < maxLength) &&              //  _htspos is currently at) into output.
          (_htspos < _htsbam->core.l_qseq)); seqLength++, _htspos++)   //  NOTE bam_seqi() is a macro
      seq[seqLength] = ".AC.G...T......N"[bam_seqi(bamseq, _htspos)];  //  and evals the 2nd arg twice!
  }
  else {
    for (uint8_t *bamseq = bam_get_seq(_htsbam);  //  For reverse sequences, start at the
         ((seqLength < maxLength) &&              //  end and copy backwards.
          (_htspos-- > 0)); seqLength++)
      seq[seqLength] = ".TG.C...A......N"[bam_seqi(bamseq, _htspos)];
  }

  endOfSequence = (_htspos >= _htsbam->core.l_qseq);   //  htspos == uint64max for reverse.

  return seqLength > 0;
}

}  //  namespace merylutil::sequence::v1
