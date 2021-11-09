
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

#include "align-parasail-driver.H"

#include "arrays.H"
#include "sequence.H"


parasailLib::parasailLib(int32 matchScore,
                         int32 mismatchScore,
                         int32 gapopenPenalty,
                         int32 gapextendPenalty) {
  setMatchScores(matchScore, mismatchScore);
  setGapPenalties(gapopenPenalty, gapextendPenalty);
}



parasailLib::~parasailLib() {

  if (_scoreMatrix)
    parasail_matrix_free(_scoreMatrix);

  delete [] _cigarCode;
  delete [] _cigarValu;
  delete [] _cigarMapBgn;
  delete [] _cigarMapEnd;
  delete [] _aMap;
}



void
parasailLib::setMatchScores(int8 match, int8 mismatch) {

  assert(match >= 0);
  assert(mismatch <= 0);

  _scoreMatrix = parasail_matrix_create("ACGTN", match, mismatch);
}



void
parasailLib::setGapPenalties(int8 open, int8 extend) {

  assert(open   >= 0);
  assert(extend >= 0);

  _gapOpen   = open;
  _gapExtend = extend;
}



bool
parasailLib::align(char const *seqA_, uint32 seqlenA_, int32 bgnA_, int32 endA_,
                   char const *seqB_, uint32 seqlenB_, int32 bgnB_, int32 endB_,
                   bool verbose_,
                   parasail_function_t *pf) {

  //  Clear the results.  Erate is set to max, for ease of discarding alignment failures.

  _bgnA = _endA = 0;
  _bgnB = _endB = 0;

  _score = 0;

  _aLen  = 0;
  _aMis  = 0;
  _aGap  = 0;
  _aMat  = 0;

  _erate = 100.0;

  //  Silently adjust input ranges if they exceed the limits of the sequence.
  //  Return failure if they make no sense.
  //
  //  After this block, forget about the NAME_ parameters.
  //
  //  Importantly, seqlenA_ and seqlenB_ are NEVER used except for checking that the bgn-end range is valid.

  if (bgnA_ < 0)         bgnA_ = 0;
  if (bgnB_ < 0)         bgnB_ = 0;

  if (seqlenA_ < endA_)  endA_ = seqlenA_;
  if (seqlenB_ < endB_)  endB_ = seqlenB_;

  if ((seqlenA_ < bgnA_) || (endA_ <= bgnA_) ||
      (seqlenB_ < bgnB_) || (endB_ <= bgnB_)) {
    fprintf(stdout, "ERROR %u < %u OR %u < %u OR %u < %u OR %u < %u\n",
            seqlenA_, bgnA_, endA_, bgnA_,
            seqlenB_, bgnB_, endB_, bgnB_);
    return(false);
  }

  _offA = bgnA_;   _seqA = seqA_;   _lenA = endA_ - bgnA_;   //  _lenA IS NOT seqlenA_.
  _offB = bgnB_;   _seqB = seqB_;   _lenB = endB_ - bgnB_;   //  It's the length we're aligning.

  //  Align.

  parasail_result_t *result = (*pf)(_seqA + _offA, _lenA,
                                    _seqB + _offB, _lenB,
                                    _gapOpen, _gapExtend, _scoreMatrix);

  //  Parse the results.
  //
  //  Coordinates seem to be 0-based inclusive.  We convert to space-based.

  parasail_cigar_t *cigar = parasail_result_get_cigar(result,
                                                      _seqA + _offA, _lenA,
                                                      _seqB + _offB, _lenB, _scoreMatrix);

  _bgnA = _offA + cigar->beg_query;
  _endA = _offA + parasail_result_get_end_query(result) + 1;

  _bgnB = _offB + cigar->beg_ref;
  _endB = _offB + parasail_result_get_end_ref(result) + 1;

  _score = parasail_result_get_score(result);

  //  Make space for the alignment, and copy it over.

  uint32  ccPos = 0;   //  Canu Cigar position

  uint32  pcPos = 0;   //  Parasail Cigar position
  uint32  pcLen = cigar->len;

  resizeArrayPair(_cigarCode, _cigarValu, 0, _cigarMax, cigar->len + 1);

  //  If the alignment begins with a gap, remove it and adjust the positions.

  if        (parasail_cigar_decode_op(cigar->seq[0]) == 'I') {
    _bgnA += parasail_cigar_decode_len(cigar->seq[0]);

    ccPos = 0;
    pcPos = 1;
  } else if (parasail_cigar_decode_op(cigar->seq[0]) == 'D') {
    _bgnB += parasail_cigar_decode_len(cigar->seq[0]);

    ccPos = 0;
    pcPos = 1;
  }

  //  Decode and copy the alignment.

  while (pcPos < pcLen) {
    _cigarCode[ccPos] = parasail_cigar_decode_op (cigar->seq[pcPos]);
    _cigarValu[ccPos] = parasail_cigar_decode_len(cigar->seq[pcPos]);

    ccPos++;
    pcPos++;
  }

  //  If the alignment ends with a gap, just remove it.  The positions are
  //  already adjusted.

  if      (_cigarCode[ccPos-1] == 'I') {
    ccPos--;
  }

  else if (_cigarCode[ccPos-1] == 'D') {
    ccPos--;
  }

  _cigarLen = ccPos;

  _cigarCode[_cigarLen] = 0;
  _cigarValu[_cigarLen] = 0;

  //  Clean up.

  parasail_cigar_free(cigar);
  parasail_result_free(result);

  //  Analyze the results.  We cleared these variables at the start.

  analyzeAlignment();

  //  Maybe emit some logging.

  if (verbose_) {
    fprintf(stdout, "\n");
    fprintf(stdout, "A: %6d-%6d  score=%d matches=%u  mismatches=%u  gaps=%u  length=%u\n", _bgnA, _endA, _score, _aMat, _aMis, _aGap, _aLen);
    fprintf(stdout, "B: %6d-%6d  %.4f%%\n", _bgnB, _endB, percentIdentity());

    for (uint32 ii=0; ii<_cigarLen; ii++)
      fprintf(stdout, " %u%c", _cigarValu[ii], _cigarCode[ii]);
    fprintf(stdout, "\n");
  }


  return(true);
}



//  Draw the alignment, truncating long match regions to at most 2 *
//  maxMatchLength letters, with dots to indicate bases dropped.
//
void
parasailLib::display(uint32 maxMatchLength) {
  char const  *aseq = _seqA + bgnA();
  char const  *bseq = _seqB + bgnB();

  char        *aaln = new char [16 + _aLen + 16 + 1];
  char        *baln = new char [16 + _aLen + 16 + 1];

  uint32 apos = 0;  //  Position in the output strings.
  uint32 bpos = 0;

  //  Start off showing how much sequence wasn't aligned.

  if (bgnA() > 0)
    sprintf(aaln, "%6u ", bgnA());
  else
    sprintf(aaln, "%6s ", "");

  if (bgnB() > 0)
    sprintf(baln, "%6u ", bgnB());
  else
    sprintf(baln, "%6s ", "");

  apos = strlen(aaln);
  bpos = strlen(baln);

  //  Walk down the cigar string, adding bases as appropriate.

  for (uint32 c=0; c<cigarLength(); c++) {
    aaln[apos++] = ' ';    aaln[apos++] = cigarCode(c);    aaln[apos++] = ' ';
    baln[bpos++] = ' ';    baln[bpos++] = '*';             baln[bpos++] = ' ';

    switch (cigarCode(c)) {
      case 'M':
      case '=':
      case 'X':
        if (cigarValu(c) < maxMatchLength) {
          for (uint32 p=0; p<cigarValu(c); p++) {
            aaln[apos++] = *aseq++;
            baln[bpos++] = *bseq++;
          }
        } else {
          uint32  dl = (maxMatchLength - 6) / 2;

          for (uint32 p=0; p < dl; p++) {
            aaln[apos++] = *aseq++;
            baln[bpos++] = *bseq++;
          }

          sprintf(aaln + apos, "..%u..", cigarValu(c) - dl - dl);
          sprintf(baln + bpos, "..%u..", cigarValu(c) - dl - dl);

          apos = strlen(aaln);
          bpos = strlen(baln);

          aseq += cigarValu(c) - dl - dl;
          bseq += cigarValu(c) - dl - dl;

          for (uint32 p=0; p < dl; p++) {
            aaln[apos++] = *aseq++;
            baln[bpos++] = *bseq++;
          }
        }
        break;

      case 'I':
        for (uint32 p=0; p<cigarValu(c); p++) {
          aaln[apos++] = *aseq++;
          baln[bpos++] = '-';
        }
        break;

      case 'D':
        for (uint32 p=0; p<cigarValu(c); p++) {
          aaln[apos++] = '-';
          baln[bpos++] = *bseq++;
        }
        break;

      default:
        assert(0);
        break;
    }

    aaln[apos] = 0;
    baln[bpos] = 0;

  }

  fprintf(stdout, "\n");
  fprintf(stdout, "A %s\n", aaln);
  fprintf(stdout, "B %s\n", baln);
  fprintf(stdout, "\n");

  delete [] aaln;
  delete [] baln;
}



void
parasailLib::analyzeAlignment(void) {

  //  Count the length of the alignment, and matches/mismatches/etc.

  for (uint32 cc=0; cc<_cigarLen; cc++) {
    char   code = _cigarCode[cc];
    uint32 valu = _cigarValu[cc];

    switch (code) {
      case 'M': _aLen += valu;                  break;  //  Non-gap, either a match or a mismatch.
      case 'I': _aLen += valu;  _aGap += valu;  break;  //  Insertion, gap in query.
      case 'D': _aLen += valu;  _aGap += valu;  break;  //  Deletion,  gap in target.
      case 'N':                                 break;  //  Skipped in reference.
      case 'S':                                 break;  //  Softmask, query doesn't appear in alignment.
      case 'H':                                 break;  //  Hardclip, query doesn't appear in alignment.
      case 'P':                                 break;  //  Padding, silent deletion from reference.
      case '=': _aLen += valu;  _aMat += valu;  break;  //  Match.
      case 'X': _aLen += valu;  _aMis += valu;  break;  //  Mismatch,
      default:
        assert(0);
        break;
    }
  }

  //  Compute the same erate as overlapper does.

  _erate = (double)(_aMis + _aGap) / std::min((_endA - _bgnA), (_endB - _bgnB));

  //  Allocate stuff for building a map between the A and B sequences and the
  //  cigar string.

  resizeArrayPair(_cigarMapBgn, _cigarMapEnd, 0, _cigarMapMax, _cigarLen);
  resizeArray    (_aMap,                      0, _aMapMax,     _aLen);

  uint32       apos = _bgnA;
  uint32       bpos = _bgnB;
  uint32       ai   = 0;

  for (uint32 cc=0; cc<_cigarLen; cc++) {
    char   code = _cigarCode[cc];
    uint32 valu = _cigarValu[cc];

    _cigarMapBgn[cc] = ai;

    switch (code) {
      case 'M':
      case '=':
      case 'X':
        for (uint32 p=0; p<valu; p++)
          _aMap[ai++].init(apos++, bpos++, cc);
        break;

      case 'I':
        for (uint32 p=0; p<valu; p++)
          _aMap[ai++].init(apos++, bpos, cc);
        break;

      case 'D':
        for (uint32 p=0; p<valu; p++)
          _aMap[ai++].init(apos, bpos++, cc);
        break;

      default:
        assert(0);
        break;
    }

    _cigarMapEnd[cc] = ai;

  }

  assert(_aLen == ai);
}
