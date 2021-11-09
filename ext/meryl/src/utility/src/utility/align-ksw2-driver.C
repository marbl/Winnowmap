
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

#include "align-ksw2-driver.H"
#include "align-ksw2.H"

#include "arrays.H"
#include "sequence.H"


ksw2Lib::ksw2Lib(int32 matchScore,
                 int32 mismatchScore,
                 int32 gapopenPenalty,
                 int32 gapextendPenalty) {
  setMatchScores(matchScore, mismatchScore);
  setGapPenalties(gapopenPenalty, gapextendPenalty);
}



ksw2Lib::~ksw2Lib() {
  delete [] _intA;
  delete [] _intB;
  delete [] _cigarCode;
  delete [] _cigarValu;
  delete [] _cigarMapBgn;
  delete [] _cigarMapEnd;
  delete [] _aMap;
}



void
ksw2Lib::setMatchScores(int8 match, int8 mismatch) {

  assert(match >= 0);
  assert(mismatch <= 0);

  fprintf(stderr, "MATCH %d MISMATCH %d\n", match, mismatch);

  for (uint32 m=0; m<25; m++)         //  Everything is a mismatch,
    _scoreMatrix[m] = mismatch;

  for (uint32 m=0; m<5; m++)          //  except for matches,
    _scoreMatrix[m + 5 * m] = match;

  for (uint32 m=0; m<5; m++) {        //  and matches to N, which are free.
    _scoreMatrix[ 4 + 5 * m] = 0;
    _scoreMatrix[m + 5 *  4] = 0;
  }

#define SHOW_MATRIX
#ifdef SHOW_MATRIX
  for (uint32 p=0, ii=0; ii<5; ii++) {
    for (uint32 jj=0; jj<5; jj++)
      fprintf(stdout, "%3d", _scoreMatrix[p++]);
    fprintf(stdout, "\n");
  }
#endif
}



void
ksw2Lib::setGapPenalties(int8 open, int8 extend) {

  assert(open   >= 0);
  assert(extend >= 0);

  fprintf(stderr, "OPEN %d EXTEND %d\n", open, extend);

  _gapOpen   = open;
  _gapExtend = extend;
}





bool
ksw2Lib::align(char const *seqA_, uint32 seqlenA_, int32 bgnA_, int32 endA_,
               char const *seqB_, uint32 seqlenB_, int32 bgnB_, int32 endB_, bool verbose_) {

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

  //  Allocate space for at least lenA (lenB) things.

  resizeArray(_intA, 0, _maxA, _lenA);
  resizeArray(_intB, 0, _maxB, _lenB);

  //  Convert the input sequences into integers.

  for (uint32 ii=0; ii<_lenA; ii++)
    _intA[ii] = encode2bitBase(_seqA[_offA + ii]);

  for (uint32 ii=0; ii<_lenB; ii++)
    _intB[ii] = encode2bitBase(_seqB[_offB + ii]);


  ksw_extz_t   ez;

  memset(&ez, 0, sizeof(ksw_extz_t));

  ez.max_q = -1;
  ez.max_t = -1;
  ez.mqe_t = -1;
  ez.mte_q = -1;

  ez.max = 0;
  ez.mqe = KSW_NEG_INF;
  ez.mte = KSW_NEG_INF;

  ez.n_cigar = 0;


  //_flags |= KSW_EZ_SCORE_ONLY;   //  Don't record align path or cigar
  _flags |= KSW_EZ_GENERIC_SC;     //  Without this flag match/mismatch only; last symbol is a wildcard
  _flags |= KSW_EZ_APPROX_MAX;     // approximate max; this is faster with sse
  _flags |= KSW_EZ_APPROX_DROP;    // approximate Z-drop; faster with sse
  //_flags |= KSW_EZ_EXTZ_ONLY;      // only perform extension (no cigar generated)

  //  ksw_extz
  //  ksw_extz2_sse

  ksw_extz2_sse(nullptr,               //  kalloc memory pool
                _lenA, _intA,          //  query
                _lenB, _intB,          //  target
                5,                     //  alphabet size, sqrt of scoreMatrix size
                _scoreMatrix,          //  scores!
                _gapOpen, _gapExtend,  //  penalties!  gap of length ; costs -(gapo + l * gape)
                _bandWidth,            //  bandwidth; negative disables
                _zDrop,                //  off-diagonal drop-off to stop extension; negative disables
                _endBonus,             //  ?? (only in the sse version)
                _flags,                //  KSW_EZ_ flags
                &ez);                  //  output scores and cigar

  //  Parse the results.
  //
  //  Coordinates seem to be 0-based inclusive.  We convert to space-based.

  _bgnA = _offA + 0;
  _endA = _offA + ez.max_q;

  _bgnB = _offB + 0;
  _endB = _offB + ez.max_t;

  _score = ez.score;

  //  Make space for the alignment, and copy it over.

  resizeArrayPair(_cigarCode, _cigarValu, 0, _cigarMax, ez.n_cigar + 1);

  for (int32 cc=0; cc<ez.n_cigar; ++cc) {
    _cigarCode[cc] = "MIDNSHP=X"[ez.cigar[cc] & 0xf];
    _cigarValu[cc] = ez.cigar[cc] >> 4;
  }
  _cigarLen = ez.n_cigar;

  _cigarCode[_cigarLen] = 0;
  _cigarValu[_cigarLen] = 0;

  //  Clean up.

  //  Analyze the results.  We cleared these variables at the start.

  analyzeAlignment();

  //  Maybe emit some logging.

  if (verbose_) {
    fprintf(stdout, "\n");
    fprintf(stdout, "A: %6d-%6d  mat %6u  mis %6u  gap %6u  len %6u\n", _bgnA, _endA, _aMat, _aMis, _aGap, _aLen);
    fprintf(stdout, "B: %6d-%6d  score %6d  erate %f%%\n", _bgnB, _endB, _score, 100.0 * _erate);
  }

  for (uint32 ii=0; ii<_cigarLen; ii++)
    fprintf(stdout, "%3u - %4u %c\n", ii, _cigarValu[ii], _cigarCode[ii]);

  return(true);
}



void
ksw2Lib::analyzeAlignment(void) {

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
