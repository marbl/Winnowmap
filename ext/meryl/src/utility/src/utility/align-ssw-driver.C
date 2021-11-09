
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

#include "align-ssw-driver.H"
#include "align-ssw.H"

#include "arrays.H"
#include "sequence.H"


//                        A   C   T   G   N
//int8 scoreMatrix[25] = {  1, -1, -1, -1,  0,    // A
//                         -1,  1, -1, -1,  0,    // C
//                         -1, -1,  1, -1,  0,    // T
//                         -1, -1, -1,  1,  0,    // G
//                          0,  0,  0,  0,  0 };  // N



sswLib::sswLib(int32 match,
               int32 mismatch,
               int32 gapopen,
               int32 gapextend) {
  setMatchScores(match, mismatch);
  setGapScores(gapopen, gapextend);
}



sswLib::~sswLib() {
  delete [] _intA;
  delete [] _intB;
  delete [] _cigarCode;
  delete [] _cigarValu;
  delete [] _cigarMapBgn;
  delete [] _cigarMapEnd;
  delete [] _aMap;
}



void
sswLib::setMatchScores(int8 match, int8 mismatch) {

  assert(match >= 0);
  assert(mismatch <= 0);

  for (uint32 m=0; m<25; m++)         //  Everything is a mismatch,
    _scoreMatrix[m] = mismatch;

  for (uint32 m=0; m<5; m++)          //  except for matches,
    _scoreMatrix[m + 5 * m] = match;

  for (uint32 m=0; m<5; m++) {        //  and matches to N, which are free.
    _scoreMatrix[ 4 + 5 * m] = 0;
    _scoreMatrix[m + 5 *  4] = 0;
  }

#ifdef SHOW_MATRIX
  for (uint32 p=0, ii=0; ii<5; ii++) {
    for (uint32 jj=0; jj<5; jj++)
      fprintf(stdout, "%3d", _scoreMatrix[p++]);
    fprintf(stdout, "\n");
  }
#endif
}



void
sswLib::setGapScores(int8 open, int8 extend) {

  assert(open   <= 0);
  assert(extend <= 0);

  _gapOpen   = open;
  _gapExtend = extend;
}



bool
sswLib::align(char const *seqA_, uint32 seqlenA_, int32 bgnA_, int32 endA_,
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
      (seqlenB_ < bgnB_) || (endB_ <= bgnB_))
    return(false);

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

  //  Create a 'profile' then align something to it.
  //    The 'profile' is A, the 'read'
  //    The 'result'  is B, the 'ref'
  //
  //  ssw_align is expecting to get gap penalties as the absolute value, but
  //  sswLib (this code) is saving them as a score penalty (i.e., negative).
  //
  //  A gap of length 1 is cost gapOpen.
  //  A gap of length 2 is cost gapOpen + gapExtend.

  s_profile  *profile = ssw_init(_intA, _lenA,
                                 _scoreMatrix,   //  Score matrix as an array
                                 5,              //  Dimension of matrix
                                 1);             //  Score estimate, 1 == high scores expected

  s_align    *result  = ssw_align(profile,
                                  _intB, _lenB, 
                                  -_gapOpen,     //  gap open (absolute value)
                                  -_gapExtend,   //  gap extend (absolute value)
                                  1,             //  flags
                                  0,             //  filter score    if flag1 and not flag0
                                  0,             //  filter distance if flag2 and not flag0
                                  1000);         //  mask len

  uint32  nMismatch = mark_mismatch(result->ref_begin1,
                                    result->read_begin1,
                                    result->read_end1,
                                    _intB,
                                    _intA,
                                    _lenA,
                                    &result->cigar,
                                    &result->cigarLen);

  //  Parse the results.
  //
  //  Coordinates seem to be 0-based inclusive.  We convert to space-based.

  _bgnA = _offA + result->read_begin1;
  _endA = _offA + result->read_end1 + 1;

  _bgnB = _offB + result->ref_begin1;
  _endB = _offB + result->ref_end1 + 1;

  _score = result->score1;

  //  Make space for the alignment, and copy it over.

  resizeArrayPair(_cigarCode, _cigarValu, 0, _cigarMax, result->cigarLen + 1);

  for (int32 cc=0; cc<result->cigarLen; ++cc) {
    _cigarCode[cc] = "MIDNSHP=X"[result->cigar[cc] & 0xf];
    _cigarValu[cc] = result->cigar[cc] >> 4;
  }
  _cigarLen = result->cigarLen;

  _cigarCode[_cigarLen] = 0;
  _cigarValu[_cigarLen] = 0;

  //  Clean up.

  init_destroy(profile);
  align_destroy(result);

  //  Analyze the results.  We cleared these variables at the start.

  analyzeAlignment();

  //  Maybe emit some logging.

  if (verbose_) {
    fprintf(stdout, "\n");
    fprintf(stdout, "A: %6d-%6d  mat %6u  mis %6u  gap %6u  len %6u\n", _bgnA, _endA, _aMat, _aMis, _aGap, _aLen);
    fprintf(stdout, "B: %6d-%6d  score %6d  erate %f%%\n", _bgnB, _endB, _score, 100.0 * _erate);
  }

  //for (uint32 ii=0; ii<_cigarLen; ii++)
  //  fprintf(stdout, "%3u - %4u %c\n", ii, _cigarValu[ii], _cigarCode[ii]);

  return(true);
}



//  Draw the alignment, truncating long match regions to at most 2 *
//  maxMatchLength letters, with dots to indicate bases dropped.
//
void
sswLib::display(uint32 maxMatchLength) {
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
sswLib::analyzeAlignment(void) {

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
