
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

#include "system.H"

#ifndef ALIGN_KSW2_DRIVER_H
#define ALIGN_KSW2_DRIVER_H

class ksw2Lib {
public:
  ksw2Lib(int32 matchScore       =  1,
          int32 mismatchScore    = -2,
          int32 gapopenPenalty   = -2,
          int32 gapextendPenalty = -1);
  ~ksw2Lib();

  void     setMatchScores(int8 match, int8 mismatch);
  void     setGapPenalties(int8 open, int8 extend);

  bool     align(char const *seqA, uint32 lenA,
                 char const *seqB, uint32 lenB, bool verbose=false) {
    return(align(seqA, lenA, 0, lenA,
                 seqB, lenB, 0, lenB, verbose));
  };

  bool     align(char const *seqA, uint32 lenA, int32 bgnA, int32 endA,
                 char const *seqB, uint32 lenB, int32 bgnB, int32 endB, bool verbose=false);

  double   percentIdentity(void)     { return(100.0 * (1.0 - _erate)); };
  double   errorRate(void)           { return(_erate);                 };
  uint32   score(void)               { return(_score);                 };

  uint32   numMatches(void)          { return(_aMat); };
  uint32   numMisMatches(void)       { return(_aMis); };
  uint32   numGaps(void)             { return(_aGap); };
  uint32   alignmentLength(void)     { return(_aLen); };

  uint32   bgnA(void)                { return(_bgnA); };
  uint32   endA(void)                { return(_endA); };

  uint32   bgnB(void)                { return(_bgnB); };
  uint32   endB(void)                { return(_endB); };

  uint32   cigarLength(void)         { return(_cigarLen);     };
  char     cigarCode(uint32 i)       { return(_cigarCode[i]); };
  uint32   cigarValu(uint32 i)       { return(_cigarValu[i]); };

  void     display(uint32 maxMatchLength=20);

private:
  void     analyzeAlignment(void);

private:  //  Inputs
  int8        _scoreMatrix[25] = { 0 };
  int8        _gapOpen   = 0;
  int8        _gapExtend = 0;

  int32       _bandWidth = -1;
  int32       _zDrop     = -1;
  int32       _endBonus  =  0;
  int32       _flags     =  0;

  uint32      _maxA = 0;         //  Space allocated in _intA.
  uint32      _lenA = 0;         //  Length of sequence we're aligning - NOT the length of _seqA.
  const char *_seqA = nullptr;   //  Pointer to input array.
  uint8      *_intA = nullptr;   //  Encoded bases we're aligning.
  uint32      _offA = 0;         //  Offset into _seqA that we're trying to align.

  uint32      _maxB = 0;
  uint32      _lenB = 0;
  const char *_seqB = nullptr;
  uint8      *_intB = nullptr;
  uint32      _offB = 0;

private:  //  Results
  uint32      _bgnA = 0;
  uint32      _endA = 0;

  uint32      _bgnB = 0;
  uint32      _endB = 0;

  int32       _score = 0;

  uint32      _aLen = 0;
  uint32      _aMis = 0;
  uint32      _aGap = 0;
  uint32      _aMat = 0;

  double      _erate = 100.0;

  uint32      _cigarMax  = 0;
  uint32      _cigarLen  = 0;
  char       *_cigarCode = nullptr;
  uint32     *_cigarValu = nullptr;

  uint32      _cigarMapMax = 0;
  uint32     *_cigarMapBgn = nullptr;   //  Index into aMap for each cigar element.
  uint32     *_cigarMapEnd = nullptr;   //  (C-style semantics)

  class sswLibMap {
  public:
    void  init(uint32 a, uint32 b, uint32 c)   { _aPos = a;  _bPos = b;  _cIdx = c; };

    uint32  _aPos;   //  Position in A, space-based
    uint32  _bPos;   //  Position in B
    uint32  _cIdx;   //  Index into cigar array.
  };

  uint32      _aMapMax = 0;
  sswLibMap  *_aMap    = nullptr;       //  Map showing A-B correspondence.
};


#endif  //  ALIGN_KSW2_DRIVER_H
