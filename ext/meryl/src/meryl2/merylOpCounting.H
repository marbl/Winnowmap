
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

#ifndef MERYLOPCOUNTING_H
#define MERYLOPCOUNTING_H

#ifndef MERYLINCLUDE
#error "Do not use merylOpCounting.H, use meryl.H instead."
#endif


//  For the simple counting algorithm, how big of a word to use for initial counts.
typedef uint16 lowBits_t;


class merylOpCounting {
public:
  merylOpCounting(char mode)  {
    switch (mode) {
      case 'f':  case 'F':  _countForward   = true;  break;
      case 'r':  case 'R':  _countReverse   = true;  break;
      default:              _countCanonical = true;  break;
    }
  };
  ~merylOpCounting() {};

public:
  void    doCounting     (std::vector<merylInput *> &inputs,
                          uint64                     allowedMemory,
                          uint32                     allowedThreads,
                          merylFileWriter           *output);
private:
  void    countSimple    (std::vector<merylInput *> &inputs,
                          uint64                     memoryAllowed,
                          uint32                     threadsAllowed,
                          merylFileWriter           *output);

  void    countThreads   (std::vector<merylInput *> &inputs,
                          uint64                     memoryAllowed,
                          uint32                     threadsAllowed,
                          merylFileWriter           *output);

  void    countSequential(std::vector<merylInput *> &inputs,
                          uint64                     memoryAllowed,
                          uint32                     threadsAllowed,
                          merylFileWriter           *output);

public:
  void    setDefaultLabel(kmlabl l)  { _lConstant = l; }
  void    setDefaultValue(kmvalu v)  { _vConstant = v; }

  void    setCountSuffix(char const *s) {
    _countSuffixLength = strlen(s);

    for (uint32 ii=0; ii<_countSuffixLength; ii++) {
      _countSuffixString[ii] = s[ii];

      _countSuffix.addR(s[ii]);
    }
  }

  void    setExpectedNumberOfKmers(uint64 n) {
    _expNumKmers  = n;
  }

private:
  uint64  guesstimateNumberOfkmersInInput_dnaSeqFile(dnaSeqFile *sequence);
  uint64  guesstimateNumberOfkmersInInput_sqStore(sqStore *store, uint32 bgnID, uint32 endID);
  uint64  guesstimateNumberOfkmersInInput(std::vector<merylInput *> &inputs);

  void    findBestValues(uint64 nKmers, uint32 bestPrefix, uint64 memoryUsed);

public:
  bool      _onlyConfig     = false;

  bool      _countCanonical = false;
  bool      _countForward   = false;
  bool      _countReverse   = false;

  char      _countSuffixString[65] = {0};
  uint32    _countSuffixLength     =  0;
  kmer      _countSuffix;

  uint64    _expNumKmers = 0;

  //  do I want to move all the wPrefix etc parameters to here?
  //  labelConstant too?
  //  why not inputs and output then?

  //  Parameters used by countThreads() and countSequential().

  //  countSimple() uses similar names, but they're set at a different place.

  uint32    _wPrefix      = 0;   //  Number of bits in the prefix (== bucket address)
  uint32    _wSuffix      = 0;   //  Number of bits in the suffix
  kmdata    _wSuffixMask  = 0;   //  Mask to get the suffix

  //  wData -> wSuffix

  uint64    _nPrefix      = 0;   //  Number of prefixes there are (== number of buckets)

  //uint64    _sMask      = 0;   //  
  //uint64    _cSuffix    = 0;   //  

  kmvalu    _vConstant    = 0;   //  If non-zero, all kmers have this value.
  kmlabl    _lConstant    = 0;   //  All kmers have this label.
};

#endif  //  MERYLOPCOUNTING_H
