
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

#ifndef MERYL_UTIL_KMER_ITERATOR_H
#define MERYL_UTIL_KMER_ITERATOR_H

#ifndef MERYL_UTIL_KMER_H
#error "include kmers.H, not this."
#endif


//  Converts a buffer of characters (or a file of characters) into kmers, one
//  kmer at a time.

class kmerIterator {
public:
  kmerIterator(void) {
    assert(kmer::merSize() > 0);
    reset();
    addSequence(NULL, 0);
  };
  kmerIterator(FILE *input);
  kmerIterator(char const *buffer, uint64 bufferLen) {
    assert(kmer::merSize() > 0);
    reset();
    addSequence(buffer, bufferLen);
  };

  void       reset(void) {
    _kmerSize   = _fmer.merSize();
    _kmerLoad   = 0;
    _kmerValid  = _fmer.merSize() - 1;
  };

  void       addSequence(char const *buffer, uint64 bufferLen) {
    _buffer    = buffer;
    _bufferLen = bufferLen;
    _bufferPos = 0;
  };

  //
  //  The primary interface.  Iterate over all valid mers, silently skipping over
  //  any invalid ones.
  //

  bool       nextMer(void) {
  nextMer_anotherBase:
    if (_bufferPos >= _bufferLen)      //  No more sequence, and not a valid kmer.
      return(false);

    if (isACGT(_bufferPos) == false) {
      _kmerLoad = 0;                   //  Not a valid base.  Clear the current
      _bufferPos++;                    //  kmer and move to the next base.
      goto nextMer_anotherBase;
    }

    _fmer.addR(_buffer[_bufferPos]);   //  A valid base, so push it onto
    _rmer.addL(_buffer[_bufferPos]);   //  the kmer.

    _bufferPos++;

    if (_kmerLoad < _kmerValid) {      //  Not a full kmer, keep loading
      _kmerLoad++;                     //  bases from the buffer.
      goto nextMer_anotherBase;
    }

    return(true);                      //  Valid kmer!
  };


  //
  //  Alternate interface.  Iterate over all bases.  Use isValid() to test if the kmer
  //  ending at this base is valid.  Use isACGTbgn() and isACGTend() to decide if the
  //  base at the start/end of the kmer is ACGT.
  //

  bool       isValid(void) {
    return(_kmerLoad == _kmerSize);
  }

  bool       isACGT(uint64 pos) {
    return((_buffer[pos] == 'A') || (_buffer[pos] == 'a') ||
           (_buffer[pos] == 'C') || (_buffer[pos] == 'c') ||
           (_buffer[pos] == 'G') || (_buffer[pos] == 'g') ||
           (_buffer[pos] == 'T') || (_buffer[pos] == 't'));
  }

  bool       isACGTbgn(void) {
    return(isACGT(_bufferPos - _kmerSize));
  }

  bool       isACGTend(void) {
    return(isACGT(_bufferPos - 1));
  }

  char       getBaseBgn(void) {return(_buffer[_bufferPos - _kmerSize]);}
  char       getBaseEnd(void) {return(_buffer[_bufferPos - 1]);}

  bool       nextBase(void) {

    //  Preload the first kmerSize-1 bases.

    while ((_bufferPos < _kmerSize - 1) &&
           (_bufferPos < _bufferLen)) {
      if (isACGT(_bufferPos) == false) {   //  Not a valid base, reset the counter.
        _kmerLoad = 0;
      } else {
        _fmer.addR(_buffer[_bufferPos]);   //  A valid base, so push it onto
        _rmer.addL(_buffer[_bufferPos]);   //  the kmer.

        _kmerLoad++;
      }

      _bufferPos++;
    }

    //  Stop if we're out of sequence.

    if (_bufferPos >= _bufferLen)        //  No more sequence, stop.
      return(false);

    //  Load another base.

    if (isACGT(_bufferPos) == false) {   //  Not a valid base, reset the counter.
      _kmerLoad = 0;
    }

    else {
      _fmer.addR(_buffer[_bufferPos]);   //  A valid base, so push it onto
      _rmer.addL(_buffer[_bufferPos]);   //  the kmer.

      if (_kmerLoad < _kmerSize)         //  Increment the loaded size,
        _kmerLoad++;                     //  if not full already.
    }

    _bufferPos++;

    return(true);                        //  A base was consumed.
  }


  kmerTiny   fmer(void)      { return(_fmer);                        };
  kmerTiny   rmer(void)      { return(_rmer);                        };
  uint64     position(void)  { return(_bufferPos - _kmerSize);       };

  uint64     bgnPosition(void)  { return(_bufferPos - _kmerSize);    };
  uint64     endPosition(void)  { return(_bufferPos);                };

private:
  uint32       _kmerSize;
  uint32       _kmerLoad;
  uint32       _kmerValid;

  char const  *_buffer;
  uint64       _bufferLen;
  uint64       _bufferPos;

  kmerTiny     _fmer;
  kmerTiny     _rmer;
};


#endif  //  MERYL_UTIL_KMER_ITERATOR_H
