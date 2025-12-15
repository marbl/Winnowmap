
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

#ifdef CANU

void
merylInput::openInputSeqStore(void) {
  _store            = new sqStore(_storeName);

  _sqBgn            = 1;                                  //  C-style, not the usual
  _sqEnd            = _store->sqStore_lastReadID() + 1;   //  sqStore semantics!

  if (_segMax > 1) {
    uint64  nBases = 0;

    for (uint32 ss=1; ss <= _store->sqStore_lastReadID(); ss++)
      nBases += _store->sqStore_getReadLength(ss);

    uint64  nBasesPerSeg = nBases / _segMax;

    _sqBgn = 0;
    _sqEnd = 0;

    nBases = 0;

    for (uint32 ss=1; ss <= _store->sqStore_lastReadID(); ss++) {
      nBases += _store->sqStore_getReadLength(ss);

      if ((_sqBgn == 0) && ((nBases / nBasesPerSeg) == _seg - 1))
        _sqBgn = ss;

      if ((_sqEnd == 0) && ((nBases / nBasesPerSeg) == _seg))
        _sqEnd = ss;
    }

    if (_seg == _segMax)                           //  Annoying special case; if the last segment,
      _sqEnd = _store->sqStore_lastReadID() + 1;   //  sqEnd is set to the last read, not N+1.

    fprintf(stderr, "merylInput-- segment %u/%u picked reads %u-%u out of %u\n",
            _seg, _segMax, _sqBgn, _sqEnd, _store->sqStore_lastReadID());
  }

  _read        = new sqRead;
  _readID      = _sqBgn - 1;       //  Incremented before loading the first read
  _readPos     = uint32max;
}


bool
merylInput::loadBasesFromCanu(char    *seq,
                              uint64   maxLength,
                              uint64  &seqLength,
                              bool    &endOfSequence) {

  //  If no read currently loaded, load one, or return that we're done.
  //  We need to loop so we can ignore the length zero reads in seqStore
  //  that exist after correction/trimming.

  while (_readPos >= _read->sqRead_length()) {
    _readID++;

    if (_readID >= _sqEnd)  //  C-style iteration, not usual sqStore semantics.
      return(false);

    _store->sqStore_getRead(_readID, _read);
    _readPos = 0;
  }

  //  How much of the read is left to return?

  uint32  len = _read->sqRead_length() - _readPos;

  assert(len > 0);

  //  If the output space is big enough to hold the rest of the read, copy it,
  //  flagging it as the end of a sequence, and setup to load the next read.

  if (len <= maxLength) {
    memcpy(seq, _read->sqRead_sequence() + _readPos, sizeof(char) * len);

    _readPos       = _read->sqRead_length();

    seqLength      = len;
    endOfSequence  = true;
  }

  //  Otherwise, only part of the data will fit in the output space.

  else {
    memcpy(seq, _read->sqRead_sequence() + _readPos, sizeof(char) * maxLength);

    _readPos      += maxLength;

    seqLength      = maxLength;
    endOfSequence  = false;
  }

  return(true);
}

#endif
