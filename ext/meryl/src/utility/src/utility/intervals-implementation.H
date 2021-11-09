
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

#include "arrays.H"
#include <algorithm>


#ifndef INTERVALS_IMPLEMENTATION
#error Include intervals.H instead of intervals-implementation.H
#else


template <class iNum>
void
intervals<iNum>::add_position(iNum bgn, iNum end) {

  if (bgn > end)
    fprintf(stderr, "intervals<iNum>::add_position()-- ERROR: bgn=%u > end=%u\n", bgn, end);
  assert(bgn <= end);

  if (_listMax == 0)
    allocateArray(_list, _listMax, 32);

  increaseArray(_list, _listLen, _listMax, _listMax / 4);

  _list[_listLen]._bgn = bgn;
  _list[_listLen]._end = end;
  _list[_listLen]._num = 1;

  _listLen++;

  _isSorted   = false;
  _isSquashed = false;
}



template <class iNum>
void
intervals<iNum>::add(intervals<iNum> const &that) {

  resizeArray(_list, _listLen, _listMax, _listLen + that._listLen);

  for (uint32 ii=0; ii<that._listLen; ii++, _listLen++)
    _list[_listLen] = that._list[ii];

  _isSorted   = false;
  _isSquashed = false;
}



template <class iNum>
void
intervals<iNum>::remove(uint32 idx) {

  assert(idx < _listLen);

  for (uint32 ii=idx; ii+1<_listLen; ii++)
    _list[ii] = _list[ii+1];

  _listLen--;
}



template <class iNum>
void
intervals<iNum>::sort(void) {

  if ((_isSorted == true) ||
      (_listLen < 2))
    return;

  auto increasing = [](_ir const &a,
                       _ir const &b) {
                      return(((a._bgn  < b._bgn)) ||
                             ((a._bgn == b._bgn) && (a._end  < b._end)));
                    };

  std::sort(_list, _list + _listLen, increasing);

  _isSorted = true;
}



template <class iNum>
void
intervals<iNum>::squash(iNum minOverlap) {
  uint32  intoI = 0;  //  Interval we're merging into.
  uint32  fromI = 1;  //  Interval we're merging from.

  if (_isSquashed == true)
    return;

  sort();

  while (fromI < _listLen) {
    assert(_list[intoI]._bgn <  _list[intoI]._end);   //  Basic checks.  Both intervals
    assert(_list[fromI]._bgn <  _list[fromI]._end);   //  cannot be empty, and intoI
    assert(_list[intoI]._bgn <= _list[fromI]._bgn);   //  must be before fromI.

    //  If the fromI intersects with intoI -- either contained in intoI, or
    //  has a thick overlap to intoI -- merge it in.  We're guaranteed that
    //  this._bgn is before next._bgn, so all we need to do is extend
    //  this._end to cover the next interval.

    if ((_list[intoI]._end >= _list[fromI]._end) ||
        (_list[intoI]._end >= _list[fromI]._bgn + minOverlap)) {
      _list[intoI]._end  = std::max(_list[fromI]._end, _list[intoI]._end);
      _list[intoI]._num +=          _list[fromI]._num;
    }

    //  Otherwise, move to the next intoI, copy the current fromI to it, and
    //  then move to the next fromI.  We should, to be pedantic, check that
    //  intoI != fromI before the copy, but no harm if we don't.

    else {
      _list[++intoI] = _list[fromI];
    }

    fromI++;
  }

  _listLen    = intoI + 1;   //  Update the length of the list,
  _isSquashed = true;        //  and note that it's now merged.
}



template <class iNum>
void
intervals<iNum>::filter(iNum minLength, iNum maxLength) {
  uint32  intoI = 0;
  uint32  fromI = 0;

  //  Over every interval, if it is long enough, copy it
  //  into the 'new' list.

  while (fromI < _listLen) {
    iNum  length = _list[fromI]._end - _list[fromI]._bgn;

    if ((minLength <= length) &&
        (length <= maxLength))
      _list[intoI++] = _list[fromI];

    fromI++;
  }

  _listLen = intoI;
}




#if 0
template <class iNum>
void
setToUnion(intervals<iNum> const &A,
           intervals<iNum> const &B) {
}


template <class iNum>
void
setToIntersection(intervals<iNum> const &A,
                  intervals<iNum> const &B) {
}


template <class iNum>
void
setToContained(intervals<iNum> const &A,
               intervals<iNum> const &B) {
}


template <class iNum>
void
setToUnion(iNum bgn, iNum end,
           intervals<iNum> const &A) {
}


template <class iNum>
void
setToIntersection(iNum bgn, iNum end,
                  intervals<iNum> const &A) {
}


template <class iNum>
void
setToContained(iNum bgn, iNum end,
               intervals<iNum> const &A) {
}
#endif


//  Helper function to invert a squashed intervals list.
template <class iNum>
void
intervals<iNum>::setToInversion1(iNum bgn, iNum end,
                                      intervals<iNum> const &A) {

  delete [] _list;

  _listLen = 0;                    //  Create a new list to store the
  _listMax = A._listLen + 1;       //  inversion.  We need at most one
  _list    = new _ir [_listMax];   //  more interval than the original.

  //  If no existing list, just add a single interval covering the universe.
  //
  //  If the inversion range falls entirely inside a gap in the original list
  //  (which would also result in the inverted list having one interval
  //  covering the whole range) we'll catch it in the loop below.

  if (A._listLen == 0) {
    _list[_listLen++] = { bgn, end, 1 };
  }

  //  For an existing list:
  //    1) Add an interval for the first gap, if it's inside the inersion
  //       range.
  //    2) Add intervals covering the middle gaps.  Threshold each endpoint
  //       by the inversion range, and only add a new interval if it is of
  //       positive length.
  //    3) Add an interval for the last gap, if it's inside the inversion
  //       range.

  else {
    if (bgn < A._list[0]._bgn)
      _list[_listLen++] = { bgn, A._list[0]._bgn, 1 };

    for (uint32 ii=1; ii<A._listLen; ii++) {
      iNum  nb = std::max(bgn, A._list[ii-1]._end);
      iNum  ne = std::min(end, A._list[ii  ]._bgn);

      if (nb < ne)
        _list[_listLen++] = { nb, ne, 1 };
    }

    if (A._list[A._listLen-1]._end < end)
      _list[_listLen++] = { A._list[A._listLen-1]._end, end, 1 };
  }

  //  Check that we didn't blow up.

  assert(_listLen <= _listMax);
}


//  Helper function to invert a non-squashed intervals list.
template <class iNum>
void
intervals<iNum>::setToInversion2(iNum bgn, iNum end,
                                      intervals<iNum> const &A) {

  delete [] _list;

  _listLen = 0;                    //  Create a new list to store the
  _listMax = A._listLen * 2;       //  inversion.  We need at most twice
  _list    = new _ir [_listMax];   //  the original size.

  //  If no existing list, just add a single interval covering the universe.

  if (A._listLen == 0) {
    _list[_listLen++] = { bgn, end, 1 };
  }

  //  For an existing list:
  //    Add two intervals for each existing interval, one on each end of the
  //    interval.  The new intervals are thresholded aginst the inversion
  //    range, and only added if they are of positive length.
  //
  //  Note the symmetrically-opposite comparisons; these prevent us from
  //  adding length=0 intervals.

  else {
    iNum   nb, ne;

    for (uint32 ii=0; ii<A._listLen; ii++) {
      if ((bgn < A._list[ii]._bgn) && (A._list[ii]._bgn <= end))
        _list[_listLen++] = { bgn, A._list[ii]._bgn, 1 };

      if ((bgn <= A._list[ii]._end) && (A._list[ii]._end < end))
        _list[_listLen++] = { A._list[ii]._end, end, 1 };
    }
  }

  //  Check that we didn't blow up.

  assert(_listLen <= _listMax);
}


template <class iNum>
void
intervals<iNum>::setToInversion(iNum bgn, iNum end,
                                     intervals<iNum> const &A) {
  if (A._isSquashed)
    setToInversion1(bgn, end, A);
  else
    setToInversion2(bgn, end, A);
}



template <class iNum>
void
intervalsDepth<iNum>::computeDepth(intervals<iNum> const &IL) {
  uint32    idplen = IL.size() * 2;
  _idp     *idp    = new _idp [idplen];

  for (uint32 ii=0; ii<IL.size(); ii++) {
    idp[2*ii  ]._pos = IL.bgn(ii);     //  Enter into an inteval, change
    idp[2*ii  ]._dlt = 1;              //  depth by +1.

    idp[2*ii+1]._pos = IL.end(ii);     //  Leave an interval, change
    idp[2*ii+1]._dlt = -1;             //  depth by -1.
  }

  delete [] _list;

  _listLen  = 0;
  _list     = nullptr;

  if (idplen > 0)
    computeDepth(idplen, idp);

  delete [] idp;
}



template <class iNum>
void
intervalsDepth<iNum>::computeDepth(uint32 idplen, _idp *idp) {

  //  Sort regions so that earlier positions are first, and so that depth
  //  increases (+1) are before decreases (-1).

  auto increasing = [](_idp const &a,
                       _idp const &b) {
                      return(((a._pos  < b._pos)) ||
                             ((a._pos == b._pos) && (a._dlt > b._dlt)));
                    };

  std::sort(idp, idp + idplen, increasing);

  //  The first thing must be an 'open' event.  If not, someone supplied a
  //  negative length to the original intervalList.  Or, possibly, two
  //  zero-length intervals.

  if (idp[0]._dlt == -1)
    for (uint32 ii=0; ii<idplen; ii++)
      fprintf(stderr, "idp[%u] pos %d dlt %d\n", ii, idp[ii]._pos, idp[ii]._dlt);

  assert(idp[0]._dlt == 1);

  //  Init first interval.

  _listLen  = 0;
  _list     = new _idr [idplen + 1];

  _list[_listLen]._bgn = idp[0]._pos;
  _list[_listLen]._end = idp[0]._pos;
  _list[_listLen]._dpt = 1;

  for (uint32 i=1; i<idplen; i++) {

    //  Update the end of the current interval to this position.

    _list[_listLen]._end = idp[i]._pos;

    //  If this position is different than the last position, make
    //  a new depth interval.

    if (idp[i-1]._pos != idp[i]._pos) {
      _listLen++;

      _list[_listLen]._bgn = idp[i]._pos;
      _list[_listLen]._end = idp[i]._pos;
      _list[_listLen]._dpt = _list[_listLen-1]._dpt;
    }

    //  Process any depth change associated with this position.

    _list[_listLen]._dpt += idp[i]._dlt;

    //  Is it safe to blindly change the depth of this region?  Yes.  If the
    //  position is different than the last, we've already made a new depth
    //  region.  And if it wasn't different, one of the previous positions
    //  must have made a new depth region.  In any case, the depth region
    //  we're currently at must always be length 0 -- we're never able to
    //  change the depth of a region when we're at the end coordinate.

    assert(_list[_listLen]._bgn == _list[_listLen]._end);
  }

  assert(_listLen < idplen + 1);
  assert(_list[_listLen]._bgn == _list[_listLen]._end);
  assert(_list[_listLen]._dpt == 0);
}

#endif   //  INTERVALS_IMPLEMENTATION
