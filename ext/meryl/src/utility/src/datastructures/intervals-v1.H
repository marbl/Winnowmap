
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-20
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef MERYLUTIL_DATASTRUCTURES_INTERVALS_V1_H
#define MERYLUTIL_DATASTRUCTURES_INTERVALS_V1_H

#include "types.H"

namespace merylutil::inline files::inline v1 {

//  The interval coordinates use the usual C semantics of [bgn, end) -
//  'x=bgn' is inside the interval, but 'x=end' is not.

template <class iNum>
class intervals {
private:
  struct _ir {
    iNum      _bgn;
    iNum      _end;
    uint32    _num;
  };

public:
  intervals()    {                  };
  ~intervals()   { delete [] _list; };

  void      clear(void) {
    _isSorted   = true;
    _isSquashed = true;
    _listLen    = 0;
  };

  //  Accessors.

  uint32    size(void) const         { return(_listLen); };

  iNum      bgn (uint32 idx) const   { return(_list[idx]._bgn); };
  iNum      end (uint32 idx) const   { return(_list[idx]._end); };
  iNum      span(uint32 idx) const   { return(_list[idx]._end - _list[idx]._bgn); };

  uint32    count(uint32 idx) const  { return(_list[idx]._num); };

  //  Modifiers.

  iNum     &bgn (uint32 idx)         { return(_list[idx]._bgn); };
  iNum     &end (uint32 idx)         { return(_list[idx]._end); };

  uint32   &count(uint32 idx)        { return(_list[idx]._num); };

  void      clear(uint32 idx) {
    _list[idx]._bgn = iNum();
    _list[idx]._end = iNum();
    _list[idx]._num = 0;
  }

  //  Creation.
  //
  //  Add a single interval to the list of intervals specified by either
  //   - the position of the end points
  //   - the position of the start and the length of the span
  //
  //  Add(intervals) will copy all the intervals from B into this object,
  //  no further processing (sorting, squashing or filtering) is performed.
  //
  //  Remove the interval at position 'idx' in our list.  Doing so will
  //  greatly screw up interation over the intervals, and it is suggested
  //  to instead change the span of the interval to zero and then filter
  //  them out after iteration is complete.

  void      add_position(iNum bgn, iNum end);
  void      add_span    (iNum bgn, iNum len) {
    if (len < 0)
      add_position(bgn+len, bgn);
    else
      add_position(bgn, bgn+len);
  };

  void      add(intervals<iNum> const &B);

  void      remove(uint32 idx);

  //  Sort intervals by increasing coordinate, breaking ties with the end
  //  coordinate.
  //
  //  Combine intervals that overlap by at least 'minOverlap' into one item.
  //
  //  Discard intervals that are smaller than minLength or larger than
  //  maxLength.

  void      sort(void);
  void      squash(iNum minOverlap=0);
  void      filter(iNum minLength, iNum maxLength);

  //  setToUnion - populate this intervals object with all the intervals in A
  //  and B.  If both A and B are squashed, this intervals object will also
  //  be squashed.
  //
  //  setToIntersection - each interval in A (B) is intersected with all
  //  intervals in B (A), and the resulting interval is added to this object.
  //
  //  setToContained - each interval in A that is contained fully in some
  //  interval in B is added to this intervals object.
#if 0
  void      setToUnion       (intervals<iNum> const &A, intervals<iNum> const &B);
  void      setToIntersection(intervals<iNum> const &A, intervals<iNum> const &B);
  void      setToContained   (intervals<iNum> const &A, intervals<iNum> const &B);
#endif
  //  setToUnion - copy the intervals in A that oveerlap with the interval
  //  bgn-end.
  //
  //  setToIntersection - copy the intervals in A that intersect with the
  //  interval bgn-end, and trim them to that range.
  //
  //  setToContained - copy the intervals in A that are contained within the
  //  interval bgn-end.
  //
  //  setToInversion
  //   - if A is squashed, intervals that fill the 'holes' in A, bounded by
  //     bgn and end) are added to this object.
  //   - if A is not squashed, each interval in A will contribute 0, 1 or 2
  //     new intervals to this object, representing the holes, bounded by bgn and end,
  //     created by only that single interval in A.
  //
  //                   bgn[               ]end
  //                --------  ---------     ----  A
  //                --------  ---------           union
  //                      --  ---------           intersection
  //                          ---------           contained
  //                        --         ----       inversion
#if 0
  void      setToUnion       (iNum bgn, iNum end, intervals<iNum> const &A);
  void      setToIntersection(iNum bgn, iNum end, intervals<iNum> const &A);
  void      setToContained   (iNum bgn, iNum end, intervals<iNum> const &A);
#endif
  void      setToInversion   (iNum bgn, iNum end, intervals<iNum> const &A);

  //  Helper functions.
private:
  void      setToInversion1(iNum bgn, iNum end, intervals<iNum> const &A);
  void      setToInversion2(iNum bgn, iNum end, intervals<iNum> const &A);

private:
  bool     _isSorted   = true;
  bool     _isSquashed = true;

  uint32   _listMax    = 0;
  uint32   _listLen    = 0;
  _ir     *_list       = nullptr;
};



template <class iNum>
class intervalsDepth {
private:
  struct _idp {     //  An intervalDepthPosition stores the position
    iNum   _pos;    //  of a change in depth, and the delta of that
    int32  _dlt;    //  change (which is either +1 or -1).
  };

  struct _idr {     //  An intervalDepthRegion has the coordinates
    iNum   _bgn;    //  of the region and the depth.
    iNum   _end;
    uint32 _dpt;
  };

public:
  intervalsDepth() {
  };
  intervalsDepth(intervals<iNum> const &IL) {
    computeDepth(IL);
  };
  ~intervalsDepth() {
    delete [] _list;
  };

  uint32    size(void)          { return(_listLen); };

  iNum      bgn (uint32 idx)    { return(_list[idx]._bgn); };
  iNum      end (uint32 idx)    { return(_list[idx]._end); };
  iNum      span(uint32 idx)    { return(_list[idx]._end - _list[idx]._bgn); };

  uint32    depth(uint32 idx)   { return(_list[idx]._dpt); };

  void      computeDepth(intervals<iNum> const &IL);

private:
  void      computeDepth(uint32 idplen, _idp *idp);

  uint32   _listLen = 0;
  _idr    *_list    = nullptr;
};




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

}  //  merylutil::files::v1

#endif  //  MERYLUTIL_DATASTRUCTURES_INTERVALS_V1_H
