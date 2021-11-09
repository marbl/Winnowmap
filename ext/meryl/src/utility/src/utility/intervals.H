
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

#ifndef INTERVALS_H
#define INTERVALS_H

#include "types.H"

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



#define INTERVALS_IMPLEMENTATION
#include "intervals-implementation.H"
#undef  INTERVALS_IMPLEMENTATION


#endif  //  INTERVALS_H
