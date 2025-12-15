
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

#ifndef MERYLUTIL_DS_SETS_V1_H
#define MERYLUTIL_DS_SETS_V1_H

#include "types.H"

#include "bits/murmur-v1.H"


namespace merylutil::inline sets::inline v1 {


//  A simple fixed-size ordered set that allows fast comparisons of inequality.
//
//  For sorted inputs, insert() is linear.  For unsorted inputs, insert() is
//  n(n+1)/2 worst-case.
//
//  Assumes, for now, that T is uint32.
//
//  computeSignature() MUST be called after all words are inserted.
//
template<typename T>
class linearset {
public:
  linearset(uint64 maxWords) {
    _lMax = maxWords;
    _l    = new T [_lMax];
  }
  linearset(linearset const &that) {
    _sig = that._sig;
    duplicateArray(     _l,      _lLen,      _lMax,
                   that._l, that._lLen, that._lMax);
  }

  ~linearset() {
    delete [] _l;
  }

  uint64   size(void)              { return _lLen; }

  T        get       (uint64 idx)  { return _l[idx]; }
  T        operator[](uint64 idx)  { return _l[idx]; }

  void     insert(T w) {
    assert(_lLen < _lMax);

    if ((_lLen == 0) || (_l[_lLen-1] < w)) {           //  If empty or w is the largest
      _l[_lLen++] = w;                                 //  element, add to the end of the.
    }                                                  //  list.  Else, over all elements,
    else {                                             //  if w is smaller, swap it into
      for (uint32 ll=0; ll<_lLen; ll++)                //  the list; if equal, do nothing (it's
        if      (w  < _l[ll])   std::swap(w, _l[ll]);  //  already in the set).  Then add
        else if (w == _l[ll])   return;                //  whatever is in w to the end of the list.

      _l[_lLen++] = w;
    }
  }

  void computeSignature(void)  {
    if (_sig == uint64zero)
      _sig = mumurmur(_l, _lLen);
  }

  bool operator==(linearset const &that) const {
    if (_sig  != that._sig)    return false;
    if (_lLen != that._lLen)   return false;

    for (uint32 ll=0; ll<_lLen; ll++)
      if (_l[ll] != that._l[ll])
        return false;

    return true;
  }

private:
  uint64  _sig  = uint64zero;

  uint64  _lMax = 0;
  uint64  _lLen = 0;
  T      *_l    = nullptr;
};


//  A set of linear sets.
template<typename T>
class linearsets {
public:
  linearsets() {
  }

  ~linearsets() {
    for (uint64 ii=0; ii<_sLen; ii++)
      delete _s[ii];

    delete [] _s;
  }

  uint64         size(void)              { return _sLen; }

  linearset<T>  *get       (uint64 idx)  { assert(idx < _sLen);  return _s[idx]; }
  linearset<T>  *operator[](uint64 idx)  { assert(idx < _sLen);  return _s[idx]; }

  uint64         find  (linearset<T> &s) {   //  Return the index of 's' if
    s.computeSignature();                    //  it exists, otherwise, uint64max.
    for (uint64 ii=0; ii<_sLen; ii++)
      if (s == *_s[ii])
        return ii;
    return uint64max;
  };

  uint64         insert(linearset<T> &s) {   //  Return the index of 's' if
    s.computeSignature();                    //  it is in the set, otherwise, add
    for (uint64 ii=0; ii<_sLen; ii++)        //  it and and return the index.
      if (s == *_s[ii])
        return ii;
    return appendToArray(new linearset(s), _s, _sLen, _sMax, 256);
#if 0
    increaseArray(_s, _sLen, _sMax, 256);

    _s[_sLen++] = new linearset(s);

    return _sLen-1;
#endif
  };

private:
  uint64          _sMax = 0;
  uint64          _sLen = 0;
  linearset<T>  **_s    = nullptr;
};



}  //  merylutil::sets::v1

#endif  //  MERYLUTIL_DS_SETS_V1_H
