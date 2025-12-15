
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

#ifndef TIME_H
#define TIME_H

#include "types.H"

#include <sys/stat.h>  //  for UTIME_OMIT

#include <time.h>      //  localtime_r

namespace merylutil::inline system::inline v1 {

//
//  The meryl-utility time objects.
//
//  muTime stores a single time value and will convert to/from printable
//  formats (eventually).
//
//  muDeltaTime stores two time values and reports the difference between
//  them.
//
//  muFileTime stores two time values, one for the modification time of a
//  file, and one for the last access time (but this isn't always supported).
//


class muTime {
public:
  muTime(int64 s, int64 n) {
    _s = s;
    _n = n;
  }
  muTime(bool init=false) {
    if (init == true)
      setTime();
  }

  muTime &setTime(void);

  double  getTime_double(void)     { return (double)_s + (double)_n / 1000000000.0;  }
  int64   getTime_int64 (void)     { return         _s;                              }

  uint64  getTime_seconds(void)    { return         _s;                              }
  uint64  getTime_nanosec(void)    { return         _n;                              }

  uint64  getTime_billionths(void) { return         _s << 32 | _n;                   } 

  //char const *getTimeString(void);

  void    getTime(uint32 &yr, uint32 &mo, uint32 &da, uint32 &hh, uint32 &mm, uint32 &ss) {
    time_t     tt = _s + ((_n < 500000000) ? 0 : 1);
    struct tm  tm;

    localtime_r(&tt, &tm);

    yr = tm.tm_year + 1900;
    mo = tm.tm_mon  + 1;
    da = tm.tm_mday;
    hh = tm.tm_hour;
    mm = tm.tm_min;
    ss = tm.tm_sec;
  }

public:
  bool    operator< (muTime &that) { return (_s <  that._s) || (_s == that._s && _n <  that._n); }
  bool    operator> (muTime &that) { return (_s >  that._s) || (_s == that._s && _n >  that._n); }
  bool    operator<=(muTime &that) { return (_s <= that._s) || (_s == that._s && _n <= that._n); }
  bool    operator>=(muTime &that) { return (_s >= that._s) || (_s == that._s && _n >= that._n); }

  bool    operator==(muTime &that) { return (_s == that._s) && (_n == that._n); }
  bool    operator!=(muTime &that) { return (_s != that._s) || (_n != that._n); }

#if 0
  muTime &operator=(const muTime &that) {
    _s = that._s;
    _n = that._n;
    return *this;
  }
#endif

private:
  int64  _s = 0;
  int64  _n = 0;

  friend class muDeltaTime;
  friend class muFileTime;
};



class muDeltaTime {
public:
  muDeltaTime(bool init=false) {
    if (init == true)
      _bgn = _end.setTime();
  }

public:
  void      startTimer(void)  {  _bgn.setTime();  }
  void      stopTimer(void)   {  _end.setTime();  }

  muTime   &beginTime(void)   {  return _bgn;  }
  muTime   &endTime(void)     {  return _end;  }

public:
  int64     years(void)   { if (_yrs == int64max) updateYDHMS(); return _yrs; }    //  These return the integer number of (unit)s
  int64     days(void)    { if (_yrs == int64max) updateYDHMS(); return _dys; };   //  after subtracting out the larger units.
  int64     hours(void)   { if (_yrs == int64max) updateYDHMS(); return _hrs; };   //  Example: for a time of 181 seconds:
  int64     minutes(void) { if (_yrs == int64max) updateYDHMS(); return _mns; };   //    minutes() = 3
  double    seconds(void) { if (_yrs == int64max) updateYDHMS(); return _scs; };   //    seconds() = 1.0

public:
  double    full_years(void)   { return ttt() / _years;   }   //  These return the absolue quantity of
  double    full_days(void)    { return ttt() / _days;    }   //  each unit.  Example: for 181 seconds:
  double    full_hours(void)   { return ttt() / _hours;   }   //    full_hours()   = 0.0502777
  double    full_minutes(void) { return ttt() / _minutes; }   //    full_minutes() = 3.01666
  double    full_seconds(void) { return ttt() / _seconds; }   //    full_seconds() = 181.0

private:
  double    ttt(void) {
    return _end.getTime_double() - _bgn.getTime_double();
  }

  void      updateYDHMS(void) {
    double t = ttt();

    _yrs = (int64)floor(t / _years);    t -= _yrs * _years;
    _dys = (int64)floor(t / _days);     t -= _dys * _days;
    _hrs = (int64)floor(t / _hours);    t -= _yrs * _hours;
    _mns = (int64)floor(t / _minutes);  t -= _yrs * _minutes;
    _scs =              t;
  }

private:
  muTime  _bgn;
  muTime  _end;

  static constexpr double _seconds =   1.0;
  static constexpr double _minutes =  60.0 * _seconds;
  static constexpr double _hours   =  60.0 * _minutes;
  static constexpr double _days    =  24.0 * _hours;
  static constexpr double _years   = 365.0 * _days;

  int64   _yrs = int64max;
  int64   _dys = int64max;
  int64   _hrs = int64max;
  int64   _mns = int64max;
  double  _scs = DBL_MAX;
};





//  setTimeOfFile()
//   - Will set the modification time to the time in the object.  If
//     getTimeOfFile() was used to set the object time, then the access
//     time is also set.
//
//  getTimeOfFile()
//   - Will set the object time to the modification time of the file.
//     The access time is also saved.
//
class muFileTime {
public:
  muFileTime() {
    _mod = muTime(0, UTIME_OMIT);
    _acc = muTime(0, UTIME_OMIT);
  }
  muFileTime(char const *prefix, char separator='.', char const *suffix=nullptr) {
    getTimeOfFile(prefix, separator, suffix);
  }
  muFileTime(FILE *F) {
    getTimeOfFile(F);
  }

  bool      setTimeOfFile(char const *prefix, char separator='.', char const *suffix=nullptr);
  muTime    getTimeOfFile(char const *prefix, char separator='.', char const *suffix=nullptr);
  muTime    getTimeOfFile(FILE *F);

  muTime   &lastModifyTime(void)  { return _mod; }
  muTime   &lastAccessTime(void)  { return _acc; }

  muDeltaTime  fileAge(void) {
    muDeltaTime  dt;

    dt.beginTime() = _mod;
    dt.endTime().setTime();

    return dt;
  }

  bool    operator< (muFileTime &that) { return (_mod <  that._mod); }
  bool    operator> (muFileTime &that) { return (_mod >  that._mod); }
  bool    operator<=(muFileTime &that) { return (_mod <= that._mod); }
  bool    operator>=(muFileTime &that) { return (_mod >= that._mod); }

  bool    operator==(muFileTime &that) { return (_mod == that._mod); }
  bool    operator!=(muFileTime &that) { return (_mod != that._mod); }

private:
  muTime  _mod;
  muTime  _acc;
};


}  //  namespace merylutil::system::v1

#endif
