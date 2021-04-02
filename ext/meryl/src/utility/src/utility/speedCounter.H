
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

#ifndef SPEEDCOUNTER_H
#define SPEEDCOUNTER_H

#include "types.H"
#include "system.H"

class speedCounter {
public:
  //  fmt specifies the status format.  An example:
  //    "    %8f [unit]things (%8.5f [unit]things/sec)\r"
  //
  speedCounter(char const   *fmt,
               double        unit,
               uint64        freq,
               bool          enabled=true);
  ~speedCounter();

  void   enableSpinner(void) { _spin = true; };
  void   enableLiner(void)   { _line = true; };

  bool   tick(void) {
    if (_enabled && ((++_count & _freq) == 0)) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      _draws++;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
      return(true);
    }
    return(false);
  };

  bool   tick(uint64 increment) {
    if (_enabled == false)
      return(false);

    _count += increment;
    if ((_count & _freq) == 0) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      _draws++;
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fflush(stderr);
      return(true);
    }
    return(false);
  };

  void   finish(void) {
    if (_enabled && (_count >= _freq)) {
      double  v = _count / _unit;
      if (_spin) fputs(_spinr[_draws %  4], stderr);
      if (_line) fputs(_liner[_draws % 19], stderr);
      fprintf(stderr, _fmt, v, v / (getTime() - _startTime));
      fprintf(stderr, "\n");
      fflush(stderr);
    }
    _count = 0;
  };

private:
  static const char  *_spinr[4];
  static const char  *_liner[19];
  uint64              _count;
  uint64              _draws;
  double              _unit;
  uint64              _freq;
  double              _startTime;
  char const         *_fmt;
  bool                _spin;
  bool                _line;
  bool                _enabled;
};


#endif  //  SPEEDCOUNTER_H
