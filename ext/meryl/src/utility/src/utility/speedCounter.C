
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

#include "speedCounter.H"

const char*
speedCounter::_spinr[4] = { "[|]", "[/]", "[-]", "[\\]" };

const char*
speedCounter::_liner[19] = { "[-         ]",
                             "[--        ]",
                             "[ --       ]",
                             "[  --      ]",
                             "[   --     ]",
                             "[    --    ]",
                             "[     --   ]",
                             "[      --  ]",
                             "[       -- ]",
                             "[        --]",
                             "[         -]",
                             "[        --]",
                             "[       -- ]",
                             "[      --  ]",
                             "[     --   ]",
                             "[    --    ]",
                             "[   --     ]",
                             "[  --      ]",
                             "[ --       ]" };


speedCounter::speedCounter(char const   *fmt,
                           double        unit,
                           uint64        freq,
                           bool          enabled) {
  _count     = 0;
  _draws     = 0;
  _unit      = unit;
  _freq      = freq;
  _startTime = getTime();
  _fmt       = fmt;
  _spin      = false;
  _line      = false;
  _enabled   = enabled;

  //  We use _draws instead of shifting _count just because it's
  //  simpler, and both methods need another variable anyway.

  //  Set all the bits below the hightest set in _freq --
  //  this allows us to do a super-fast test in tick().
  //
  _freq |= _freq >> 1;
  _freq |= _freq >> 2;
  _freq |= _freq >> 4;
  _freq |= _freq >> 8;
  _freq |= _freq >> 16;
  _freq |= _freq >> 32;
}

speedCounter::~speedCounter() {
  finish();
}
