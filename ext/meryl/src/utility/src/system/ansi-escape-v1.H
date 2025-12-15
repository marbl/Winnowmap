
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

#ifndef ANSI_ESCAPE_H
#define ANSI_ESCAPE_H

#include "types.H"

namespace merylutil::inline system::inline v1 {

enum class ansiCode {
  Normal           =   0,
  Bold             =   1,
  Faint            =   2,
  Underline        =   4,
  Blink            =   5,
  Inverse          =   7,
  CrossedOut       =   9,

  Black            =  30,
  Red              =  31,
  Green            =  32,
  Yellow           =  33,
  Blue             =  34,
  Magenta          =  35,
  Cyan             =  36,
  White            =  37,

  BlackBG          =  40,
  RedBG            =  41,
  GreenBG          =  42,
  YellowBG         =  43,
  BlueBG           =  44,
  MagentaBG        =  45,
  CyanBG           =  46,
  WhiteBG          =  47,

  BrightBlack      =  90,
  BrightRed        =  91,
  BrightGreen      =  92,
  BrightYellow     =  93,
  BrightBlue       =  94,
  BrightMagenta    =  95,
  BrightCyan       =  96,
  BrightWhite      =  97,

  BrightBlackBG    = 100,
  BrightRedBG      = 101,
  BrightGreenBG    = 102,
  BrightYellowBG   = 103,
  BrightBlueBG     = 104,
  BrightMagentaBG  = 105,
  BrightCyanBG     = 106,
  BrightWhiteBG    = 107,

  END              = 255
};


//  Build a string to set the above modes.
//  Can be called as:
//    makeAnsiEscapeSequence( (const ansiCode[]){ ansiCode::Green, ansiCode::Bold, ansiCode::END } );
//    makeAnsiEscapeSequence( (const ansiCode[]){ ansiCode::Normal,                ansiCode::END } );
//
char const *
makeAnsiEscapeSequence(const ansiCode *codes);

}

#endif  //  ANSI_ESCAPE_H
