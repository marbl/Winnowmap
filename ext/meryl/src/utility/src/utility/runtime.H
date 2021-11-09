
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

//  This is the global include file that all C files in the AS
//  subsystem should include.

#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include "types.H"

//  AS_configure() will set some run-time stuff like FPU modes and initialize
//  accounting (like time since the start of execution).
//
//  Most importantly, it will query the environment to decide how many
//  threads to enable.  The maxThreads parameter can be used to restrict this
//  to something _lower_.  To set the number of threads to something
//  arbitrary, call AS_configure() then setNumThreads().
//
//  Lastly, AS_configure() will handle the --version switch.
//
int AS_configure(int argc, char **argv, uint32 maxThreads=uint32max);

#endif
