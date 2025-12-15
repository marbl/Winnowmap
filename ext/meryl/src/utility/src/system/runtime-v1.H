
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
#include "version.H"

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


//  A helper for, usually, processing command line errors.
//  Lines longer than 1024 are truncated.
//
int sprintf(std::vector<char const *> &ev, char const *fmt, ...);
int sprintf(std::vector<char const *> &ev, char const *fmt, va_list ap);

int sprintf(std::vector<char const *> *ev, char const *fmt, ...);
int sprintf(std::vector<char const *> *ev, char const *fmt, va_list ap);


//  A helper.
//
//  If 'fatal', print message to stderr and exit(1).
//  Else, return false.
//
//  To be used in functions that can switch from returning fatal errors
//  (which was the default behavior) to returning true/false.
//
//    bool doSomething(bool crash-on-failure) {
//      [...]
//      if (something-failed)
//        return fatalError(crash-on-failure, "no good");
//
//  If crash-on-failure is true, we print 'no good' and exit the program.
//  If not true, we (silently) return false.
//
//  The vector variant behaves the same if no vector is supplied.
//
//  If a vector is supplied, the error message is appended to it
//  and false is always returned.
//
//  fmt is expected to end WITHOUT a new line.  One will be added
//  when printed to the screen / when the error list is reported.
//
bool fatalError(bool fatal,                                char const *fmt, ...);
bool fatalError(bool fatal, std::vector<char const *> *ev, char const *fmt, ...);
bool fatalError(            std::vector<char const *> *ev, char const *fmt, ...);


//  Return true of an external command is exists and is executable.
//
//  This will popen("cmd > /dev/null 2>&1") and return true if both
//  the popen() and pclose() succeed, which will occur if the command
//  can be executed and it returns a success error code.
//
//    'lsd' cannot be executed.
//    'ls /...' will (usually) return failure as file '/...' doesn't (usually) exist.
//    'ls -d /bin' returns success as directory '/bin' (usually) exists.
//
bool
commandAvailable(char const *cmd);

#endif
