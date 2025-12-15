
/******************************************************************************
 *
 *  This file is part of meryl, a genomic k-kmer counter with nice features.
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

#include "system.H"

#include "matchToken.H"

int
main(int argc, char **argv) {
  char const *optout = nullptr;

  //  True cases: <opt> is a prefix of <pat> and the separator matches.

  assert(matchToken("output:",  optout, "output:") == true);   assert(optout != nullptr);   assert(optout[0] == 0);
  assert(matchToken("out:",     optout, "output:") == true);   assert(optout != nullptr);   assert(optout[0] == 0);
  assert(matchToken("o:",       optout, "output:") == true);   assert(optout != nullptr);   assert(optout[0] == 0);

  assert(matchToken("output:b", optout, "output:") == true);   assert(optout != nullptr);   assert(strcmp(optout, "b") == 0);
  assert(matchToken("out:b",    optout, "output:") == true);   assert(optout != nullptr);   assert(strcmp(optout, "b") == 0);
  assert(matchToken("o:b",      optout, "output:") == true);   assert(optout != nullptr);   assert(strcmp(optout, "b") == 0);

  assert(matchToken("output=",  optout, "output=") == true);   assert(optout != nullptr);   assert(optout[0] == 0);
  assert(matchToken("out=",     optout, "output=") == true);   assert(optout != nullptr);   assert(optout[0] == 0);
  assert(matchToken("o=",       optout, "output=") == true);   assert(optout != nullptr);   assert(optout[0] == 0);

  assert(matchToken("output=b", optout, "output=") == true);   assert(optout != nullptr);   assert(strcmp(optout, "b") == 0);
  assert(matchToken("out=b",    optout, "output=") == true);   assert(optout != nullptr);   assert(strcmp(optout, "b") == 0);
  assert(matchToken("o=b",      optout, "output=") == true);   assert(optout != nullptr);   assert(strcmp(optout, "b") == 0);

  //  True cases: opt does not end with a separator.

  assert(matchToken("output", optout, "output:") == true);   assert(optout == nullptr);
  assert(matchToken("out",    optout, "output:") == true);   assert(optout == nullptr);
  assert(matchToken("o",      optout, "output:") == true);   assert(optout == nullptr);

  assert(matchToken("output", optout, "output=") == true);   assert(optout == nullptr);
  assert(matchToken("out",    optout, "output=") == true);   assert(optout == nullptr);
  assert(matchToken("o",      optout, "output=") == true);   assert(optout == nullptr);

  assert(matchToken("output", optout, "output")  == true);   assert(optout == nullptr);
  assert(matchToken("out",    optout, "output")  == true);   assert(optout == nullptr);
  assert(matchToken("o",      optout, "output")  == true);   assert(optout == nullptr);

  //  False cases: opt ends with a separator and pat does not.

  assert(matchToken("output:",  optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("out:",     optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("o:",       optout, "output") == false);   assert(optout == nullptr);

  assert(matchToken("output:b", optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("out:b",    optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("o:b",      optout, "output") == false);   assert(optout == nullptr);

  assert(matchToken("output=",  optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("out=",     optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("o=",       optout, "output") == false);   assert(optout == nullptr);

  assert(matchToken("output=b", optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("out=b",    optout, "output") == false);   assert(optout == nullptr);
  assert(matchToken("o=b",      optout, "output") == false);   assert(optout == nullptr);

  //  False cases: opt is not a prefix of pat.

  assert(matchToken("outwit",   optout, "output:") == false);   assert(optout == nullptr);
  assert(matchToken("outwit:",  optout, "output:") == false);   assert(optout == nullptr);
  assert(matchToken("outwit=",  optout, "output:") == false);   assert(optout == nullptr);
  assert(matchToken("outwit:a", optout, "output:") == false);   assert(optout == nullptr);
  assert(matchToken("outwit=a", optout, "output:") == false);   assert(optout == nullptr);

  assert(matchToken("outwit",   optout, "output=") == false);   assert(optout == nullptr);
  assert(matchToken("outwit:",  optout, "output=") == false);   assert(optout == nullptr);
  assert(matchToken("outwit=",  optout, "output=") == false);   assert(optout == nullptr);
  assert(matchToken("outwit:a", optout, "output=") == false);   assert(optout == nullptr);
  assert(matchToken("outwit=a", optout, "output=") == false);   assert(optout == nullptr);

  assert(matchToken("outwit",   optout, "output")  == false);   assert(optout == nullptr);
  assert(matchToken("outwit:",  optout, "output")  == false);   assert(optout == nullptr);
  assert(matchToken("outwit=",  optout, "output")  == false);   assert(optout == nullptr);
  assert(matchToken("outwit:a", optout, "output")  == false);   assert(optout == nullptr);
  assert(matchToken("outwit=a", optout, "output")  == false);   assert(optout == nullptr);

  //  Exact cases.

  assert(matchToken("output",    optout, "output:", true) == true);   assert(optout == nullptr);
  assert(matchToken("output:",   optout, "output:", true) == true);   assert(optout != nullptr);   assert(optout[0] == 0);
  assert(matchToken("output:b",  optout, "output:", true) == true);   assert(optout != nullptr);   assert(optout[0] == 'b');
  assert(matchToken("out",       optout, "output:", true) == false);
  assert(matchToken("out:b",     optout, "output:", true) == false);

  assert(matchToken("output",    optout, "output",  true) == true);    assert(optout == nullptr);
  assert(matchToken("output:b",  optout, "output",  true) == false);   assert(optout == nullptr);
  assert(matchToken("out",       optout, "output",  true) == false);   assert(optout == nullptr);
  assert(matchToken("out:b",     optout, "output",  true) == false);   assert(optout == nullptr);

  assert(matchToken("outwit:b",  optout, "output:", true) == false);   assert(optout == nullptr);
  assert(matchToken("outwit",    optout, "output:", true) == false);   assert(optout == nullptr);

  assert(matchToken("outwit:b",  optout, "output",  true) == false);   assert(optout == nullptr);
  assert(matchToken("outwit",    optout, "output",  true) == false);   assert(optout == nullptr);

  //  Let the user try to break it.

  if (argc == 1) {
    fprintf(stderr, "usage: %s <option>[:=]<=value> <pattern>[:=]\n", argv[0]);
    fprintf(stderr, "  Tests if 'option' matches 'pattern', reports the matching pieces.\n");
    fprintf(stderr, "  (Also tests a bunch of built-in examples to ensure\n");
    fprintf(stderr, "   basic sanity, and those tests HAVE passed.)\n");
    return 0;
  }

  if (matchToken(argv[1], optout, argv[2]) == false) {
    fprintf(stderr, "NEGATIVE. '%s' doesn't match '%s'.\n", argv[1], argv[2]);
    return 1;
  }

  fprintf(stderr, "POSITIVE! '%s' matches '%s' with optout '%s'.\n",
          argv[1],
          argv[2],
          (optout) ? optout : "<nullptr>");

  return 0;
}
