
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

#include "strings.H"



int
main(int argc, char **argv) {
  std::vector<char const *> *errors = nullptr;

  if (argc == 1) {
    fprintf(stderr, "usage: %s [-e] pattern string string ...\n", argv[0]);
    fprintf(stderr, "  -e  use vector error handler\n");
    return 1;
  }

  int arg = 1;
  if (strcmp(argv[arg], "-e") == 0) {
    errors = new std::vector<char const *>;
    arg++;
  }

  merylutil::RegEx  *re = new merylutil::RegEx;

  if (re->compile(argv[arg++], errors) == false) {
    for (uint32 ii=0; errors != nullptr && ii<errors->size(); ii++)
      fprintf(stderr, "ERROR: %s\n", (*errors)[ii]);
    return 1;
  }

  for (; arg < argc; arg++) {
    if (re->match(argv[arg], errors) == false) {
      fprintf(stderr, "Failed to match string '%s'.\n", argv[arg]);
      for (uint32 ii=0; errors != nullptr && ii<errors->size(); ii++)
        fprintf(stderr, "ERROR: %s\n", (*errors)[ii]);
      return 1;
    }

    else {
      for (uint32 ii=0; ii<re->matchesLen(); ii++)
        fprintf(stderr, "MATCH: '%s'\n", re->getMatch(ii));
    }
  }

  delete re;

  return 0;
}
