
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

void
showSplit(char *str, splitType type) {
  splitToWords  W;

  W.split(str, type);

  fprintf(stderr, "SPLIT: '%s'\n", str);

  for (uint32 ii=0; ii<W.numWords(); ii++)
    fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
}


void
showSplit(char  sep, char *str) {
  splitToWords  W;

  W.split(str, sep);

  fprintf(stderr, "SPLIT: '%s'  (on letter '%c')\n", str, sep);

  for (uint32 ii=0; ii<W.numWords(); ii++)
    fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
}


void
showSplit(char *sep, char *str) {
  splitToWords  W;

  W.split(str, sep);

  fprintf(stderr, "SPLIT: '%s'  (on letters '%s')\n", str, sep);

  for (uint32 ii=0; ii<W.numWords(); ii++)
    fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
}



void
showKeyValue(char const *str) {
  KeyAndValue  kv;

  if (kv.find(str) == true) {
    fprintf(stderr, "KV '%s'\n", str);
    fprintf(stderr, "K  '%s'\n", kv.key());
    fprintf(stderr, "V  '%s'\n", kv.value());
  } else {
    fprintf(stderr, "KV '%s' - no key=value pair found.\n", str);
  }
}



bool
checkKeyValue(char const *line, char const *key, char const *val, bool pass) {
  KeyAndValue  kv;
  bool         success = kv.find(line);

  if (success != pass) {
    fprintf(stderr, "ERROR:  line '%s' returned '%s' but should have returned '%s'\n",
            line,
            (success) ? "true" : "false",
            (pass)    ? "true" : "false");

    return(false);
  }

  if (success == false)   //  If KeyAndValue failed, the test was successful (unless
    return(true);         //  it was supposed to succeed, then we'd fail above).

  if (strcmp(kv.key(), key) != 0) {
    fprintf(stderr, "ERROR:  line '%s' failed to match key '%s'; returned '%s'\n", line, key, kv.key());
    success = false;
  }

  if (strcmp(kv.value(), val) != 0) {
    fprintf(stderr, "ERROR:  line '%s' failed to match value '%s'; returned '%s'\n", line, val, kv.value());
    success = false;
  }

  if (success == false)
    exit(1);

  //if (success == false) {
  //  fprintf(stderr, "\n");
  //  fprintf(stderr, "KeyValue test failed.\n");
  //  exit(1);
  //}

  return(success);
}



bool
checkKeyValue(void) {
  KeyAndValue  kv;
  bool         pass = true;

  pass &= checkKeyValue(nullptr,                        "",  "",   false);
  pass &= checkKeyValue("",                             "",  "",   false);

  pass &= checkKeyValue("=",                            "",  "",   false);
  pass &= checkKeyValue("#",                            "",  "",   false);
  pass &= checkKeyValue(" =",                           "",  "",   false);
  pass &= checkKeyValue("  #",                          "",  "",   false);

  pass &= checkKeyValue("a===b",                        "a", "b",  true);
  pass &= checkKeyValue("a::: b",                       "a", "b",  true);
  pass &= checkKeyValue("a =:= b",                      "a", "b",  true);

  pass &= checkKeyValue("a = =b",                       "a", "=b", true);

  pass &= checkKeyValue("a=",                           "a", "",   true);
  pass &= checkKeyValue("a= # comment",                 "a", "",   true);
  pass &= checkKeyValue(" a=",                          "a", "",   true);
  pass &= checkKeyValue(" a= ",                         "a", "",   true);
  pass &= checkKeyValue(" a=\r\n",                      "a", "",   true);
  pass &= checkKeyValue(" a= #comment",                 "a", "",   true);
  pass &= checkKeyValue("a = ! comment",                "a", "",   true);
  pass &= checkKeyValue("a  = # comment",               "a", "",   true);
  pass &= checkKeyValue("  a = ! comment",              "a", "",   true);
  pass &= checkKeyValue("  a  = # comment",             "a", "",   true);

  pass &= checkKeyValue("a=b",                          "a", "b",  true);
  pass &= checkKeyValue("a=b    # comment",             "a", "b",  true);
  pass &= checkKeyValue(" a=b ",                        "a", "b",  true);
  pass &= checkKeyValue(" a=b   # comment",             "a", "b",  true);
  pass &= checkKeyValue("  a=b  ",                      "a", "b",  true);
  pass &= checkKeyValue("  a=b  # comment",             "a", "b",  true);

  pass &= checkKeyValue("a=#b",                         "a", "#b", true);
  pass &= checkKeyValue("a=#b # comment",               "a", "#b", true);
  pass &= checkKeyValue(" a=#b ",                       "a", "#b", true);
  pass &= checkKeyValue(" a=#b # comment",              "a", "#b", true);
  pass &= checkKeyValue("  a=#b  ",                     "a", "#b", true);
  pass &= checkKeyValue("  a=#b  # comment",            "a", "#b", true);

  pass &= checkKeyValue("a=b#",                         "a", "b#", true);
  pass &= checkKeyValue("a=b# # comment # # ",          "a", "b#", true);
  pass &= checkKeyValue(" a=b# ",                       "a", "b#", true);
  pass &= checkKeyValue(" a=b#  # comment # #",         "a", "b#", true);
  pass &= checkKeyValue("  a=b#  ",                     "a", "b#", true);
  pass &= checkKeyValue("  a=b#    # comment #      ",  "a", "b#", true);

  pass &= checkKeyValue("a=!b",                         "a", "!b", true);
  pass &= checkKeyValue("a=!b !comment",                "a", "!b", true);
  pass &= checkKeyValue(" a=!b ",                       "a", "!b", true);
  pass &= checkKeyValue(" a=!b   !comment # !",         "a", "!b", true);
  pass &= checkKeyValue("  a=!b  ",                     "a", "!b", true);
  pass &= checkKeyValue("  a=!b  !  comment",           "a", "!b", true);

  pass &= checkKeyValue(" a = b ",                      "a", "b",  true);
  pass &= checkKeyValue("  a  =  b  ",                  "a", "b",  true);
  pass &= checkKeyValue("   a   =   b",                 "a", "b",  true);

  pass &= checkKeyValue(" a = b ",                      "a", "b",  true);
  pass &= checkKeyValue(" a = b # comment",             "a", "b",  true);
  pass &= checkKeyValue("  a  =  b  ",                  "a", "b",  true);
  pass &= checkKeyValue("  a  =  b  ! comment",         "a", "b",  true);
  pass &= checkKeyValue("   a   =   b",                 "a", "b",  true);
  pass &= checkKeyValue("   a   =   b#",                "a", "b#", true);
  pass &= checkKeyValue("   a   =   #b#",               "a", "",   true);
  pass &= checkKeyValue("   a   =   #b#    ",           "a", "",   true);
  pass &= checkKeyValue("   a   =   #b#    # comment",  "a", "",   true);

  return(pass);
}



int
main(int argc, char **argv) {
  bool  pass = true;

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-p") == 0) {
      showSplit(argv[++arg], splitPaths);
    }

    else if (strcmp(argv[arg], "-w") == 0) {
      showSplit(argv[++arg], splitWords);
    }

    else if (strcmp(argv[arg], "-c") == 0) {
      showSplit(argv[arg+1][0], argv[arg+2]);
      arg += 2;
    }

    else if (strcmp(argv[arg], "-s") == 0) {
      showSplit(argv[arg+1], argv[arg+2]);
      arg += 2;
    }

    else if (strcmp(argv[arg], "-k") == 0) {
      showKeyValue(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-K") == 0) {
      bool  p = checkKeyValue();

      if (p == true)
        fprintf(stderr, "KeyValue tests passed.\n");
      else
        fprintf(stderr, "KeyValue tests FAILED.\n");

      pass &= p;
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }


  if ((argc == 1) || (err.size() > 0)) {
    fprintf(stderr, "usage: %s [...]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, " -p /path/to/split      - show splitToWords operating on a path\n");
    fprintf(stderr, " -w 'string to split'   - show splitToWords operating on a string\n");
    fprintf(stderr, " -c n splitnatnletter   - show splitToWords operating on a specific letter\n");
    fprintf(stderr, " -s [] split[at]letters - show splitToWrords operating on specific letters\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " -k key=value         - show KeyAndValue operating on a key=value pair\n");
    fprintf(stderr, "\n");
    fprintf(stderr, " -K                   - run tests on keyAndValue\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    return(1);
  }

  return(pass == false);
}
