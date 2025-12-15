
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

using merylutil::splitToWords;
using merylutil::KeyAndValue;
using merylutil::stringList;
using merylutil::trimString;



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

  return(success);
}



bool
checkKeyValue(void) {
  KeyAndValue  kv;
  bool         pass = true;

  fprintf(stderr, "Testing 'KeyAndValue'.\n");

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



bool
trimTest(char *s, char const *r) {
  char   o[16];
  bool   p = true;

  memcpy(o, s, 16);

  trimString(s);

  if (memcmp(s, r, 16) != 0) {
    fprintf(stderr, "FAILED to trim '%s'\n", o);
    fprintf(stderr, "          into '%s'\n", r);
    fprintf(stderr, "           got '%s' instead.\n", s);
    p = false;
  }

  return(p);
}



bool
trimTest(void) {
  char s[16] = {0};
  char r[16] = {0};
  bool pass = true;

  fprintf(stderr, "Testing 'trimString'.\n");

  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "               "); strcpy(r, "");      r[ 0] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "AaaaaaaaaaaaaaB"); memcpy(r, s, 16);   r[15] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "A             B"); memcpy(r, s, 16);   r[15] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, " AbbbbbbbbbbbB "); memcpy(r, s+1, 14); r[13] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, " A           B "); memcpy(r, s+1, 14); r[13] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "  AcccccccccB  "); memcpy(r, s+2, 12); r[11] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "  A         B  "); memcpy(r, s+2, 12); r[11] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "AB             "); strcpy(r, "AB");    r[ 3] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "     AB        "); strcpy(r, "AB");    r[ 3] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "             AB"); strcpy(r, "AB");    r[ 3] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "A              "); strcpy(r, "A");     r[ 2] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "      A        "); strcpy(r, "A");     r[ 2] = 0;  pass &= trimTest(s, r);
  memset(s, 0, 16); memset(r, 0, 16); strcpy(s, "              A"); strcpy(r, "A");     r[ 2] = 0;  pass &= trimTest(s, r);

  return(pass);
}



int
main(int argc, char **argv) {
  splitToWords    W;
  KeyAndValue     kv;
  bool            pass = true;


  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {

    //
    //  Manual tests.
    //

    if      (strcmp(argv[arg], "-p") == 0) {
      fprintf(stderr, "SPLIT: '%s'  (as paths)\n", argv[arg+1]);
      W.split(argv[arg+1], merylutil::splitPaths);
      for (uint32 ii=0; ii<W.numWords(); ii++)
        fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
      arg+=1;
    }

    else if (strcmp(argv[arg], "-w") == 0) {
      fprintf(stderr, "SPLIT: '%s'  (as words)\n", argv[arg+1]);
      W.split(argv[arg+1], merylutil::splitWords);
      for (uint32 ii=0; ii<W.numWords(); ii++)
        fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
      arg+=1;
    }

    else if (strcmp(argv[arg], "-c") == 0) {
      fprintf(stderr, "SPLIT: '%s'  (on letter '%c')\n", argv[arg+2], argv[arg+1][0]);
      W.split(argv[arg+2], argv[arg+1][0]);
      for (uint32 ii=0; ii<W.numWords(); ii++)
        fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
      arg+=2;
    }

    else if (strcmp(argv[arg], "-s") == 0) {
      fprintf(stderr, "SPLIT: '%s'  (on letters '%s')\n", argv[arg+2], argv[arg+1]);
      W.split(argv[arg+2], argv[arg+1]);
      for (uint32 ii=0; ii<W.numWords(); ii++)
        fprintf(stderr, "%02u - '%s'\n", ii, W[ii]);
      arg+=2;
    }

    else if (strcmp(argv[arg], "-k") == 0) {
      if (kv.find(argv[arg+1]) == true) {
        fprintf(stderr, "KV '%s'\n", argv[arg+1]);
        fprintf(stderr, "K  '%s'\n", kv.key());
        fprintf(stderr, "V  '%s'\n", kv.value());
      } else {
        fprintf(stderr, "KV '%s' - no key=value pair found.\n", argv[arg+1]);
      }
      arg+=1;
    }

    else if (strcmp(argv[arg], "-f") == 0) {
      stringList sl(argv[arg+1]);
      for (uint32 ii=0; sl[ii]; ii++)
        fprintf(stdout, "%s\n", sl[ii]);
      arg+=1;
    }

    //
    //  Automated tests.
    //

    else if (strcmp(argv[arg], "-K") == 0) {
      pass &= checkKeyValue();
    }

    else if (strcmp(argv[arg], "-t") == 0) {
      pass &= trimTest();
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
    fprintf(stderr, "  Tests that require user inspection:\n");
    fprintf(stderr, "    -p /path/to/split      - show splitToWords operating on a path\n");
    fprintf(stderr, "    -w 'string to split'   - show splitToWords operating on a string\n");
    fprintf(stderr, "    -c n splitnatnletter   - show splitToWords operating on a specific letter\n");
    fprintf(stderr, "    -s [] split[at]letters - show splitToWrords operating on specific letters\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -f <file>              - use stringList to print a file to stdout\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -k key=value           - show KeyAndValue operating on a key=value pair\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Automatic tests:\n");
    fprintf(stderr, "  -K                     - run tests on keyAndValue\n");
    fprintf(stderr, "  -t                     - test white-space trimming\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    return(1);
  }

  if (pass)
    fprintf(stdout, "Pass!\n");

  return(pass == false);
}
