
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

#include "types.H"



void
pps(char const *str, uint64 bgn, uint64 end, char const *val, std::vector<char const *> &log) {
  uint64  ee = (end == 0) ? strlen(str) : end;

  fprintf(stdout, "'%s'\n", str);

  for (uint64 ii=0; ii<=bgn; ii++)
    fprintf(stdout, " ");

  for (uint64 ii=bgn; ii<ee; ii++)   //  To handle 'end=0' being strlen
    fprintf(stdout, "^");

  if (log.empty() == true)
    fprintf(stdout, "\n --> %s\n", val);
  else {
    fprintf(stdout, "\n --> FAIL!\n");

    for (uint32 ii=0; ii<log.size(); ii++)
      if (log[ii] != NULL)
        fprintf(stderr, " --> %s\n", log[ii]);
  }

  log.clear();
}



void
selftest(int32 testNum, int32 desired, int32 bgn, int32 end, bool shouldPass) {
  std::vector<char const *>  log;
  int32                      i32 = 0;

  int32 result = decodeInteger("d-4ki10hg45dif", bgn, end, i32, log);

  if ((shouldPass == true) && (result == desired) && (log.empty() == true))
    return;

  if ((shouldPass == false) && (log.empty() == false)) {
    for (uint32 ii=0; ii<log.size(); ii++)   //  Just to make valgrind happy.
      delete [] log[ii];
    return;
  }

  fprintf(stderr, "TEST %d failed:  result %d != desired %d\n", testNum, result, desired);

  for (uint32 ii=0; ii<log.size(); ii++) {
    if (log[ii] != NULL)
      fprintf(stderr, "  %s\n", log[ii]);
    delete [] log[ii];
  }
}



int
main(int argc, char **argv) {
  bool    test = false;
  bool    pass = true;
  uint64  bgn  = 0;
  uint64  end  = 0;

  std::vector<char const *>  log;
  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-b") == 0) {
      bgn = strtouint64(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-e") == 0) {
      end = strtouint64(argv[++arg]);
    }

    else if (strcmp(argv[arg], "-s8") == 0) {
      int8  val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-u8") == 0) {
      uint8 val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-s16") == 0) {
      int16  val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-u16") == 0) {
      uint16 val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-s32") == 0) {
      int32  val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-u32") == 0) {
      uint32 val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-s64") == 0) {
      int64  val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-u64") == 0) {
      uint32 val = 88;

      decodeInteger(argv[++arg], bgn, end, val, log);
      pps(argv[arg], bgn, end, toDec(val), log);
    }

    else if (strcmp(argv[arg], "-t") == 0) {
      test = true;
    }

    else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  if ((err.size() > 0) || (argc == 1)) {
    fprintf(stderr, "usage: %s ..\n", argv[0]);
    fprintf(stderr, "  -b bgn    decode starting at position bgn\n");
    fprintf(stderr, "  -e end    decode ending at position end (strlen if end=0)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -s8  str  decode string as  int8\n");
    fprintf(stderr, "  -u8  str  decode string as uint8\n");
    fprintf(stderr, "  -s16 str  decode string as  int16\n");
    fprintf(stderr, "  -u16 str  decode string as uint16\n");
    fprintf(stderr, "  -s32 str  decode string as  int32\n");
    fprintf(stderr, "  -u32 str  decode string as uint32\n");
    fprintf(stderr, "  -s64 str  decode string as  int64\n");
    fprintf(stderr, "  -u64 str  decode string as uint64\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t        run self tests\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    return(1);
  }

  if (test) {                            //  Test against string "d-4ki10hg45dif"
    selftest(1,     0,  0,  1,  true);   //                       ^      .        - 'd'
    selftest(2,    -4,  1,  3,  true);   //                        ^^    .        - '-4'
    selftest(3, -4000,  1,  4,  true);   //                        ^^^   .        - '-4k'
    selftest(4, -4096,  1,  5,  true);   //                        ^^^^  .        - '-4ki'
    selftest(5,  4096,  2,  5,  true);   //                         ^^^  .        - '4ki'
    selftest(6,    16,  5,  8,  true);   //                            ^^^        - '10h'
    selftest(7,    16,  5,  9, false);   //                            ^^^^       - '10hg' should fail
    selftest(8,    45,  9, 12,  true);   //                                ^^^    - '45d'
    selftest(9,    45,  9, 13,  true);   //                                ^^^^   - '45di'
  }

  if (log.size() > 0) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Decoding Errors:\n");

    for (uint32 ii=0; ii<log.size(); ii++)
      if (log[ii] != NULL)
        fprintf(stderr, "  %s\n", log[ii]);
  }


  return(pass == false);
}


