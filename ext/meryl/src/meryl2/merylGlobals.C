
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

#include "meryl.H"

merylGlobals  globals;



bool
merylGlobals::processDebugOption(int &arg, char **argv,
                                 std::vector<char const *> &err) {

  if (strcmp(argv[arg], "dumpIndex") == 0) {         //  Report the index for the dataset.
    arg++;                                           //  It's just the parameters used for encoding.
    delete new merylFileReader(argv[arg++], true);   //  Expects a meryl db directory as a parameter.
    exit(0);
  }

  if (strcmp(argv[arg], "dumpFile") == 0) {          //  Dump the index for a single data file.
    arg++;                                           //  Expects a meryl file prefix as a parameter.
    dumpMerylDataFile(argv[arg++]);                  //  (e.g., db.meryl/0x000000)
    exit(0);
  }

  return false;
}


bool
merylGlobals::processGlobalOption(int &arg, char **argv,
                                  std::vector<char const *> &err) {

  if (strcmp(argv[arg], "-k") == 0) {
    kmerTiny::setSize(strtouint32(argv[++arg]));
    return true;
  }

  if (strncmp(argv[arg], "k=", 2) == 0) {
    fprintf(stderr, "WARNING: obsolete '%s' supplied; use '-k %s' instead.\n",
            argv[arg], argv[arg]+2);
    kmerTiny::setSize(strtouint32(argv[arg]+2));
    return true;
  }

  if ((strcmp(argv[arg],  "-c")      == 0) ||
      (strcmp(argv[arg], "--hpc") == 0) ||
      (strcmp(argv[arg], "--homopoly-compress") == 0)) {
    sprintf(err, "Global '%s' not yet supported; use 'compress' instead.", argv[arg]);
    //_doCompression = true;
    return true;
  }

  if (strcmp(argv[arg], "--ssr") == 0) {
    sprintf(err, "Global '%s %s' not yet supported; only local homopoly compression with 'compress' supported.", argv[arg], argv[arg+1]);
    arg++;
    return true;
  }

  if (strcmp(argv[arg], "-l") == 0) {
    kmerTiny::setLabelSize(strtouint32(argv[++arg]));
    return true;
  }

  if ((strcmp(argv[arg],  "-m")      == 0) ||
      (strcmp(argv[arg],  "-memory") == 0) ||
      (strcmp(argv[arg], "--memory") == 0)) {
    allowedMemory() = getAllowedMemory(argv[++arg], err);
    return true;
  }


  if (strncmp(argv[arg], "memory=", 7) == 0) {
    fprintf(stderr, "WARNING: obsolete '%s' supplied; use '-m %s' instead.\n",
            argv[arg], argv[arg]+7);
    allowedMemory() = getAllowedMemory(argv[arg]+7, err);
    return true;
  }

  if ((strcmp(argv[arg],  "-t")       == 0) ||
      (strcmp(argv[arg],  "-threads") == 0) ||
      (strcmp(argv[arg], "--threads") == 0)) {
    allowedThreads() = getAllowedThreads(argv[++arg], err);
    return true;
  }

  if (strncmp(argv[arg], "threads=", 8) == 0) {
    fprintf(stderr, "WARNING: obsolete '%s' supplied; use '-t %s' instead.\n",
            argv[arg], argv[arg]+8);
    allowedThreads() = getAllowedThreads(argv[arg]+8, err);
    return true;
  }

  if (strncmp(argv[arg], "-V", 2) == 0) {          //  Anything that starts with -V
    for (uint32 vv=1; vv<strlen(argv[arg]); vv++)  //  increases verbosity by the
      increaseVerbosity();                         //  number of letters.
    return true;
  }

  if (strcmp(argv[arg], "-Q") == 0) {
    beQuiet();
    return true;
  }

  if (strcmp(argv[arg], "-P") == 0) {
    enableProgressReport();
    return true;
  }

  if (strcmp(argv[arg], "-C") == 0) {
    enableConfigureOnly();
    return true;
  }

  if ((strcmp(argv[arg], "-h")     == 0) ||
      (strcmp(argv[arg], "-help")  == 0) ||
      (strcmp(argv[arg], "--help") == 0)) {
    err.push_back(nullptr);
    return true;
  }

  return false;
}
