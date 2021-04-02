
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

#include "meryl-lookup.H"

void   dumpExistence(lookupGlobal *G);
void   reportExistence(lookupGlobal *G);
void   filter(lookupGlobal *G);


void
lookupGlobal::initialize(void) {
  omp_set_num_threads(nThreads);
}



void
lookupGlobal::loadLookupTables(void) {
  std::vector<merylFileReader *>    merylDBs;    //  Input meryl database.
  std::vector<double>               minMem;      //  Estimated min memory for lookup table.
  std::vector<double>               optMem;      //  Estimated max memory for lookup table.

  //  Open input meryl databases, initialize lookup.

  for (uint32 ii=0; ii<lookupDBname.size(); ii++) {
    merylDBs .push_back(new merylFileReader(lookupDBname[ii]));
    minMem   .push_back(0.0);
    optMem   .push_back(0.0);
    lookupDBs.push_back(new merylExactLookup());
  }

  //  Estimate memory needed for each lookup table.

  double   minMemTotal = 0.0;
  double   optMemTotal = 0.0;

  for (uint32 ii=0; ii<lookupDBname.size(); ii++) {
    fprintf(stderr, "--\n");
    fprintf(stderr, "-- Estimating memory usage for '%s'.\n", lookupDBname[ii]);
    fprintf(stderr, "--\n");

    double  minm, optm;
    lookupDBs[ii]->estimateMemoryUsage(merylDBs[ii], maxMemory, minm, optm, minV, maxV);

    minMemTotal += minm;
    optMemTotal += optm;
  }

  //  Use either the smallest or 'fastest' table, or fail, depending on how
  //  much memory the use lets us use.

  bool  useOpt = (optMemTotal <= maxMemory);
  bool  useMin = (minMemTotal <= maxMemory) && (useOpt == false);

  fprintf(stderr, "--\n");
  fprintf(stderr, "-- Minimal memory needed: %.3f GB%s\n", minMemTotal, (useMin) ? "  enabled" : "");
  fprintf(stderr, "-- Optimal memory needed: %.3f GB%s\n", optMemTotal, (useOpt) ? "  enabled" : "");
  fprintf(stderr, "-- Memory limit           %.3f GB\n",   maxMemory);
  fprintf(stderr, "--\n");

  if ((useMin == false) &&
      (useOpt == false)) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Not enough memory to load databases.  Increase -memory.\n");
    fprintf(stderr, "\n");
    exit(1);
  }

  if (doEstimate == true) {
    fprintf(stderr, "-- Stopping after memory estimated reported; -estimate option enabled.\n");
    exit(0);
  }

  //  Now load the data and forget about the input databases.

  for (uint32 ii=0; ii<lookupDBname.size(); ii++) {
    fprintf(stderr, "--\n");
    fprintf(stderr, "-- Loading kmers from '%s' into lookup table.\n", lookupDBname[ii]);
    fprintf(stderr, "--\n");

    if (lookupDBs[ii]->load(merylDBs[ii], maxMemory, useMin, useOpt, minV, maxV) == false)
      exit(1);

    delete merylDBs[ii];
  }
}



//  Open input sequences.
void
lookupGlobal::openInputs(void) {

  if (seqName1) {
    fprintf(stderr, "-- Opening input sequences '%s'.\n", seqName1);
    seqFile1 = new dnaSeqFile(seqName1);
  }

  if (seqName2) {
    fprintf(stderr, "-- Opening input sequences '%s'.\n", seqName2);
    seqFile2 = new dnaSeqFile(seqName2);
  }
}



//  Open output writers.
void
lookupGlobal::openOutputs(void) {

  if (outName1) {
    fprintf(stderr, "-- Opening output file '%s'.\n", outName1);
    outFile1 = new compressedFileWriter(outName1);
  }

  if (outName2) {
    fprintf(stderr, "-- Opening output file '%s'.\n", outName1);
    outFile2 = new compressedFileWriter(outName2);
  }
}






int
main(int argc, char **argv) {
  lookupGlobal  *G = new lookupGlobal;

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if        (strcmp(argv[arg], "-sequence") == 0) {
      G->seqName1 = argv[++arg];

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->seqName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-mers") == 0) {
      while ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->lookupDBname.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-labels") == 0) {
      while ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->lookupDBlabel.push_back(argv[++arg]);

    } else if (strcmp(argv[arg], "-output") == 0) {
      G->outName1 = argv[++arg];

      if ((arg + 1 < argc) && (argv[arg + 1][0] != '-'))
        G->outName2 = argv[++arg];

    } else if (strcmp(argv[arg], "-min") == 0) {
      G->minV = (kmvalu)strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-max") == 0) {
      G->maxV = (kmvalu)strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-threads") == 0) {
      G->nThreads = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-memory") == 0) {
      G->maxMemory = strtodouble(argv[++arg]);

    } else if (strcmp(argv[arg], "-dump") == 0) {
      G->reportType = lookupOp::opDump;

    } else if (strcmp(argv[arg], "-existence") == 0) {
      G->reportType = lookupOp::opExistence;

    } else if (strcmp(argv[arg], "-include") == 0) {
      G->reportType = lookupOp::opInclude;

    } else if (strcmp(argv[arg], "-exclude") == 0) {
      G->reportType = lookupOp::opExclude;

    } else if (strcmp(argv[arg], "-estimate") == 0) {
      G->doEstimate = true;

    } else if (strcmp(argv[arg], "-V") == 0) {
      G->showProgress = true;

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }
  }

  //  Check for invalid usage.

  if (G->reportType == lookupOp::opNone) {
    err.push_back("No report-type (-existence, -dump, -include, -exclude) supplied.\n");
  }

  if (G->reportType == lookupOp::opDump) {
    if (G->seqName1 == nullptr)  err.push_back("No input sequences (-sequence) supplied.\n");
    if (G->seqName2 != nullptr)  err.push_back("Only one input sequence (-sequence) supported for -dump.\n");

    if (G->outName1 == nullptr)  err.push_back("No output file (-output) supplied.\n");
    if (G->outName2 != nullptr)  err.push_back("Only one output file (-output) supported for -dump.\n");

    if (G->lookupDBname.size() == 0) err.push_back("No meryl database (-mers) supplied.\n");
    if (G->lookupDBname.size()  > 1) err.push_back("Only one meryl database (-mers) supported for -dump.\n");
  }

  if (G->reportType == lookupOp::opExistence) {
    if (G->seqName1 == nullptr)  err.push_back("No input sequences (-sequence) supplied.\n");
    if (G->seqName2 != nullptr)  err.push_back("Only one input sequence (-sequence) supported for -existence.\n");

    if (G->outName1 == nullptr)  err.push_back("No output file (-output) supplied.\n");
    if (G->outName2 != nullptr)  err.push_back("Only one output file (-output) supported for -existence.\n");

    if (G->lookupDBname.size() == 0) err.push_back("No meryl database (-mers) supplied.\n");
  }

  if ((G->reportType == lookupOp::opInclude) ||
      (G->reportType == lookupOp::opExclude)) {
    if (G->seqName1 == nullptr)  err.push_back("No input sequences (-sequence) supplied.\n");
    if (G->outName1 == nullptr)  err.push_back("No output file (-output) supplied.\n");

    if ((G->seqName2 != nullptr) &&
        (G->outName2 == nullptr)) err.push_back("No second output file (-output) supplied for second input (-input) file.\n");

    if ((G->seqName2 == nullptr) &&
        (G->outName2 != nullptr)) err.push_back("No second input file (-input) supplied for second output (-output) file.\n");

    if (G->lookupDBname.size() == 0) err.push_back("No meryl database (-mers) supplied.\n");
    if (G->lookupDBname.size()  > 1) err.push_back("Only one meryl database (-mers) supported for -include or -exclude.\n");
  }


  if (err.size() > 0) {
    fprintf(stderr, "usage: %s <report-type> \\\n", argv[0]);
    fprintf(stderr, "        [-estimate] \\\n");
    fprintf(stderr, "         -sequence <input1.fasta> [<input2.fasta>] \\\n");
    fprintf(stderr, "         -output   <output1>      [<output2>]\n");
    fprintf(stderr, "         -mers     <input1.meryl> [<input2.meryl>] [...] \\\n");
    fprintf(stderr, "         -labels   <input1name>   [<input2name>]   [...]\n");
    fprintf(stderr, "  Query the kmers in meryl database(s) <input.meryl> with the sequences\n");
    fprintf(stderr, "  in <input.fasta>.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Multiple databases are supported.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Up to two input sequences are supported (only for -include / -exclude).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Input files can be FASTA or FASTQ; uncompressed, gz, bz2 or xz compressed\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Output from each input is sent to the associated output file.  Files will be\n");
    fprintf(stderr, "  compressed if the appropriate extension is supplied (gz, bz2 or xz).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Each input database can be filtered by value.  More advanced filtering\n");
    fprintf(stderr, "  requires a new database to be constructed using meryl.\n");
    fprintf(stderr, "    -min   m    Ignore kmers with value below m\n");
    fprintf(stderr, "    -max   m    Ignore kmers with value above m\n");
    fprintf(stderr, "    -threads t  Number of threads to use when constructing lookup table.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Memory usage can be limited, within reason, by sacrificing kmer lookup\n");
    fprintf(stderr, "  speed.  If the lookup table requires more memory than allowed, the program\n");
    fprintf(stderr, "  exits with an error.\n");
    fprintf(stderr, "    -memory m   Don't use more than m GB memory\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  If -estimate is supplied, processing will stop after a (quick) estimate\n");
    fprintf(stderr, "  of memory needed to load the databases is written to stdout.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Exactly one report type must be specified.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -existence");
    fprintf(stderr, "    Report a tab-delimited line for each sequence showing the number of kmers\n");
    fprintf(stderr, "    in the sequence, in the database, and in both.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    Multiple input -mers may be supplied.  If no output is supplied, output is written\n");
    fprintf(stderr, "    to stdout.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output:  seqName <tab> mersInSeq <tab> mersInDB1 <tab> mersInSeq&DB1 [ <tab> mersInDB2 <tab> mersInSeq&DB2 ... ]\n");
    fprintf(stderr, "      seqName      - name of the sequence\n");
    fprintf(stderr, "      mersInSeq    - number of mers in the sequence\n");
    fprintf(stderr, "      mersInDB     - number of mers in the meryl database\n");
    fprintf(stderr, "      mersInSeq&DB - number of mers in the sequence that are\n");
    fprintf(stderr, "                     also in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -dump\n");
    fprintf(stderr, "    Report a tab-delimited line reporting each kmer in the input sequences, in\n");
    fprintf(stderr, "    order, annotated with the value of the kmer in the input database.  If the kmer\n");
    fprintf(stderr, "    does not exist in the database its value will be reported as zero.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    Only one input may be supplied.  If no output is supplied, output is written\n");
    fprintf(stderr, "    to stdout.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output:  seqName <tab> seqId <tab> seqPos <tab> exists <tab> fwd-mer <tab> fwd-val <tab> rev-mer <tab> rev-val\n");
    fprintf(stderr, "      seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "      seqId      - numeric version of the seqName (0-based)\n");
    fprintf(stderr, "      seqPos     - start position (0-based) of the kmer in the sequence\n");
    fprintf(stderr, "      exists     - 'T' if the kmer exists in the database, 'F' if it does not\n");
    fprintf(stderr, "      fwd-mer    - forward mer sequence\n");
    fprintf(stderr, "      fwd-val    - value of the forward mer in the database\n");
    fprintf(stderr, "      rev-mer    - reverse mer sequence\n");
    fprintf(stderr, "      rev-val    - value of the reverse mer in the database\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -include / -exclude\n");
    fprintf(stderr, "    Extract sequences containing (-include) or not containing (-exclude) kmers in\n");
    fprintf(stderr, "    any input database.  Output sequences are written in the same format as the input\n");
    fprintf(stderr, "    sequences, with the number of kmers found added to the name.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    If two input files are supplied, the corresponding sequences are treated as a pair,\n");
    fprintf(stderr, "    and two output files MUST be supplied.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output:  sequence given format (fasta or fastq) with the number of overlapping kmers appended\n");
    fprintf(stderr, "             if pairs of sequences are given, R1 will be stdout and R2 be named as <output.r2>\n");
    fprintf(stderr, "              <output.r2> will be automatically compressed if ends with .gz, .bz2, or xs\n");
    fprintf(stderr, "      seqName    - name of the sequence this kmer is from\n");
    fprintf(stderr, "      mersInBoth - number of mers in both sequence and in the database\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  G->initialize();
  G->loadLookupTables();
  G->openInputs();
  G->openOutputs();

  switch (G->reportType) {
    case lookupOp::opNone:                               break;
    case lookupOp::opDump:         dumpExistence(G);     break;
    case lookupOp::opExistence:    reportExistence(G);   break;
    case lookupOp::opInclude:      filter(G);            break;
    case lookupOp::opExclude:      filter(G);            break;
    default:                                             break;
  }

  delete G;
  fprintf(stderr, "Bye!\n");

  return(0);
}




