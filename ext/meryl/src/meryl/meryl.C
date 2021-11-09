
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

#include "runtime.H"
#include "strings.H"
#include "system.H"



int
main(int argc, char **argv) {
  merylCommandBuilder  *B = new merylCommandBuilder;

  argc = AS_configure(argc, argv);

  std::vector<char *>  err;
  for (int32 arg=1; arg < argc; arg++) {

    //
    //  Scan for debug options and requests for help.
    //

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

    if ((strcmp(argv[arg], "-h")   == 0) ||
        (strcmp(argv[arg], "help") == 0)) {
      err.push_back(NULL);
      continue;
    }

    //
    //
    //

    B->initialize(argv[arg]);

    if (B->processOptions() == true)       //  Process and add options to the current command.
      continue;

    if (B->processOperation() == true)     //  Detect a new operation.
      continue;

    if (B->isOutput() == true)             //  Handle 'output' and 'print' flags, and their
      continue;                           //  (possibly optional) output path.
    if (B->isPrinter() == true)
      continue;

    if (B->isMerylInput() == true)         //  Last, try to make an input.  These must come
      continue;                           //  after 'print' so we can do 'print some.meryl'
    if (B->isCanuInput(err) == true)       //  detect that as an input.
      continue;
    if (B->isSequenceInput() == true)
      continue;

    //
    //  Otherwise, we have no idea what this word is.
    //

    char *s = new char [1024];
    snprintf(s, 1024, "Can't interpret '%s': not a meryl command, option, or recognized input file.", argv[arg]);
    err.push_back(s);
  }

  B->finalize();   //  Finalize the command tree.

  //  If any errors, fail.

  if ((argc == 1) ||                //  No commands
      (B->numRoots() == 0) ||   //  No actions
      (err.size() > 0)) {           //  Errors
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  A meryl command line is formed as a series of commands and files, possibly\n");
    fprintf(stderr, "  grouped using square brackets.  Each command operates on the file(s) that\n");
    fprintf(stderr, "  are listed after it.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  COMMANDS:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    statistics           display total, unique, distnict, present number of the kmers on the screen.  accepts exactly one input.\n");
    fprintf(stderr, "    histogram            display kmer frequency on the screen as 'frequency<tab>count'.  accepts exactly one input.\n");
    fprintf(stderr, "    print                display kmers on the screen as 'kmer<tab>count'.  accepts exactly one input.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    count                Count the occurrences of canonical kmers in the input.  must have 'output' specified.\n");
    fprintf(stderr, "    count-forward        Count the occurrences of forward kmers in the input.  must have 'output' specified.\n");
    fprintf(stderr, "    count-reverse        Count the occurrences of reverse kmers in the input.  must have 'output' specified.\n");
    fprintf(stderr, "      k=<K>              create mers of size K bases (mandatory).\n");
    fprintf(stderr, "      n=<N>              expect N mers in the input (optional; for precise memory sizing).\n");
    fprintf(stderr, "      memory=M           use no more than (about) M GB memory.\n");
    fprintf(stderr, "      threads=T          use no more than T threads.\n");
    fprintf(stderr, "      compress           compress homopolymer runs to a single letter.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    less-than N          return kmers that occur fewer than N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "    greater-than N       return kmers that occur more than N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "    equal-to N           return kmers that occur exactly N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "    not-equal-to N       return kmers that do not occur exactly N times in the input.  accepts exactly one input.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    increase X           add X to the count of each kmer.\n");
    fprintf(stderr, "    decrease X           subtract X from the count of each kmer.\n");
    fprintf(stderr, "    multiply X           multiply the count of each kmer by X.\n");
    fprintf(stderr, "    divide X             divide the count of each kmer by X.\n");
    fprintf(stderr, "    divide-round X       divide the count of each kmer by X and round results. count < X will become 1.\n");
    fprintf(stderr, "    modulo X             set the count of each kmer to the remainder of the count divided by X.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    union                return kmers that occur in any input, set the count to the number of inputs with this kmer.\n");
    fprintf(stderr, "    union-min            return kmers that occur in any input, set the count to the minimum count\n");
    fprintf(stderr, "    union-max            return kmers that occur in any input, set the count to the maximum count\n");
    fprintf(stderr, "    union-sum            return kmers that occur in any input, set the count to the sum of the counts\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    intersect            return kmers that occur in all inputs, set the count to the count in the first input.\n");
    fprintf(stderr, "    intersect-min        return kmers that occur in all inputs, set the count to the minimum count.\n");
    fprintf(stderr, "    intersect-max        return kmers that occur in all inputs, set the count to the maximum count.\n");
    fprintf(stderr, "    intersect-sum        return kmers that occur in all inputs, set the count to the sum of the counts.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    subtract             return kmers that occur in the first input, subtracting counts from the other inputs\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    difference           return kmers that occur in the first input, but none of the other inputs\n");
    fprintf(stderr, "    symmetric-difference return kmers that occur in exactly one input\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  MODIFIERS:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    output O             write kmers generated by the present command to an output  meryl database O\n");
    fprintf(stderr, "                         mandatory for count operations.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  EXAMPLES:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Example:  Report 22-mers present in at least one of input1.fasta and input2.fasta.\n");
    fprintf(stderr, "            Kmers from each input are saved in meryl databases 'input1' and 'input2',\n");
    fprintf(stderr, "            but the kmers in the union are only reported to the screen.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "            meryl print \\\n");
    fprintf(stderr, "                    union \\\n");
    fprintf(stderr, "                      [count k=22 input1.fasta output input1] \\\n");
    fprintf(stderr, "                      [count k=22 input2.fasta output input2]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Example:  Find the highest count of each kmer present in both files, save the kmers to\n");
    fprintf(stderr, "            database 'maxCount'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "            meryl intersect-max input1 input2 output maxCount\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Example:  Find unique kmers common to both files.  Brackets are necessary\n");
    fprintf(stderr, "            on the first 'equal-to' command to prevent the second 'equal-to' from\n");
    fprintf(stderr, "            being used as an input to the first 'equal-to'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "            meryl intersect [equal-to 1 input1] equal-to 1 input2\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    exit(1);
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "Found %u command tree%s.\n", B->numRoots(), (B->numRoots() == 1) ? "" : "s");

  //  opHistogram is limited to showing only histograms already stored in a database.
  //  opHistogram cannot take input from anything but a database.
  //  opHistogram does not generate kmer outputs.
  //  So, if the top op is histogram, all we can do is load the histogram and dump it.
  //
  //  Eventully, maybe, opHistogram will pass through mers (but then we need to figure out
  //  where to report the histogram).
  //
  //  Eventually, maybe, opHistogram will allow input from a kmer stream.

  if (B->getOperation(0)->getOperation() == opHistogram) {
    B->getOperation(0)->initialize();
    B->getOperation(0)->reportHistogram();
    exit(0);
  }

  if (B->getOperation(0)->getOperation() == opStatistics) {
    B->getOperation(0)->initialize();
    B->getOperation(0)->reportStatistics();
    exit(0);
  }

  //  Counting operations are a big headache.  They don't fit into the
  //  tree nicely:
  //   - they do their own threading, so only one thread can start the operation
  //   - when done, they transform themselves into a pass-through operation that
  //     simply reads the (just counted) input and passes kmers through.
  //
  //  So, we special case them here.  Process in order, counting, writing the
  //  output, and converting to a pass-through operation.

  for (uint32 oo=0; oo<B->numOperations(); oo++) {
    merylOperation *op = B->getOperation(oo);

    if (op->isCounting() == true) {
      op->initialize();
      op->doCounting();
    }
  }

  //  Initialize all the root nodes.  This initializes constants and,
  //  importantly, opens outputs.

  for (uint32 rr=0; rr<B->numRoots(); rr++)
    B->getRoot(rr)->initialize(true);

  //  Initialize nodes for all the threads.  All the root nodes need to be
  //  initialized before we spawn, so we get thresholds set correctly.

  B->spawnThreads();

  //  Process each file, in parallel.  Just keep getting the next mer and let
  //  each op do their work.

  for (uint32 rr=0; rr<B->numRoots(); rr++) {
    merylOperation *root = B->getRoot(rr);

    if (root->getOperation() == opPassThrough)
      continue;

    fprintf(stderr, "\n");
    fprintf(stderr, "PROCESSING TREE #%u using %u thread%s.\n", rr+1, getMaxThreadsAllowed(), getMaxThreadsAllowed() == 1 ? "" : "s");
    B->printTree(root, 2);

#pragma omp parallel for schedule(dynamic, 1)
    for (uint32 ff=0; ff<64; ff++) {
      merylOperation *op = B->getRoot(rr, ff);

      if (op->initialize() == true)
        while (op->nextMer() == true)
          ;

      op->finalize();
    }
  }

  //  Now that everything is done, delete!

  fprintf(stderr, "\n");
  fprintf(stderr, "Cleaning up.\n");

  delete B;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  return(0);
}
