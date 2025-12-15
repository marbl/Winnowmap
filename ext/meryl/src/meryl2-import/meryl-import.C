
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

#include "kmers.H"
#include "sequence.H"
#include "strings.H"
#include "bits.H"

#include "merylCountArray.H"

using namespace merylutil;
using namespace merylutil::kmers::v2;

int
main(int argc, char **argv) {
  char   *inputName    = nullptr;
  char   *outputDBname = nullptr;
  uint32  kLen         = 0;

  uint32  valueWidth   = 0;
  uint32  labelWidth   = 0;

  bool    doMultiSet   = false;

  bool    useC         = true;
  bool    useF         = false;

  uint32  threads      = getMaxThreadsAllowed();

  argc = AS_configure(argc, argv);

  std::vector<char const *>  err;
  int                        arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-kmers") == 0) {
      inputName = argv[++arg];

    } else if (strcmp(argv[arg], "-output") == 0) {
      outputDBname = argv[++arg];

    } else if (strcmp(argv[arg], "-k") == 0) {
      kLen = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-maxvalue") == 0) {
      valueWidth = countNumberOfBits64(strtouint64(argv[++arg]));

    } else if (strcmp(argv[arg], "-valuewidth") == 0) {
      valueWidth =strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-labelwidth") == 0) {
      labelWidth = strtouint32(argv[++arg]);
      kmerTiny::setLabelSize(labelWidth);

    } else if (strcmp(argv[arg], "-multiset") == 0) {
      doMultiSet = true;

    } else if (strcmp(argv[arg], "-forward") == 0) {
      useC = false;
      useF = true;

    } else if (strcmp(argv[arg], "-reverse") == 0) {
      useC = false;
      useF = false;

    } else if (strcmp(argv[arg], "-threads") == 0) {
      threads = strtouint32(argv[++arg]);

    } else {
      char *s = new char [1024];
      snprintf(s, 1024, "Unknown option '%s'.\n", argv[arg]);
      err.push_back(s);
    }

    arg++;
  }

  if (inputName == nullptr)
    err.push_back("No input kmer file (-kmers) supplied.\n");
  if (outputDBname == nullptr)
    err.push_back("No output database name (-output) supplied.\n");
  if (kLen == 0)
    err.push_back("No kmer size (-k) supplied.\n");

  if (err.size() > 0) {
    fprintf(stderr, "usage: %s [...] -k <kmer-size> -kmers <input-kmers> -output <db.meryl>\n", argv[0]);
    fprintf(stderr, "  Loads the kmers and values listed in <input-kmers> into a meryl kmer database.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "INPUTS and OUTPUTS\n");
    fprintf(stderr, "  -kmers <input-kmers>  A file consisting of kmers and values, one per line, separated\n");
    fprintf(stderr, "                        by white space ('AGTTGCC 4').  Order of kmers is not important.\n");
    fprintf(stderr, "                        Duplicate kmers will be handled according to the -multiset\n");
    fprintf(stderr, "                        option.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "                        A persistent value can be specified as '#<value>' (e.g., '#3')\n");
    fprintf(stderr, "                        All kmers with no value after this line will use this value.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -k <size>             The size of a kmer, in bases.  Setting this larger than the\n");
    fprintf(stderr, "                        kmers in the input will probably lead to a crash.  Setting it\n");
    fprintf(stderr, "                        smaller will result in only the left-most bases being used.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -output <db.meryl>    Create (or overwrite) meryl database 'database.meryl'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONS\n");
    fprintf(stderr, "  -multiset             Write duplicate kmers in the input to the database as individual\n");
    fprintf(stderr, "                        entries.  A kmer AGTTGCC in the input twice with values 4 and 7\n");
    fprintf(stderr, "                        will be listed in the output database twice, once with value 4,\n");
    fprintf(stderr, "                        and once with value 7.  Without this option, the values will be\n");
    fprintf(stderr, "                        summed: AGTTGCC will be listed once with value 11.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -maxvalue <value>     An optional memory and time optimization, useful if your values\n");
    fprintf(stderr, "                        are randomly distributed and below some known maximum value.\n");
    fprintf(stderr, "                        For data whose values are the counts from actual data, it is\n");
    fprintf(stderr, "                        probably best to not use this option.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -valuewidth <vw>      Explicitly set the width (in bits) of values/labels to store\n");
    fprintf(stderr, "  -labelwidth <lw>      with each kmer.  If vw is zero, values are stored using the default\n");
    fprintf(stderr, "                        encoding (recommended); if lw is zero, labels are NOT stored.\n");
    fprintf(stderr, "                        The default value for both is zero.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -forward              By default, the canonical kmer is loaded into the database.  These\n");
    fprintf(stderr, "  -reverse              options force either the forward or reverse-complement kmer to be\n");
    fprintf(stderr, "                        loaded instead.  These options are mutually exclusive.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -threads <t>          Use <t> compute threads when sorting and writing data.\n");
    fprintf(stderr, "                        (most of the time is NOT spent here though)\n");
    fprintf(stderr, "\n");

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii])
        fputs(err[ii], stderr);

    exit(1);
  }

  omp_set_num_threads(threads);

  //  Figure out the kmer size.  We need this to set up encoding parameters.

  kmerTiny::setSize(kLen);

  //  Decide on some parameters.

  uint32  wPrefix   = 6;
  uint32  nPrefix   = 1 << wPrefix;

  uint32  wData     = 2 * kmerTiny::merSize() - wPrefix;
  uint64  wDataMask = buildLowBitMask<uint64>(wData);

  //  Open the input kmer file, allocate space for reading kmer lines.

  FILE         *K    = merylutil::openInputFile(inputName);
  uint32        Llen = 0;
  uint32        Lmax = 1023;
  char         *L    = new char [Lmax + 1];

  splitToWords  W;

  kmerTiny      kmerF;
  kmerTiny      kmerR;

  kmvalu        val, pval;
  kmlabl        lab, plab;

  uint64        nKmers = 0;
  uint64        nErrrs = 0;

  //  Allocate a bunch of counting arrays.  The naming here follows merylOp-count.C.

  fprintf(stderr, "Allocating %u count arrays.\n", nPrefix);

  merylCountArray  *data  = new merylCountArray [nPrefix];

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 pp=0; pp<nPrefix; pp++) {
    data[pp].initialize(pp, wData);
    data[pp].initializeValues(valueWidth, labelWidth);

    data[pp].enableMultiSet(doMultiSet);
  }

  //  Read each kmer/value/label, stuff into a merylCountArray, writing when the array is full.


  while (merylutil::readLine(L, Llen, Lmax, K) == true) {
    splitToWords   W(L);
    char const    *word = W[0];

    //  If an empty line, continue.

    if (word == nullptr)
      continue;

    //  Got a line.  Set any persistent value or label, then continue.

    if (strncmp(word, "value=", 6) == 0) {
      pval = decodeInteger(W[0], 6, 0, pval, err);
      continue;
    }

    if (strncmp(word, "label=", 6) == 0) {
      plab = decodeInteger(W[0], 6, 0, plab, err);
      continue;
    }

    //  Make a kmer.  Decode the value/label and encode the kmer.

    for (uint32 ii=0; word[ii]; ii++)
      kmerF.addR(word[ii]);

    if (W.numWords() > 1)
      kmerF._val = decodeInteger(W[1], 0, 0, kmerF._val, err);
    else
      kmerF._val = pval;

    if (W.numWords() > 2)
      kmerF._lab = decodeInteger(W[2], 0, 0, kmerF._lab, err);
    else
      kmerF._lab = plab;

    //  Add the kmer - canonical, forward or reverse as per -forward/-reverse switches.

    if ((useC == true) || (useF == false)) {
      kmerR = kmerF;
      kmerR.reverseComplement();
    }

    if (useC == true)
      useF = (kmerF < kmerR) ? true : false;

    kmdata  pp = (useF == true) ? ((kmdata)kmerF >> wData)     : ((kmdata)kmerR >> wData);
    kmdata  mm = (useF == true) ? ((kmdata)kmerF  & wDataMask) : ((kmdata)kmerR  & wDataMask);

    assert(pp < nPrefix);

    uint64 ak = data[pp].add(mm);
    uint64 av = data[pp].addValue(kmerF._val);
    uint64 al = data[pp].addLabel(kmerF._lab);

    if ((kmerF._val > 0) && (av == 0)) {
      char *a = new char [1024];
      snprintf(a, 1024, "failed to correctly set value %u\n", kmerF._val);
      err.push_back(a);
    }

    if ((kmerF._lab > 0) && (al == 0)) {
      char *a = new char [1024];
      snprintf(a, 1024, "failed to correctly set label %lu\n", kmerF._lab);
      err.push_back(a);
    }

    if (err.size() == 0) {
      nKmers++;
    }
    else {
      nErrrs++;

      fprintf(stderr, "ERROR: kmer '%s' failed:\n", word);
      for (uint32 ii=0; ii<err.size(); ii++)
        fprintf(stderr, "ERROR:   %s\n", err[ii]);
    }
  }

  //  All data loaded, cleanup.

  fprintf(stderr, "Added %lu kmers; incorrectly added %lu kmers.\n", nKmers, nErrrs);

  //  And dump to the output.

  merylFileWriter   *output = new merylFileWriter(outputDBname);

  output->initialize(wPrefix);

  merylBlockWriter  *writer = output->getBlockWriter();

#pragma omp parallel for schedule(dynamic, 1)
  for (uint32 ff=0; ff<output->numberOfFiles(); ff++) {
    for (uint64 pp=output->firstPrefixInFile(ff); pp <= output->lastPrefixInFile(ff); pp++) {
      data[pp].countKmers();                   //  Convert the list of kmers into a list of (kmer, count).
      data[pp].dumpCountedKmers(writer, 0);    //  Write that list to disk.
      data[pp].removeCountedKmers();           //  And remove the in-core data.
    }
  }

  writer->finish();

  delete    writer;
  delete    output;
  delete [] data;

  fprintf(stderr, "\n");
  fprintf(stderr, "Bye.\n");

  exit(0);
}
