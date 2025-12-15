
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

void
merylOperation::reportHistogram(void) {

  if (_inputs.size() > 1)
    fprintf(stderr, "ERROR: told to dump a histogram for more than one input!\n"), exit(1);

  if (_inputs[0]->_operation)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  if (_inputs[0]->_sequence)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  //  Tell the stream to load and return the histogram.

  merylHistogram *stats = _inputs[0]->_stream->stats();

  //  Now just dump it.

  for (uint32 ii=0; ii<stats->histogramLength(); ii++)
    fprintf(stdout, F_U64 "\t" F_U64 "\n",
            stats->histogramValue(ii),
            stats->histogramOccurrences(ii));
}



void
merylOperation::reportStatistics(void) {

  if (_inputs.size() > 1)
    fprintf(stderr, "ERROR: told to dump a histogram for more than one input!\n"), exit(1);

  if (_inputs[0]->_operation)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  if (_inputs[0]->_sequence)
    fprintf(stderr, "ERROR: told to dump a histogram from input '%s'!\n", _inputs[0]->_name), exit(1);

  //  Tell the stream to load and return the histogram.

  merylHistogram *stats = _inputs[0]->_stream->stats();

  //  Now just dump it.

  uint64  nUniverse = buildLowBitMask<uint64>(kmer::merSize() * 2) + 1;
  uint64  sDistinct = 0;
  uint64  sTotal    = 0;

  fprintf(stdout, "Number of %u-mers that are:\n", kmer::merSize());
  fprintf(stdout, "  unique   %20" F_U64P "  (exactly one instance of the kmer is in the input)\n", stats->numUnique());
  fprintf(stdout, "  distinct %20" F_U64P "  (non-redundant kmer sequences in the input)\n", stats->numDistinct());
  fprintf(stdout, "  present  %20" F_U64P "  (...)\n", stats->numTotal());
  fprintf(stdout, "  missing  %20" F_U64P "  (non-redundant kmer sequences not in the input)\n", nUniverse - stats->numDistinct());
  fprintf(stdout, "\n");
  fprintf(stdout, "             number of   cumulative   cumulative     presence\n");
  fprintf(stdout, "              distinct     fraction     fraction   in dataset\n");
  fprintf(stdout, "frequency        kmers     distinct        total       (1e-6)\n");
  fprintf(stdout, "--------- ------------ ------------ ------------ ------------\n");

  for (uint32 ii=0; ii<stats->histogramLength(); ii++) {
    uint64  value = stats->histogramValue(ii);
    uint64  occur = stats->histogramOccurrences(ii);

    sDistinct  += occur;
    sTotal     += occur * value;

    fprintf(stdout, "%9" F_U64P " %12" F_U64P " %12.4f %12.4f %12.6f\n",
            value,
            occur,
            (double)sDistinct / stats->numDistinct(),
            (double)sTotal    / stats->numTotal(),
            (double)value     / stats->numTotal() * 1000000.0);
  }

  assert(sDistinct == stats->numDistinct());
  assert(sTotal    == stats->numTotal());
}






void
merylOperation::reportPloidy(void) {

  if (_inputs.size() > 1)
    fprintf(stderr, "ERROR: told to dump ploidy for more than one input!\n"), exit(1);

  if (_inputs[0]->_operation)
    fprintf(stderr, "ERROR: told to dump ploidy from input '%s'!\n", _inputs[0]->_name), exit(1);

  if (_inputs[0]->_sequence)
    fprintf(stderr, "ERROR: told to dump ploidy from input '%s'!\n", _inputs[0]->_name), exit(1);

  merylHistogram  *stats        = nullptr;
  merylHistogram  *statsPrivate = nullptr;
  bool             debugPloidy  = (_verbosity >= sayDetails);

  if      (_inputs[0]->isFromDatabase() == true) {
    stats = _inputs[0]->_stream->stats();
    stats->computePloidyPeaks(debugPloidy);
  }
  else if (_inputs[0]->isFromFile() == true) {
    stats = statsPrivate = new merylHistogram();

    stats->load(_inputs[0]->_filename);
    stats->computePloidyPeaks(debugPloidy);
  }
  else {
    fprintf(stderr, "'ploidy' can't run from input type '%s'.\n", _inputs[0]->inputType()), exit(1);
  }

  double no = stats->getNoiseTrough();

  double c1 = stats->getCoverage(1);
  double c2 = stats->getCoverage(2);
  double c3 = stats->getCoverage(3);
  double c4 = stats->getCoverage(4);

  double p1 = stats->getDepth(1);
  double p2 = stats->getDepth(2);
  double p3 = stats->getDepth(3);
  double p4 = stats->getDepth(4);

  fprintf(stderr, "\n");
  fprintf(stderr, "Noise/genomic trough: %6.3f\n", no);
  fprintf(stderr, "%4.2fx coverage peak:   %6.3f\n", c1, p1);
  fprintf(stderr, "%4.2fx coverage peak:   %6.3f\n", c2, p2);
  fprintf(stderr, "%4.2fx coverage peak:   %6.3f\n", c3, p3);
  fprintf(stderr, "%4.2fx coverage peak:   %6.3f\n", c4, p4);

  if (isatty(fileno(stdout)) == 0)
    fprintf(stdout, "noise-trough\t%.3f\tploidy-peaks\t%.3f\t%.3f\t%.3f\t%.3f\n", no, p1, p2, p3, p4);

  delete statsPrivate;
}
