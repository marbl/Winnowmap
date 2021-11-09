
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


void
helpBED(char const *progname) {

  if (progname) {
    fprintf(stderr, "usage: %s [-bed | -bed-runs] \\\n", progname);
    fprintf(stderr, "         -sequence input.fasta \\\n");
    fprintf(stderr, "         -output   output.bed\\\n");
    fprintf(stderr, "         -mers     <input1.meryl> [<input2.meryl>] [...] [-estimate] \\\n");
    fprintf(stderr, "         -labels   <input1name>   [<input2name>]   [...]\n");
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "  -bed:\n");
  fprintf(stderr, "     Generate a BED format file showing the location of kmers in\n");
  fprintf(stderr, "     any input database on each sequence in 'input1.fasta'.\n");
  fprintf(stderr, "     Each kmer is reported in a separate bed record.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -bed-runs:\n");
  fprintf(stderr, "     Generate a BED format file showing the location of kmers in\n");
  fprintf(stderr, "     any input database on each sequence in 'input1.fasta'.\n");
  fprintf(stderr, "     Overlapping kmers are combined into a single bed record.\n");
  fprintf(stderr, "\n");

  if (progname == nullptr)
    return;

  fprintf(stderr, "     If multiple databases are supplied, the output file reports\n");
  fprintf(stderr, "     the location of kmers in any database.  If the -labels option\n");
  fprintf(stderr, "     is supplied, each line will be annotated with the label of the\n");
  fprintf(stderr, "     databse the kmer is found in.  If no -labels are supplied, a\n");
  fprintf(stderr, "     single line will be emitted regardless of how many databases the\n");
  fprintf(stderr, "     kmer is found in.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     For example, with two input databases (-mers A.meryl B.meryl) and\n");
  fprintf(stderr, "     labels supplied (-labels A B), with the first two kmers in the input\n");
  fprintf(stderr, "     sequence found in both databases, the output will be:\n");
  fprintf(stderr, "       sequence1 <tab> 0 <tab> 21 <tab> A\n");
  fprintf(stderr, "       sequence1 <tab> 0 <tab> 21 <tab> B\n");
  fprintf(stderr, "       sequence1 <tab> 1 <tab> 22 <tab> A\n");
  fprintf(stderr, "       sequence1 <tab> 1 <tab> 22 <tab> B\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Without the -labels option, the output will be:\n");
  fprintf(stderr, "       sequence1 <tab> 0 <tab> 21\n");
  fprintf(stderr, "       sequence1 <tab> 1 <tab> 22\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     If -bed-runs is used, the output will bed:\n");
  fprintf(stderr, "       sequence1 <tab> 0 <tab> 22 <tab> A\n");
  fprintf(stderr, "       sequence1 <tab> 0 <tab> 22 <tab> B\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Output lines are written in the order sequences appear in the input\n");
  fprintf(stderr, "     file, and are in increasing position within the sequence itself.\n");
  fprintf(stderr, "\n");
}

void
helpWIGcount(char const *progname) {

  if (progname) {
    fprintf(stderr, "usage: %s -wig-count \\\n", progname);
    fprintf(stderr, "         -sequence input.fasta \\\n");
    fprintf(stderr, "         -output   output.wig \\\n");
    fprintf(stderr, "         -mers     <input1.meryl> [<input2.meryl>] [...] [-estimate] \\\n");
  }

  fprintf(stderr, "  -wig-count:\n");
  fprintf(stderr, "     Generate a WIGGLE format file showing the multiplicity of the\n");
  fprintf(stderr, "     kmer starting at each position in the sequence, if it exists in\n");
  fprintf(stderr, "     an input kmer database.\n");
  fprintf(stderr, "\n");

  if (progname == nullptr)
    return;

  fprintf(stderr, "     If multiple databases are supplied, the reported multiplicity\n");
  fprintf(stderr, "     is the sum of multiplicities in all databases.  The -labels\n");
  fprintf(stderr, "     option is not used.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Exactly one input -sequence must be supplied.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     If no -output path is supplied, output is written to stdout.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     The output file has format:\n");
  fprintf(stderr, "         variableStep chrom=<sequence_name>\n");
  fprintf(stderr, "         <position> <tab> <sum_of_multiplicities>\n");
  fprintf(stderr, "         <position> <tab> <sum_of_multiplicities>\n");
  fprintf(stderr, "\n");

}

void
helpWIGdepth(char const *progname) {

  if (progname) {
    fprintf(stderr, "usage: %s -wig-depth \\\n", progname);
    fprintf(stderr, "         -sequence input.fasta \\\n");
    fprintf(stderr, "         -output   output.wig \\\n");
    fprintf(stderr, "         -mers     <input1.meryl> [<input2.meryl>] [...] [-estimate] \\\n");
  }

  fprintf(stderr, "  -wig-depth:\n");
  fprintf(stderr, "     Generate a WIGGLE format file showing the number of kmers in\n");
  fprintf(stderr, "     any input database that cover each position in the sequence.\n");
  fprintf(stderr, "\n");

  if (progname == nullptr)
    return;

  fprintf(stderr, "     If multiple databases are supplied, the depth does not change\n");
  fprintf(stderr, "     when the kmer is present in more than one database.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Exactly one input -sequence must be supplied.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     If no -output path is supplied, output is written to stdout.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     The output file has format:\n");
  fprintf(stderr, "         variableStep chrom=<sequence_name>\n");
  fprintf(stderr, "         <position> <tab> <kmer_depth>\n");
  fprintf(stderr, "         <position> <tab> <kmer_depth>\n");
  fprintf(stderr, "\n");
}

void
helpExistence(char const *progname) {

  if (progname) {
    fprintf(stderr, "usage: %s -wig-depth \\\n", progname);
    fprintf(stderr, "         -sequence input.fasta \\\n");
    fprintf(stderr, "         -output   output.wig \\\n");
    fprintf(stderr, "         -mers     <input1.meryl> [<input2.meryl>] [...] [-estimate] \\\n");
  }

  fprintf(stderr, "  -existence:\n");
  fprintf(stderr, "     Generate a tab-delimited line for each input sequence with the\n");
  fprintf(stderr, "     number of kmers in the sequence, in the database and common to both.\n");
  fprintf(stderr, "\n");

  if (progname == nullptr)
    return;

  fprintf(stderr, "     The number of kmers in common is counted individually for each\n");
  fprintf(stderr, "     database, but still reported on a single output line.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Exactly one input -sequence must be supplied.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     The -labels option is not used.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     The output file has format:\n");
  fprintf(stderr, "         <sequence_name> <tab> <kmers_in_sequence> <tab> <kmers_in_db> <tab> <kmers_shared> [...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     With one input database:\n");
  fprintf(stderr, "         sequence1 <tab> 8415 <tab> 12856825 <tab> 8145\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     With two input databases:\n");
  fprintf(stderr, "         sequence1 <tab> 8415 <tab> 12856825 <tab> 8145 <tab> 575757256 <tab> 8354\n");
  fprintf(stderr, "\n");
}

void
helpIncludeExclude(char const *progname) {

  if (progname) {
    fprintf(stderr, "usage: %s [-include | -exclude] \\\n", progname);
    fprintf(stderr, "         -sequence <input1.fasta> [<input2.fasta>] \\\n");
    fprintf(stderr, "         -output   <output1>      [<output2>] \\\n");
    fprintf(stderr, "         -mers     <input1.meryl> [-estimate] \\\n");
    fprintf(stderr, "         -10x\n");
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "  -include:\n");
  fprintf(stderr, "  -exclude:\n");
  fprintf(stderr, "     Copy sequences from 'input1.fasta' (and 'input2.fasta') to the\n");
  fprintf(stderr, "     corresponding output file if the sequence has at least one kmer\n");
  fprintf(stderr, "     present (include) or no kmers present (exclude) in 'input1.meryl'.\n");
  fprintf(stderr, "\n");

  if (progname == nullptr)
    return;

  fprintf(stderr, "  -10x:\n");
  fprintf(stderr, "     When -10x is supplied, the first 23 bp of every sequence in input1.fasta\n");
  fprintf(stderr, "     will be ignored while looking up for kmer existence.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     Exactly one input database must be supplied.  The -labels option is\n");
  fprintf(stderr, "     not used.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     When one input sequence is supplied, each sequence is copied to the\n");
  fprintf(stderr, "     output file if kmers exist (-include) or do not exist (-exclude) in\n");
  fprintf(stderr, "     the input databse.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "     When two input sequence is supplied, the pair of sequences are copied\n");
  fprintf(stderr, "     the the output files if kmers from either sequence exist (-include) or\n");
  fprintf(stderr, "     do not exist (-exclude) in the input databse.\n");
  fprintf(stderr, "\n");
}



void
help(char const *progname) {
  fprintf(stderr, "usage: %s <report-type> \\\n", progname);
  fprintf(stderr, "         -sequence <input1.fasta> [<input2.fasta>] \\\n");
  fprintf(stderr, "         -output   <output1>      [<output2>] \\\n");
  fprintf(stderr, "         -mers     <input1.meryl> [<input2.meryl>] [...] [-estimate] \\\n");
  fprintf(stderr, "         -labels   <input1name>   [<input2name>]   [...]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Compare kmers in input sequences against kmers in input meryl databases.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Input sequences (-sequence) can be FASTA or FASTQ, uncompressed, or\n");
  fprintf(stderr, "  compressed with gzip, xz, or bzip2.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  To compute and report only estimated memory usage, add option '-estimate'.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  Report types:\n");
  fprintf(stderr, "    Run `%s <report-type> -help` for details on each method.\n", progname);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  helpBED();
  helpWIGcount();
  helpWIGdepth();
  helpExistence();
  helpIncludeExclude();
}
