.. _quick-start:

Quick Start
===========

To give a gentle introduction to meryl, we'll use three random Escherichia
coli assemblies pulled from GenBank.  They're not `quite` random; they're
flagged as ``complete`` and have the highest number of ``scaffolds`` - which,
we hope, means they have a complete nuclear genome and a few handfuls of
plasmids.

We'll use
`EC931 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA608094>`_,
`SCEC020022 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA418674>`_ and
`MSB1_4I-sc-2280412 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA646837>`_,
but give them shorter names for convenience.

.. code-block:: none

  mkdir data

  curl -LRO ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/522/225/GCA_014522225.1_ASM1452222v1/GCA_014522225.1_ASM1452222v1_genomic.fna.gz
  curl -LRO ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/165/095/GCA_002165095.2_ASM216509v2/GCA_002165095.2_ASM216509v2_genomic.fna.gz
  curl -LRO ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/071/835/GCA_905071835.1_MSB1_4I/GCA_905071835.1_MSB1_4I_genomic.fna.gz

  mv -i GCA_014522225.1_ASM1452222v1_genomic.fna.gz data/ec.fna.gz
  mv -i GCA_002165095.2_ASM216509v2_genomic.fna.gz  data/sc.fna.gz
  mv -i GCA_905071835.1_MSB1_4I_genomic.fna.gz      data/ms.fna.gz

First, lets count the kmers in each and save the results in meryl databases.
 
.. code-block:: shell

  % meryl count k=42 data/ec.fna.gz output ec.meryl

This command does just as it reads (left to right): count 42-mers in a FASTA input and output the results to a meryl database.

We can show the k-mers in the database with:

.. code-block:: shell

  % meryl threads=1 print ec.meryl | head

  AAAAAAAAAAACGGATCTGCATCACAACGAGAGTGATCCCAC      1
  AAAAAAAAAACTTCCACCAATATGATGGGTGGCGTACAGTAA      1
  AAAAAAAAAACGGATCTGCATCACAACGAGAGTGATCCCACC      1
  AAAAAAAAAATAGACAAAAAATACTTTATCAAAACATACATA      1
  AAAAAAAAACCATCCAAATCTGGATGGCTTTTCATAATTCTG      1
  AAAAAAAAACCCGATTTTATCAGGCGTCAACCTCTGAATTGT      1
  AAAAAAAAACCCGCTGATTAAGCGGGTTTTGAATTCTTGCTG      1
  AAAAAAAAACCTGAAAAAACGGCCTGACGTGAATCAAGCAAT      1
  AAAAAAAAACCTGCCATCGCTGGCAGGTTTTTTATGACTAAA      1
  AAAAAAAAACCGACCAAGGGTCGGGGCAAGAATCAGAGTCTG      1

(Detail: Option 'threads=1' is supplied to prevent meryl from processing its
64 data chunks in parallel; for the print operation this results in the
output kmers being unsorted.)

Let's now count the other two databases and combine the kmers from all three
genomes into one database, all in one command.

.. code-block:: shell

  % meryl \
      union-sum \
      output union.meryl \
      [count data/sc.fna.gz output sc.meryl] \
      [count output ms.meryl data/ms.fna.gz] \
      ec.meryl

There's a lot going on here.  There are three 'actions' and three 'output'
files.  Each action starts a new 'nesting-level'.  A nesting-level has
exactly one action, at most one output, and any number of inputs.

The union-sum action is at nesting-level 1, and writes output to
'union.meryl'.  Continuing to read the command line left to right, we find a
new action, the counting of kmers in sc.fna.gz.  The output of this action is
fed as input to the union-sum action.  Since this is a new action, it makes a
new nesting-level (2).  The square brackets serve to cordon off a nesting
level.  When the right square bracket is enountered, the nesting-level is
decremented back to 1.  Then another count action is encountered, which also
feeds it's output to the union-sum action.  Finally, a meryl database is
enountered.  Since we're back to nesting-level 1, this, too, becomes input to
the union-sum command.

The end result of this is to count kmers in sc.fna.gz and ms.fna.gz, writing
the kmers from each to s.meryl and ms.meryl, respectively.  Finally, the kmers
from all three genomes are combined into one database, summing the counts for
kmers that appear in multiple databases.

Notice also that the 'k=42' option is not present.  Meryl will determine the
kmer size from the 'ec.meryl' input database and use that for all operations.
If a kmer size is supplied, it must match the database, and all databases must
be of the same size.

Before meryl starts processing, it will output a summary of the processing it
will perform.  Note, however, that counting operations occur before this
summary is reported.

.. code-block:: none

  PROCESSING TREE #1 using 16 threads.
    opUnionSum
      opPassThrough
        sc.meryl
      opPassThrough
        mc.meryl
      ec.meryl
      output to union.meryl

Detail: The 'count' actions are converted to 'pass-through' actions once the
kmers are counted.

A meryl database also stores the histogram of kmer values.  This can be displayed:

.. code-block:: shell

  % meryl histogram ec.meryl

  1       4911809
  2       37336
  3       7632
  4       1217
  5       2705
  6       2232
  7       4544
  8       384
  9       862
  11      3
  12      4
  15      967
  16      230
  18      1
  19      1
  21      81
  22      3
  27      39
  29      4
  30      21
  31      19
  32      5
  37      39
  48      3
  49      42
  50      38
  51      15
  52      493
  53      13
  
Which hints there is a 52 copy repeat of around 500 bases in Escherichia coli
EC931.  Histograms from the other two genomes show either no high copy repeat
(Escherichia coli SCEC020022, 'sc.meryl') or a potential 64 copy repeat
(Escherichia coli MSB1_4I-sc-2280412, 'ms.meryl').  Lets now extract those
kmers and see where they are on the genomes.

.. code-block:: shell

  % meryl print at-least 48 ec.meryl output ec-repeats.meryl > ec-repeats.dump

This command does two things: it creates a new meryl database of just the
repeat kmers, and creates a text file of those kmers.

The `meryl-lookup` tool compares FAST/FASTQ sequences against a meryl
database (or several databases).  We'll use it to generate a bed file
of the bases covered by kmers in a database.

.. code-block:: shell

  % meryl-lookup -sequence data/ec.fna -mers ec-repeats.meryl -bed-runs > ec-repeats.bed

From this we see that there is not a single 52-copy repeat, but several shorter repeats.
There are several instances of a 447 base repeat with single base differences:

.. code-block:: none

  CP049118.1      3782667 3782794
  CP049118.1      3782794 3783114

  CP049118.1      3762937 3763257
  CP049118.1      3763257 3763384

Not really part of meryl, the high-count kmers can be passed to a greedy
assembler with nice results (the greedy assembler is included in the meryl
source code, but isn't installed in the binary directory).

.. code-block:: shell

  % perl $MERYL/scripts/greedy-assemble-kmers.pl < ec-repeats.dump

  >1
  AGCCTGTCATACGCGTAAAACAGCCAGCGCTGGCGCGATTTAGCCCCGACATAGCCCCACTGTTCGTCCATTTCCGCGCAGACGATGACGTCACTGCCCG
  GCTGTATGCGCGAGGTTACCGACTGCGGCCTGAGTTTTTTAAGTGACGTAAAATCGTGTTGAGGCCAACGCCCATAATGCGGGCTGTTGCCCGGCATCCA
  ACGCCATTCATGGCCATATCAATGATTTTCTGGTGCGTACCGGGTTGAGAAGCGGTGTAAGTGAACTGCAGTTGCCATGTTTTACGGCAGTGAGAGCAGA
  GATAGCGCTGATGTCCGGC
  >2
  ATGGCGACGCTGGGGCGTCTTATGAGCCTGCTGTCACCCTTTGACGTGGTGATATGGATGACGGATGGCTGGCCGCTGTATGAATCCCGCCTGAAGGGAA
  AGCTGCACGTAATCAGCAAGCGATATACGCAGCGAATTGAGCGGCATAACCTGAATCTGAGGCAGCACCTGGCACGGCTGGGACGGAAGTCGCTGTCGTT
  CTCAAAATCGGTGGAGCTGCATGACAAAGTCATCGGGCATTATCTGAACATAAAACACTATCAATAAGTTGGAGTCATTACC
  >3
  GTGCTTTTGCCGTTACGCACCACCCCGTCAGTAGCTGAACAGGAGGGACAGCTGATAGAAACAGAAGCCACTGGAGCACCTCAAAAACACCATCATACAC
  TAAATCAGTAAGTTGGCAGCATCACC

Dropping the kmer threshold to 10 and assembling those kmers finds 11 repeat
sequences, one of length 770 bp and one of length 1195 bp.

Let's now find the high-count kmers common to EC931 and MSB1_4I-sc-2280412
and assemble those.

.. code-block:: shell

  % meryl \
      print \
      intersect \
        [at-least 48 ec.meryl] \
        [at-least 55 ms.meryl] \
    | \
    perl $MERYL/scripts/greedy-assemble-kmers.pl 

  >1
  GCCGGACATCAGCGCTATCTCTGCTCTCACTGCCGTAAAACATGGCAACTGCAGTTCACTTACACCGCTTCTCAACCCGGTACGCACCAGAAAATCATTG
  ATATGGCCATGAATGGCGTTGGATGCCGGGCAACAGCCCGCATTATGGGCGTTGGCCTCAACACGATTTTACGTCACTTAAAAAACTCAGGCCGCAGTCG
  GTAACCTCGCGCATACAGCCGGGCAGTGACGTCATCGTCTGCGCGGAAATGGACGAACAGTGGGGCTATGTCGGGGCTAAATCGCGCCAGCGCTGGCTGT
  TTTACGCGTATGACAGGCT
  >2
  GGTGATGCTGCCAACTTACTGATTTAGTGTATGATGGTGTTTTTGAGGTGCTCCAGTGGCTTCTGTTTCTATCAGCTGTCCCTCCTGTTCAGCTACTGAC
  GGGGTGGTGCGTAACGGCAAAAGCAC
  >3
  ATGGCGACGCTGGGGCGTCTTATGAGCCTGCTGTCACCCTTTGACGTGGTGATATGGATGACGGATGGCTGGCCGCTGTATGAATCCCGCCTGAAGGGAA
  AGCTGCACGTAATCAGCAAGCGATATACGCAGCGAATTGAGCGGCATAACCTGAATCTGAGGCAGCACCTGGCACGGCTGGGACGGAAGTCGCTGTCGTT
  CTCAAAATCGGTGGAGCTGCATGACAAAGTCATCGGGCATTATCTGAACATAAAACACTATCAATAAGTTGGAGTCATTACC

Aside from a strand difference caused by the assembler, they're the same!

We'll stop the quick start here, before we barrel out of control into
repeats.
