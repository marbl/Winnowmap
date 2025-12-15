.. _quick-start:

Quick Start
===========

To give a quick introduction to meryl2, we'll use three random Escherichia
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
  :caption: Downloading data.
  :linenos:

  curl -LRO ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/522/225/GCA_014522225.1_ASM1452222v1/GCA_014522225.1_ASM1452222v1_genomic.fna.gz
  curl -LRO ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/165/095/GCA_002165095.2_ASM216509v2/GCA_002165095.2_ASM216509v2_genomic.fna.gz
  curl -LRO ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/071/835/GCA_905071835.1_MSB1_4I/GCA_905071835.1_MSB1_4I_genomic.fna.gz

  mkdir data

  mv -i GCA_014522225.1_ASM1452222v1_genomic.fna.gz data/ec.fna.gz
  mv -i GCA_002165095.2_ASM216509v2_genomic.fna.gz  data/sc.fna.gz
  mv -i GCA_905071835.1_MSB1_4I_genomic.fna.gz      data/ms.fna.gz

First, lets count the kmers in each and save the results in meryl databases.

.. code-block:: none
  :caption: Counting k-mers in a single file.
  :linenos:

  % meryl2 -k 42 count data/ec.fna.gz output=ec.meryl

  Found 1 command tree.

  |- TREE 0: action #1 count kmers
     |> database 'ec.meryl'
     |- SET value to that of the kmer in the first input
     |- SET label to first -- constant 0
     ^- INPUT @1: sequence file 'data/ec.fna.gz'
  |- TREE 0 ends.

  Counting 4475 (estimated) thousand canonical 42-mers from 1 input file:
      sequence-file: data/ec.fna.gz

  [...]

We can show the k-mers in the database by printing the output to the screen
(or to a file with `print=42-mers.txt`).

.. code-block:: none
  :caption: Displaying k-mers.
  :linenos:

  % meryl2 -t 1 -Q print ec.meryl | head

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

(Messy detail: Option `-Q` turns off the report of the command tree.  Option
`-t 1` is supplied to prevent meryl from processing its 64 data chunks in
parallel; for the print operation this results in the output kmers being
unsorted.  See {INSERT LINK TO THREADS AND PRINT HERE}.)

Let's now count the other two databases and combine the kmers from all three
genomes into one database, all in one command.

.. code-block:: none
  :caption: Counting two FASTA files and merging with a third database.
  :linenos:

  % meryl2 \
      union-sum \
      output=union.meryl \
      [count data/sc.fna.gz output=sc.meryl] \
      [count output=ms.meryl \
             data/ms.fna.gz] \
      ec.meryl

  Found 1 command tree.

  |- TREE 0: action #1 filter kmers
     |> database 'union.meryl'
     |- SET value to the sum of all kmers and constant 0
     |- SET label to or -- constant 0
     |- FILTER 1
        |- EMIT if <index filter not described>
     ^- INPUT @1: action #2 count kmers
        |> database 'sc.meryl'
        |- SET value to that of the kmer in the first input
        |- SET label to first -- constant 0
        ^- INPUT @1: sequence file 'data/sc.fna.gz'
     ^- INPUT @2: action #3 count kmers
        |> database 'ms.meryl'
        |- SET value to that of the kmer in the first input
        |- SET label to first -- constant 0
        ^- INPUT @1: sequence file 'data/ms.fna.gz'
     ^- INPUT @3: meryl database 'ec.meryl'
  |- TREE 0 ends.

  [...]

There's a lot going on here.  There are three `actions` (lines 11, 17 and
22), three `output` files (lines 12, 18 and 23), and three `inputs` (lines
17, 22, and 27).  Each action decribes how to combine the k-mers in its
inputs into an output k-mer, then passes those k-mers to its destination
actions.  Here, `action #1` takes input k-mers from `INPUT @1` (which itself
is `action #2`), `INPUT @2` (`action #3`) and `INPUT @3` (which is a
pre-computed database of k-mers).

The `union-sum` action in the command line is `action #1` in the tree.  Lines
13 and 14 describe how this action is computing the output value and label of
each k-mer.  Line 16 will {EVENTUALLY} describe the conditions that must be
met for a k-mer to be output.  For `union`, there is only one condition: the
k-mer must be in at least one input.

The end result of this is to independently count kmers in `sc.fna.gz` and
`ms.fna.gz`, writing the kmers from each to outpout databases {SEE COUNTING}
`sc.meryl` and `ms.meryl`, respectively.  When the counting operations are
done, those two new databases and the third pre-computed database are sent as
input to the first action which will combine all k-mers, summing their
values, into one output.

Notice that the `-k 42` option is not present.  Meryl will determine the
k-mer size from the `ec.meryl` input database and use that for all
operations.  However, if a k-mer size is supplied, it must match the size of
ALL input databases, and ALL input databases must have k-mers of the same
size.

A meryl database also stores the histogram of kmer values.  This can be
displayed:

.. code-block:: none
  :caption: A k-mer value histogram.
  :linenos:

  % meryl2 -Q histogram ec.meryl

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
(Escherichia coli SCEC020022, `sc.meryl`) or a potential 64 copy repeat
(Escherichia coli MSB1_4I-sc-2280412, `ms.meryl`).  Let's now extract those
kmers and see where they are on the genomes.

.. code-block:: none
  :caption: Extracting high-value k-mers.
  :linenos:

  % meryl2 print=ec-repeats.dump at-least 48 ec.meryl output=ec-repeats.meryl

  Found 1 command tree.

  |- TREE 0: action #1 filter kmers
     |> database 'ec-repeats.meryl'
     |> text file 'ec-repeats.dump'
     |- SET value to that of the kmer in the first input
     |- SET label to first -- constant 0
     |- FILTER 1
        |- EMIT if output kmer value     is-more-or-equal constant value 48
     ^- INPUT @1: meryl database 'ec.meryl'
  |- TREE 0 ends.

  [...]

  % wc -l ec-repeats.dump
       604 ec-repeats.dump

The `meryl-lookup` tool compares FASTA/FASTQ sequences against a meryl
database (or several databases) and generates various reports about how the
k-mers in the database(s) "paint" onto the input sequences.  We'll use it to
generate a bed file of the bases covered by kmers in our E.coli  database.

.. code-block:: none
  :caption: Finding runs of high-value k-mers in a genome.
  :linenos:

  % meryl2-lookup -sequence data/ec.fna.gz -mers ec-repeats.meryl -bed-runs > ec-repeats.bed
  --
  -- Estimating memory usage for 'ec-repeats.meryl'.
  --
  
   p       prefixes             bits gigabytes (allowed: 63 GB)
  -- -------------- ---------------- ---------
   2              4            53408     0.000
   3              8            53060     0.000
   4             16            52968     0.000
   5             32            53388     0.000
   6             64            54832     0.000 (smallest)
   7            128            58324     0.000
   8            256            65912     0.000
   9            512            81692     0.000 (faster)
  10           1024           113856     0.000
  11           2048           178788     0.000
  12           4096           309256     0.000
  13           8192           570796     0.000
  -- -------------- ---------------- ---------
                604 total kmers
  
  --
  -- Minimal memory needed: 0.000 GB
  -- Optimal memory needed: 0.000 GB  enabled
  -- Memory limit           63.907 GB
  --
  --
  -- Loading kmers from 'ec-repeats.meryl' into lookup table.
  --
  
  For 604 distinct 42-mers (with 9 bits used for indexing and 75 bits for tags):
      0.000 GB memory for kmer indices -          512 elements 64 bits wide)
      0.000 GB memory for kmer tags    -          604 elements 75 bits wide)
      0.000 GB memory for kmer values  -          604 elements  6 bits wide)
      0.000 GB memory
  
  Will load 604 kmers.  Skipping 0 (too low) and 0 (too high) kmers.
  Allocating space for 16732 suffixes of 75 bits each -> 1254900 bits (0.000 GB) in blocks of 32.000 MB
                       16732 values   of 6 bits each -> 100392 bits (0.000 GB) in blocks of 32.000 MB
  Loaded 604 kmers.  Skipped 0 (too low) and 0 (too high) kmers.
  -- Opening input sequences 'data/ec.fna.gz'.
  -- Opening output file '-'.
  Bye!

  % head ec-repeats.bed 
  CP049118.1      46087   46370
  CP049118.1      46409   46729
  CP049118.1      46729   46856
  CP049118.1      140430  140592
  CP049118.1      140592  140713
  CP049118.1      140752  141072
  CP049118.1      141072  141199
  CP049118.1      270969  271131
  CP049118.1      271131  271252
  CP049118.1      271291  271611

From this we see that there is not a single 52-copy repeat, but several
shorter repeats.  Curiously, there are several instances of a 447 base repeat
with single base differences:

.. code-block:: none
  :caption: Curiously similar repeats.
  :linenos:

  CP049118.1      3782667 3782794
  CP049118.1      3782794 3783114

  CP049118.1      3762937 3763257
  CP049118.1      3763257 3763384

Though this isn't really part of meryl, the high-count kmers can be passed to
a greedy assembler with nice results (the greedy assembler is included in the
meryl source code, but isn't installed in the binary directory).

.. code-block:: none
  :caption: Greedily assembling high-value k-mers.
  :linenos:

  % perl $MERYL/scripts/greedy-assemble-kmers.pl < ec-repeats.dump
  >1
  AGCCTGTCATACGCGTAAAACAGCCAGCGCTGGCGCGATTTAGCCCCGACATAGCCCCACTGTTCGTCCATTTCCGC
  GCAGACGATGACGTCACTGCCCGGCTGTATGCGCGAGGTTACCGACTGCGGCCTGAGTTTTTTAAGTGACGTAAAAT
  CGTGTTGAGGCCAACGCCCATAATGCGGGCTGTTGCCCGGCATCCAACGCCATTCATGGCCATATCAATGATTTTCT
  GGTGCGTACCGGGTTGAGAAGCGGTGTAAGTGAACTGCAGTTGCCATGTTTTACGGCAGTGAGAGCAGAGATAGCGC
  TGATGTCCGGC
  >2
  ATGGCGACGCTGGGGCGTCTTATGAGCCTGCTGTCACCCTTTGACGTGGTGATATGGATGACGGATGGCTGGCCGCT
  GTATGAATCCCGCCTGAAGGGAAAGCTGCACGTAATCAGCAAGCGATATACGCAGCGAATTGAGCGGCATAACCTGA
  ATCTGAGGCAGCACCTGGCACGGCTGGGACGGAAGTCGCTGTCGTTCTCAAAATCGGTGGAGCTGCATGACAAAGTC
  ATCGGGCATTATCTGAACATAAAACACTATCAATAAGTTGGAGTCATTACC
  >3
  GTGCTTTTGCCGTTACGCACCACCCCGTCAGTAGCTGAACAGGAGGGACAGCTGATAGAAACAGAAGCCACTGGAGC
  ACCTCAAAAACACCATCATACACTAAATCAGTAAGTTGGCAGCATCACC

Dropping the kmer threshold to 10 (`at-least 10`) and assembling those repeat
k-mers finds 11 repeat sequences, one of length 770 bp and one of length 1195
bp.

For our final example, let's find the high-value k-mers common to EC931 and
MSB1_4I-sc-2280412 and assemble those.

.. code-block:: none
  :caption: Greedily assembling medium-and-high-value k-mers.
  :linenos:

  % meryl2 -Q \
      print \
      intersect \
        [at-least 48 ec.meryl] \
        [at-least 55 ms.meryl] \
    | \
    perl $MERYL/scripts/greedy-assemble-kmers.pl 
  >1
  GGTAATGACTCCAACTTATTGATAGTGTTTTATGTTCAGATAATGCCCGATGACTTTGTCATGCAGCTCCACCGATT
  TTGAGAACGACAGCGACTTCCGTCCCAGCCGTGCCAGGTGCTGCCTCAGATTCAGGTTATGCCGCTCAATTCGCTGC
  GTATATCGCTTGCTGATTACGTGCAGCTTTCCCTTCAGGCGGGATTCATACAGCGGCCAGCCATCCGTCATCCATAT
  CACCACGTCAAAGGGTGACAGCAGGCTCATAAGACGCCCCAGCGTCGCCAT
  >2
  AGCCTGTCATACGCGTAAAACAGCCAGCGCTGGCGCGATTTAGCCCCGACATAGCCCCACTGTTCGTCCATTTCCGC
  GCAGACGATGACGTCACTGCCCGGCTGTATGCGCGAGGTTACCGACTGCGGCCTGAGTTTTTTAAGTGACGTAAAAT
  CGTGTTGAGGCCAACGCCCATAATGCGGGCTGTTGCCCGGCATCCAACGCCATTCATGGCCATATCAATGATTTTCT
  GGTGCGTACCGGGTTGAGAAGCGGTGTAAGTGAACTGCAGTTGCCATGTTTTACGGCAGTGAGAGCAGAGATAGCGC
  TGATGTCCGGC
  >3
  GTGCTTTTGCCGTTACGCACCACCCCGTCAGTAGCTGAACAGGAGGGACAGCTGATAGAAACAGAAGCCACTGGAGC
  ACCTCAAAAACACCATCATACACTAAATCAGTAAGTTGGCAGCATCACC

Aside from a strand and sequence order difference caused by the assemlber,
they're the same!

.. code-block:: none
  :caption: The same!
  :linenos:

  % minimap2 --eqx -Y -c ec.fasta ec+ms.fasta
  1 282 0 282 - 2 282 0 282 282 282 60 NM:i:0 ms:i:564 AS:i:564 nn:i:0 tp:A:P cm:i:49 s1:i:274 s2:i:0 de:f:0 rl:i:0 cg:Z:282=
  2 319 0 319 + 1 319 0 319 319 319 60 NM:i:0 ms:i:638 AS:i:638 nn:i:0 tp:A:P cm:i:58 s1:i:311 s2:i:0 de:f:0 rl:i:0 cg:Z:319=
  3 126 0 126 + 3 126 0 126 126 126 60 NM:i:0 ms:i:252 AS:i:252 nn:i:0 tp:A:P cm:i:16 s1:i:115 s2:i:0 de:f:0 rl:i:0 cg:Z:126=

.. code-block:: none
  :caption: But rearranged.
  :linenos:

  % minimap2 -x sr --eqx -Y -c data/ec.fna.gz ec.fasta
  1 319 0 319 + CP049118.1 4767973 2542965 2543284 319 319 0 NM:i:0 ms:i:638 AS:i:638 nn:i:0 tp:A:P cm:i:44 s1:i:304 s2:i:304 de:f:0 rl:i:0 cg:Z:319=
  2 282 0 282 - CP049118.1 4767973 3430238 3430520 282 282 0 NM:i:0 ms:i:564 AS:i:564 nn:i:0 tp:A:P cm:i:39 s1:i:266 s2:i:266 de:f:0 rl:i:0 cg:Z:282=
  3 126 0 126 - CP049118.1 4767973 4539117 4539243 126 126 0 NM:i:0 ms:i:252 AS:i:252 nn:i:0 tp:A:P cm:i:19 s1:i:123 s2:i:123 de:f:0 rl:i:0 cg:Z:126=

  % minimap2 -x sr --eqx -Y -c data/ms.fna.gz ec.fasta
  1 319 0 319 - LR898874.1 5278133   11725   12044 319 319 0 NM:i:0 ms:i:638 AS:i:638 nn:i:0 tp:A:P cm:i:44 s1:i:304 s2:i:304 de:f:0 rl:i:0 cg:Z:319=
  2 282 0 282 - LR898874.1 5278133 1418500 1418782 282 282 0 NM:i:0 ms:i:564 AS:i:564 nn:i:0 tp:A:P cm:i:39 s1:i:266 s2:i:266 de:f:0 rl:i:0 cg:Z:282=
  3 126 0 126 - LR898874.1 5278133 1510654 1510780 126 126 0 NM:i:0 ms:i:252 AS:i:252 nn:i:0 tp:A:P cm:i:19 s1:i:123 s2:i:123 de:f:0 rl:i:0 cg:Z:126=

We'll stop the quick-start here, before we barrel out of control into
repeats.
