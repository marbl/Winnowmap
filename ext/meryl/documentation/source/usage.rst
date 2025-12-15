.. _usage:


====================
Command Line Options
====================

DEBUG OPTIONS
-------------

'dumpIndex <merylDB>'
 - reports encoding parameters
        > meryl2 dumpIndex basic-dh-1.meryl
        Opened 'basic-dh-1.meryl'.
          magic          0x646e496c7972656d33302e765f5f7865 'merylIndex__v.03'
          prefixSize     10
          suffixSize     104
          numFilesBits   6 (64 files)
          numBlocksBits  4 (16 blocks)

'dumpFile <merylDB/0x######>>'
 - reports encoding blocks

        > meryl2 dumpFile basic-dh-1.meryl/0x000000

            prefix    blkPos    nKmers
        ---------- --------- ---------
        0x00000000         0     22363
        0x00000001    345448     16486
        0x00000002    601248     17345
        [...]

                    prefix   nKmers kCode uBits bBits                 k1 cCode                 c1                 c2
        ------------------ -------- ----- ----- ----- ------------------ ----- ------------------ ------------------
        0x0000000000000000    22363     1    15    89 0x0000000000000000     1 0x0000000000000000 0x0000000000000000
        0x0000000000000001    16486     1    15    89 0x0000000000000000     1 0x0000000000000000 0x0000000000000000
        0x0000000000000002    17345     1    15    89 0x0000000000000000     1 0x0000000000000000 0x0000000000000000
        [...]

         kmerIdx prefixDelta      prefix |--- suffix-size and both suffixes ---|    value
        -------- ----------- ----------- -- ---------------- -- ---------------- --------
               0          41 00000000029 25 00000000005026f2 64 f6a920a6cc0a36dd        1
               1           4 0000000002d 25 00000000010001f5 64 b1ec24342aaa4faa        1
               2           0 0000000002d 25 00000000015276f4 64 faaa2c6000702827        1
        [...]



----------------------------------------

Take whole command line and merge into a single string,
then pass that to regex to pull out options with spaces.

But caution!  Files with spaces mess this up.
  meryl count "input one.fasta" "input two.fasta"

So we need to use some different separator in the regex.
Or somehow remember that specific words are possibly
files?
  meryl \1 count \1 input one.fasta \1 input two.fasta
  meryl \1 s:v \1 greater-than \1 1
But need to allow matches to both \1 and space, so
that an quoted terms ("s:v > 1") work too.

Or maybe just surround everything with quotes?  But then
what about files with quotes in them?  Nah, that won't work.

  meryl \1 s:v:greater-than 1 \1 input database.meryl

but still left with the confusion over
  meryl output:statistics     in.meryl
  meryl output:statistics out in.meryl
so need to enfore an = before parameters
  meryl output:statistics=out   in.meryl
  meryl output:statistics= out  in.meryl
  meryl output:statistics = out in.meryl

  meryl s:v greater-than 1     in.meryl output x
  meryl s:v:greater-than 1 i:d=in.meryl o:d=x
  meryl s:v greater-than 1     in.meryl output=x

 "o:s" "=" "stats out" "i:d=in db.meryl"
  o:s\1=\1stats\2out\1i:d=in\2db.meryl

 "o:s = stats out" "i:d=in db.meryl"
  o:s\2=\2stats\2out\1i:d=in\2db.meryl

  o:s = "stats out" i:d= "in db.meryl"
  o:s\1=\1stats\2out\1i:d=\1in\2db.meryl

  output:statistics \s* = \s* \w+ \W+   #  out stats to file
  \w+ \W+                               #  input database
  o:d \s* = \s* \w+ \W+                 #  out database

  \1 == \W == word boundary (or newline)
  \2 ==       space internal to token
  \s == internal space or word boundary - \2 or \1
  \w == anything but word boundary

  s:v:@1:lt:@4
  s:v @1 lt @4
  s:v:@1lt@4
  s:v:@1-less-than-@4

  if input program from file, need to find quotes and parse
      input  file  .fasta  -> \Winput\Wfile\W.fasta\W
     "input  file  .fasta" -> \Winput\2file\2.fasta\W
     "input\"file\".fasta" -> \Winput"file".fasta\W

  if fail to match, report rest of string up to first \W
   - last case is inout database or file, which must exist
     in filesystem

  really need to define a grammar for this.

----------------------------------------

OPTIONS
  'compress'
  'n=N'
  'count-suffix=ACGT'
  'd=D'
  'distinct=D'
  'f=F'
  'word-frequency=F'
  't=T'
  'threshold=T'
  'segment=a/b'

ALIASES that take parameters
  distinct=
  word-freq=
  word-frequency=
  threashold=

ALIASES
  union
  union-min
  union-max
  union-sum

  intersect
  intersect-min
  intersect-max
  intersect-sum

  subtract
  difference

  less-than
  greater-than
  at-least
  at-most
  equal-to
  not-equal-to
  increase
  decrease
  multiply
  divide
  divide-round
  modulo

assign:value=
  #X       - constant X
  @X       - the value of the k-mer in the Xth input
  first    - the value of the k-mer in the 1st input
  selected - the value of the k-mer used to set the label;
             if multiple are k-mers were used, the value of the first
  count    - the number of inputs the k-mer is present in

  min(#X)
  max(#X)
  add(#X) and sum(#X)
  sub(#X) and dif(#X)
  mul(#X)
  div(#X) and divzero(#X)
  mod(#X)
  rem(#X)

assign:label=
  #X         - constant X
  @X         - the label of the k-mer in the Xth input
  first      - the label of the k-mer in the 1st input
  selected   - the label of the k-mer used to set the value;
               if multiple are k-mers were used, the label of the first

  min(#X)    - the minimum of all labels interpreted as an unsigned integer
  max(#X)    - the maximum of all labels interpreted as an unsigned integer

  lightest(#X) - the label with the fewest bits set
  heaviest(#X) - the label with the most bits set

  and(#X)    - the bitwise AND of all input labels
  or(#X)     - the bitwise OR  of all input labels
  xor(#X)    - the bitwise XOR of all input labels

  complement - the bitwise complement of the first input label
               equivalent to xor(#1111...11)

  difference(#X) - the label in the first input bitwise minus
                   all other labels

  shift-left(#X)
  shift-right(#X)
  rotate-left(#X)
  rotate-right(#X)

  GENERAL LABEL ASSIGNMENT

----------
select:value:
  ARG1      OP    ARG2
  @n              @n         - value in the nth input
  first           first      - value in the 1st input
  #n or n         #n or n    - constant n
  <omitted>                  - output kmer
  output                     - output kmer
                  distinct=f
                  word-freq=f
                  threshold=n

  OPERATORS
            ==, =, eq, equal
            !=, <>, ne, neq
            <=, le, less-than-or-equal, at-most
            >=, ge, more-than-or-equal, at-least
            <,  lt, less-than
            >,  gt, greater-than

  allow spaces, or underscores or dashes?

  s:v:output<4     - ok
  s:v:outputlt4    - this is awkward
  s:v:@1:lt:@4
  s:v:@1-lt-@4
  s:v:output-<-4   - this is awkward
  s:v:output-lt-4  - 
  s:v:<4           - shell problems, needs quotes anyway
  s:v:less-than4   - this is awkward
  s:v:less-than-4
  s:v less-than 4
  s:v first == 4
  s:v first-==-4   - this is awkward
  s:v first equal 4
  s:v first-equal-to-4
  s:v "first equal 4"
  s:v "first <= 4"

----------
select:label:













====================
Command Line Options
====================

meryl2 Global Options
---------------------

In meryl2, global options apply to every action.  Global options are expressed as
command line switches.  K-mer size, label size, memory usage, number of
threads, and simple sequence reductions are all global options.

Processing is specified as a list of ``actions`` and are described elsewhere.

.. code-block:: none
  :caption: meryl2 usage.
  :linenos:

  KMER and LABEL SIZE:
    -k K             Set kmer size to K (6 <= K <= 64).  Legacy format "k=K".
    -l L             Set label size to L bits (0 <= L <= 64).

  SIMPLE SEQUENCE REDUCTION:
    --compress       Homopolymer compress input sequence before forming kmers.
    --ssr ssr-spec   Apply 'ssr-spec' to input sequence before forming kmers.

  RESOURCE USAGE:
    -m M             Use up to M GB memory for counting.  Legacy format "memory=M".
    --memory M       Use up to M GB memory for counting.

    -t        T      Use T threads for counting and processing.
    --threads T      Use T threads for counting and processing.

  LOGGING:
    -V[V[V[...]]]    Increae verbosity by the length of this option.
    -Q               Be absolutely silent.
    -P               Show progress.
    -C               Show processing tree and stop.

  USAGE:
    -h               Display command line help.
    --help           Display command line help.
    help             Display command line help.

Obsolete forms of some of these are still allowed.  The kmer size can be set
with `k=<kmer-size>`, memory limits with `memory=<memory-in-gigabytes>`,
thread usage with `threads=<thread-count>`.

Obsolete option `-E` has been removed.  This used to estimate the size of an
imput that could be counted in some given memory size.

.. code-block:: none
  :caption: meryl2 debugging actions.
  :linenos:

  dumpIndex <meryl-database>
    Report the parameters used for storing data in this meryl database.

  dumpFile <meryl-database-file>
    Dump the raw data from a merylData file in a meryl database.


meryl2-lookup Usage
-------------------


.. code-block:: none
  :caption: meryl2-lookup usage.
  :linenos:

  usage: meryl2-lookup <report-type> \
           -sequence <input1.fasta> [<input2.fasta>] \
           -output   <output1>      [<output2>] \
           -mers     <input1.meryl> [<input2.meryl>] [...] [-estimate] \
           -labels   <input1name>   [<input2name>]   [...]

    Compare kmers in input sequences against kmers in input meryl databases.

    Input sequences (-sequence) can be FASTA or FASTQ, uncompressed, or
    compressed with gzip, xz, or bzip2.

    Report types:

    -bed:
       Generate a BED format file showing the location of kmers in
       any input database on each sequence in 'input1.fasta'.
       Each kmer is reported in a separate bed record.

    -bed-runs:
       Generate a BED format file showing the location of kmers in
       any input database on each sequence in 'input1.fasta'.
       Overlapping kmers are combined into a single bed record.

    -wig-count:
       Generate a WIGGLE format file showing the multiplicity of the
       kmer starting at each position in the sequence, if it exists in
       an input kmer database.

    -wig-depth:
       Generate a WIGGLE format file showing the number of kmers in
       any input database that cover each position in the sequence.

    -existence:
       Generate a tab-delimited line for each input sequence with the
       number of kmers in the sequence, in the database and common to both.

    -include:
    -exclude:
       Copy sequences from 'input1.fasta' (and 'input2.fasta') to the
       corresponding output file if the sequence has at least one kmer
       present (include) or no kmers present (exclude) in 'input1.meryl'.

  Run `meryl2-lookup <report-type> -help` for details on each method.

