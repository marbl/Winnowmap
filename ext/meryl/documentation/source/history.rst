.. _history:

Software Background
====================

Meryl was first imagined in late-2000/early-2001/mid-2001 while BPW as at `Celera Genomics
<https://en.wikipedia.org/wiki/Celera_Corporation>`_.

A 3 June 2001 journal entry is titled "merConting the genome", but makes no
mention of motivation for doing so.  The same journal has an 18 June entry
outlining the algorithm, and this entry also strongly hints the application
was to build a fast 20-mer to assembly coordinate lookup table.  The "Meryl"
name is mentioned on 16 August 2002, desiring to increment/decrement the
counts.  There is even *source code* dating back to 12 February 2001 deep in
the dusty archives; as that's the oldest archive found, meryl is doubtless a
little bit older yet.  It was originally called "merMaid" or "merCounter",
the "Meryl" name first appears on 11 June 2001, and the code looks reasonably
mature.

Meryl was definitely used twice in the assembly of Anopheles gambiae [`Holt,
et al. 2002 <https://science.sciencemag.org/content/298/5591/129>`_], once to
find k-mers that occur frequently and exclude them from seeding overlaps
(which is how it is still used in Canu [`Koren and Walenz 2017
<http://doi.org/10.1101/gr.215087.116>`_]) and once to estimate repeat
content of the assembly.  It is mentioned, anonymously, the paper:

  By counting the number of times each 20-nucleotide oligomer in the Anopheles and Drosophila assemblies appeared in its corre- sponding whole-genome shotgun data, we con- firmed that simple repeats are not expanded in Anopheles

The supplement elaborates a bit:

  **METHODS FOR ESTIMATING REPEAT CONTENT:**

  The consensus sequence of scaffolded contigs was analyzed for repeat content as follows. For each 20mer of the consensus, we counted the number of times that that 20mer appeared in the set of approximately 4.5 million sequence reads. Out of 262.8 M consensus 20mers, 26.3M (10.0%) occurred more than 50 times (5 times what is expected for unique sequence, given 10- fold sequence coverage). For comparison, of 133.5M consensus 20mers from a recent reassembly of the D. melanogaster genome that was done using the same version of the assembly software that was used for mosquito , 13.0M (9.7%) were observed more than 60 times (5 times more than expected, given 12-fold sequence coverage of the Drosophila genome). Using a count cutoff of 1.5 times expected or 50 times expected gives a similar picture: the repetative fractions of the A. gambiae and D. melanogaster genomes are not strikingly different.

Unfortunately, the earliest versions of the software are not online, likely
because they were never stored in a version control system.  The earliest
online version is found in the `initial commit
<https://sourceforge.net/p/kmer/code/4/>`_ of the `kmer subversion repository
<http://kmer.sourceforge.net>`_ on 2 January 2003.  The `main function in
meryl.C <https://sourceforge.net/p/kmer/code/4/tree//trunk/meryl/meryl.C>`_
shows it supporting the "binary and mathematical" operations, including
"increment/decrement" from the 16 August 2002 journal entry.

The `Celera Assembler <http://wgs-assembler.sourceforge.net/wiki/index.php?title=Main_Page>`_
subversion repository also has a relatively early version.  CA was publicly
released to SourceForge on 14 April 2004 (`revision 4 has the good stuff
<https://sourceforge.net/p/wgs-assembler/svn/4>`_).  Use of meryl is
mentioned in `wga.pl <https://sourceforge.net/p/wgs-assembler/svn/4/tree//trunk/example/wga.pl>`_.
The meryl code and a CA-specific function to access sequence data are `nicely
organized <https://sourceforge.net/p/wgs-assembler/svn/4/tree/trunk/src/AS_MER/>`_, the
bulk of meryl in one file (`AS_MER_meryl.cc <https://sourceforge.net/p/wgs-assembler/svn/4/tree/trunk/src/AS_MER/AS_MER_meryl.cc>`_).

The offline BPW archives contain several early versions, including one just nine days
after the first journal entry mentioned above:

.. code-block:: none

  -rw-r--r--  1 bri  bri  5964 Jun 10  2001 ne-20010611-1812/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri   180 Jul  3  2001 ne-20010703-1747/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  1322 Jul  6  2001 ne-20010706-1813/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4132 Jul  9  2001 ne-20010709-1647/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4294 Jul 10  2001 ne-20010710-1842/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4675 Jul 11  2001 ne-20010712-1907/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4775 Jul 18  2001 ne-20010718-1745/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4706 Jul 23  2001 ne-20010723-1800/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4850 Aug 15  2001 ne-20010817-1800/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4988 Feb 27  2005 ne-20010831-1745/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  4988 Aug 30  2001 ne-20010912-1829/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  5130 Feb 27  2005 ne-20010928-1826/near-identity/meryl/meryl.C
  -rw-r--r--  1 bri  bri  5707 Feb 27  2005 ne-20011109-1837/near-identity/meryl/meryl.C

(I'm not sure where the 2005 time stamps came from; I vaguely remember
replacing duplicate files with symlinks before some of these directories were
tar'd and compressed.)

Peeking into that first version, we find:

.. code-block:: none

  d:.../ne-20010611-1812/near-identity/meryl> ls -l
  total 88
  -rw-r--r--  1 bri  bri   2187 Jun 10  2001 1.C
  -rw-r--r--  1 bri  bri   1232 Jun 10  2001 2.C
  -rw-r--r--  1 bri  bri   1582 Jun 10  2001 3.C
  -rw-r--r--  1 bri  bri   2496 Jun 10  2001 4.C
  -rw-r--r--  1 bri  bri    538 Jun 11  2001 Makefile
  -rwxr-xr-x  1 bri  bri  50432 Jun 10  2001 meryl
  -rw-r--r--  1 bri  bri   5964 Jun 10  2001 meryl.C
  -rw-r--r--  1 bri  bri   1725 Jun 10  2001 meryl.H

Is that a binary?

.. code-block:: none

  d:.../ne-20010611-1812/near-identity/meryl> file meryl
  meryl: ELF 32-bit LSB executable, \
         Intel 80386, \
         version 1 (FreeBSD), \
         dynamically linked, \
         interpreter /usr/libexec/ld-elf.so.1, \
         for FreeBSD 4.1, \
         not stripped

The algorithm is quite simple, but the details are complicated: it simply
builds a list of all the kmers in the input sequences, sorts, then counts how
many times each kmer is in the list and write that to the output file.

The details are to split the kmer into a prefix and a suffix.  The prefix
points to a bucket into which all the suffixes are listed.  Once all kmers
are stored, each bucket is sorted then scanned to count the kmers.
Bit-packed integers are used throughout to minimize memory usage.

In all it's embarassing glory, here is the first ever version of meryl
(omitting the functions that do the actual work).

.. code-block:: c++

  #include "meryl.H"

  //  theSeq is a compressed sequence.  Three bits per character, the first bit
  //  telling us if the next two bits are valid sequence.
  //
  u64bit  *theSeq       = 0L;
  u64bit   theSeqLen    = 0;
  u64bit   numberOfMers = 0;

  //  theCounts is a packed-word array of the number of mers
  //  that all have the same TABLEBITS bits on the left side of
  //  the mer.
  //
  //  theMers is a list of the other bits in the mers.
  //
  u64bit  *theCounts       = 0L;
  u64bit   theCountsWords  = 0;
  u64bit  *theMers         = 0L;
  u64bit   theNumberOfMers = 0;

  u64bit  approxMers = 1073741824;
  u64bit  merSize    = 20;
  u64bit  merMask    = u64bitMASK(merSize);

  u64bit  tblBits = 20;
  u64bit  tblMask = u64bitMASK(tblBits);

  u64bit  chkBits = 0;
  u64bit  chkMask = 0;

  u64bit  bitsPerIndex = 0;

  bool    beVerbose           = true;
  bool    doReverseComplement = false;

  char   *seqFileName = 0L;
  FILE   *seqFile     = 0L;

  void
  usage(char *name) {
    fprintf(stderr, "usage: %s [-m merSize] [-v] [-r] [-f seqFile] [-s seqFileSize] [-h]\n", name);
    fprintf(stderr, "       -m   Sets the size of mer to count.  Must be less than 32.\n");
    fprintf(stderr, "       -v   Be noisy about it.\n");
    fprintf(stderr, "       -r   Also counts the reverse complement.\n");
    fprintf(stderr, "       -f   File to count mers in.  MultiFastA\n");
    fprintf(stderr, "       -s   Upper limit on the number of mers.\n");
    fprintf(stderr, "       -h   This.\n");
  }


  int
  main(int argc, char **argv) {

    int   arg = 1;
    while (arg < argc) {
      switch (argv[arg][1]) {
        case 'v':
          beVerbose = true;
          break;
        case 'r':
          doReverseComplement = true;
          break;
        case 'f':
          arg++;
          seqFileName = argv[arg];
          break;
        case 'm':
          arg++;
          merSize = atol(argv[arg]);
          break;
        case 's':
          arg++;
          approxMers = atol(argv[arg]);
          break;
        case 'h':
          usage(argv[0]);
          exit(1);
          break;
        default:
          usage(argv[0]);
          fprintf(stderr, "Unknown option '%s'\n", argv[arg]);
          exit(1);
          break;
      }
      arg++;
    }

    if ((seqFileName == 0L) ||
        ((seqFileName[0] == '-') && (seqFileName[1] == 0)))
      seqFile = stdin;
    else
      seqFile = fopen(seqFileName, "r");
    if (seqFile == 0L) {
      fprintf(stderr, "Couldn't open the sequence file '%s'.\n", seqFileName);
      exit(1);
    }


    ////////////////////////////////////////////////////////////
    //
    //  Based on the approximate number of mers given, set the
    //  size of the table and number of bits per each table entry.
    //

    //  Not the greatest way to find the number of bits needed
    //  to encode approxMers, but easy
    //
    bitsPerIndex = 1;
    while (approxMers > u64bitMASK(bitsPerIndex))
      bitsPerIndex++;

    //  Determine the size of the table to reduce the memory footprint.
    //  This is computed backwards, so that we pick the best smallest
    //  table size.
    //
    u64bit   bestMem = ~u64bitZERO;
    u64bit   bestSiz =  u64bitZERO;
    u64bit   mem;

    for (u64bit siz=32; siz>18; siz-=2) {
      mem   = bitsPerIndex * (u64bitONE << siz) + approxMers * (2*merSize - siz);
      mem >>= 23;

      if (mem < bestMem) {
        bestMem = mem;
        bestSiz = siz;
      }
    }

    tblBits = bestSiz;
    tblMask = u64bitMASK(tblBits);

    chkBits = 2 * merSize - tblBits;
    chkMask = u64bitMASK(chkBits);


    fprintf(stderr, "You told me there will be no more than %lu mers.\n", approxMers);
    fprintf(stderr, "Using %2lu bits to encode positions.\n", bitsPerIndex);
    fprintf(stderr, "Using %2lu bits of the mer for the count.\n", tblBits);
    fprintf(stderr, "Using %2lu bits of the mer for the check.\n", chkBits);





    ////////////////////////////////////////////////////////////
    //
    //  Allocate and clear the counts
    //
    theCountsWords = bitsPerIndex * (u64bitONE << tblBits) / 64 + 2;

    if (beVerbose) {
      fprintf(stderr, "Allocating and clearing the counting table.\n");
      fprintf(stderr, "  %lu entries at %lu bits == %lu MB.\n",
              u64bitONE << tblBits,
              bitsPerIndex,
              theCountsWords * 8 >> 20);
    }

    theCounts      = new u64bit [ theCountsWords ];
    for (u64bit i=theCountsWords; i--; )
      theCounts[i] = u64bitZERO;


    //////////////////////////////////////////////////////////////
    //
    //
    theSeqLen    = 0;
    numberOfMers = 0;
    theSeq       = new u64bit [ 3 * (u64bitONE << bitsPerIndex) / 64 + 1 ];

    if (beVerbose)
      fprintf(stderr, "Computing the bucket sizes (pass one through the sequence).\n");
    readSequenceAndComputeBucketSize();

    if (theNumberOfMers > approxMers) {
      fprintf(stderr, "Sorry.  You need to increase your estimate of the number of mers!\n");
      fprintf(stderr, "I found %lu mers.\n", theNumberOfMers);
      exit(1);
    }


    //////////////////////////////////////////////////////////////
    //
    //
    if (beVerbose)
      fprintf(stderr, "Converting counts to offsets (%lu bits or %lu Mbuckets).\n",
              tblBits,
              (u64bitONE << tblBits) >> 20);
    convertCounts();


    //////////////////////////////////////////////////////////////
    //
    //  Allocate theMers, but don't clear them.
    //
    if (beVerbose) {
      fprintf(stderr, "Allocating space for the mers.\n");
      fprintf(stderr, "  %lu mers at %lu bits each == %lu MB.\n",
              theNumberOfMers,
              chkBits,
              chkBits * theNumberOfMers >> 23);
    }
    theMers = new u64bit [ ((chkBits * theNumberOfMers) >> 6) + 1 ];


    //////////////////////////////////////////////////////////////
    //
    //
    if (beVerbose)
      fprintf(stderr, "Filling the buckets with mers (pass two through the sequence).\n");
    fill();


    //////////////////////////////////////////////////////////////
    //
    //  Sort each bucket
    //
    if (beVerbose)
      fprintf(stderr, "Sorting buckets (%lu total).\n", u64bitONE << tblBits);
    sort();


    //
    //  Yeah, I didn't close the stupid input file.  So what?
    //
  }

The first code block reads input sequence from disk, converts each base into
a 3-bit code, the counts the size of each bucket.

.. code-block:: c++

  #include "meryl.H"

  void
  readSequenceAndComputeBucketSize(void) {
    double  startTime = getTime() - 0.1;
    u64bit  count     = 0;

    theNumberOfMers = u64bitZERO;

    s32bit  timeUntilValid      = 0;
    u64bit  substring           = u64bitZERO;
    u64bit  partialEncoding     = u64bitZERO;
    u32bit  partialEncodingSize = u32bitZERO;

    for (unsigned char ch=nextCharacter(seqFile); ch != 255; ch=nextCharacter(seqFile)) {
      if (ch > 127) {
        timeUntilValid = merSize;
        ch = 0x07;
      } else {
        ch = compressSymbol[ch];

        substring <<= 2;
        substring  |= ch;
        timeUntilValid--;

        if (beVerbose && ((++count & (u64bit)0x3fffff) == u64bitZERO)) {
          fprintf(stderr, " %4ld Mbases -- %8.5f Mb per second\r",
                  count >> 20,
                  count / (getTime() - startTime) / (1000000.0));
          fflush(stderr);
        }

        if (timeUntilValid <= 0) {
          timeUntilValid = 0;

          theNumberOfMers++;

          u64bit I = ((substring >> chkBits) & tblMask) * bitsPerIndex;

          preIncrementDecodedValue(theCounts,
                                   I >> 6,
                                   I & 0x000000000000003f,
                                   bitsPerIndex);
        }
      }

      //  Place the character into the encoded sequence.  If the character was
      //  not a valid base, ch will have the "I'm a crappy character" flag set,
      //  otherwise, ch is a valid base encoding.
      //
      partialEncoding <<= 3;
      partialEncoding  |= ch;

      partialEncodingSize++;

      //  We can fit 21 bases (at 3 bits per base, 63 bits) per 64 bit word
      //  without going to any trouble.
      //
      if (partialEncodingSize == 21) {
        theSeq[theSeqLen++] = partialEncoding;
        theSeq[theSeqLen]   = ~u64bitZERO;

        partialEncoding     = u64bitZERO;
        partialEncodingSize = u32bitZERO;
      }
    }

    //  Push on the last partialEncoding, if one exists.
    //
    if (partialEncodingSize > 0) {
      while (partialEncodingSize != 21) {
        partialEncoding <<= 3;
        partialEncoding  |= 0x07;
        partialEncodingSize++;
      }

      theSeq[theSeqLen++] = partialEncoding;
    }

    if (beVerbose)
      fprintf(stderr, "\n");
  }

The second block of code converts the bucket sizes into an offset to the start of each bucket.

.. code-block:: none

  #include "meryl.H"

  void
  convertCounts(void) {
    double  startTime = getTime() - 0.1;
    u32bit  count     = 0;

    //  Convert the counts into offsets into theMers
    //
    //  We use the trick pointed out by Liliana of setting the offset to the LAST
    //  entry in the bucket, then decrementing when filling.  This also makes
    //  the code below easier -- we don't need to store the size of the bucket
    //  so we can increment the sum after we set the bucket offset.
    //

    u64bit   i = 0;
    u64bit   m = (u64bitONE << tblBits) + 1;
    u64bit   s = 0;
    while (i < m) {

      if (beVerbose && ((++count & (u64bit)0x3fffff) == u64bitZERO)) {
        fprintf(stderr, " %4ld Mbuckets -- %8.5f Mbuckets per second\r",
                count >> 20,
                count / (getTime() - startTime) / (1000000.0));
        fflush(stderr);
      }

      u64bit I = i * bitsPerIndex;

      s = sumDecodedValue(theCounts,
                          I >> 6,
                          I & 0x000000000000003f,
                          bitsPerIndex,
                          s);

      i++;
    }

    if (s != theNumberOfMers) {
      fprintf(stderr, "internal error in stage 2: s != theNumberOfMers.\n");
      exit(1);
    }

    if (beVerbose)
      fprintf(stderr, "\n");
  }

The third block of code makes a second pass through all the kmers in the sequence, adding
suffixes to the bucket indicated by the prefix.

.. code-block:: none

  #include "meryl.H"

  void
  fill(void) {
    double  startTime = getTime() - 0.1;
    u32bit  count     = 0;

    s32bit  timeUntilValid      = 0;
    u64bit  substring           = u64bitZERO;

    //  Stream again, filling theMers
    //

    for (u64bit i=0; i<theSeqLen; i++) {
      u64bit  t = theSeq[i];
      u64bit  I;
      u64bit  J;
      u64bit  ch;

      for (u32bit j=0; j<21; j++) {
        ch   = t;
        ch >>= 60;
        ch  &= 0x07;

  #if 0
        ch   = t & 0x7000000000000000;
        ch >>= 60;
  #endif

        t <<= 3;

        if (ch & 0x04) {
          timeUntilValid = merSize;
        } else {
          substring <<= 2;
          substring  |= ch;
          timeUntilValid--;

          if (beVerbose && ((++count & (u64bit)0x3fffff) == u64bitZERO)) {
            fprintf(stderr, " %4ld Mbases -- %8.5f Mb per second\r",
                    count >> 20,
                    count / (getTime() - startTime) / (1000000.0));
            fflush(stderr);
          }

          if (timeUntilValid <= 0) {
            timeUntilValid = 0;

            I = ((substring >> chkBits) & tblMask) * bitsPerIndex;

            J = chkBits * preDecrementDecodedValue(theCounts,
                                                   I >> 6,
                                                   I & 0x000000000000003f,
                                                   bitsPerIndex);

            setDecodedValue(theMers,
                            J >> 6,
                            J & 0x3f,
                            chkBits,
                            substring & chkMask);
          }
        }
      }
    }
    if (beVerbose)
      fprintf(stderr, "\n");
  }

The fourth block of code sorts the kmers and writes output.  It is using its
own implementation of `heap sort
<https://dl.acm.org/doi/10.1145/512274.512284>_`, likely beacuse the C sort()
implementation is not in place.  Output is written to ``stdout``.

.. code-block:: none

  #include "meryl.H"

  typedef u64bit heapbit;

  void
  adjustHeap(heapbit *M, s64bit i, s64bit n) {
    heapbit     m = M[i];

    s64bit     j = (i << 1) + 1;  //  let j be the left child

    while (j < n) {

      if (j<n-1 && M[j] < M[j+1])
        j++;       //  j is the larger child

      if (m >= M[j])   //  a position for M[i] has been found
        break;

      //  Move larger child up a level
      //
      M[(j-1)/2]      = M[j];

      j = (j << 1) + 1;
    }

    M[(j-1)/2]      = m;
  }


  void
  sortBucket(u64bit b) {
    u64bit  st = getCount(b);
    u64bit  ed = getCount(b+1);
    u64bit  sz = ed - st;

    if (sz == 0)
      return;

    heapbit *a = new heapbit [sz + 1];

    for (u64bit i=st; i<ed; i++) {
      u64bit J = i * chkBits;

      a[i-st] = getDecodedValue(theMers,
                                J >> 6,
                                J & 0x3f,
                                chkBits);
    }

    u64bit  substring;
    char    mer[merSize+1];
    u64bit  count     = u64bitONE;
    u32bit  i;
    u32bit  m;

    //  Sort if there is more than one item
    //
    if (sz > 1) {
      s64bit  i;

      //  Create the heap of lines.
      //
      for (i=(sz-2)/2; i>=0; i--)
        adjustHeap(a, i, sz);

      //  Interchange the new maximum with the element at the end of the tree
      //
      for (i=sz-1; i>0; i--) {
        heapbit m  = a[i];
        a[i]       = a[0];
        a[0]       = m;

        adjustHeap(a, 0, i);
      }
    }

    for (i=1; i<sz; i++) {
      if (a[i] == a[i-1]) {
        count++;
      } else {
        if (count > 100) {
          substring = b << chkBits | a[i-1];
          for (m=0; m<merSize; m++)
            mer[merSize-m-1] = decompressSymbol[(substring >> (2*m)) & 0x03];
          mer[m] = 0;
          printf("%s %lu\n", mer, count);
        }
        count = 1;
      }
    }

    //  Output the last one
    //
    if (count > 100) {
      substring = b << chkBits | a[i-1];
      for (m=0; m<merSize; m++)
        mer[merSize-m-1] = decompressSymbol[(substring >> (2*m)) & 0x03];
      mer[m] = 0;
      printf("%s %lu\n", mer, count);
    }

    delete [] a;
  }


  void
  sort(void) {
    double  startTime = getTime() - 0.1;
    u32bit  count     = 0;
    u64bit  m         = u64bitONE << tblBits;

    for (u64bit i=0; i<m; i++) {
      if (beVerbose && ((++count & (u64bit)0x3fff) == u64bitZERO)) {
        fprintf(stderr, "             %4lu buckets -- %8.5f buckets per second\r",
                count,
                count / (getTime() - startTime));
        fflush(stderr);
      }

      sortBucket(i);
    }
    if (beVerbose)
      fprintf(stderr, "\n");
  }
