.. _reference:

Parameter Reference
====================



A **meryl database** stores kmers, values and labels in a
compressed binary format.

The **value** of a kmer is typically the number of times -- its **count** --
it occurs in some input sequence.

The **label** of a kmer is a 64-bit binary string assigned to each kmer by
the user.  Labels are envisioned to be most useful as a collection of yes/no
or true/false flags, for example, "kmer occurred in file 1".  **Labels are
not fully implemented yet and will appear in a future version.**

The database is stored as 64 independent files, each file storing some subset
of kmers.  Each file can be processed independently, allowing meryl to use up
to 64 threads.

Each kmer is stored by breaking the binary represntation into three pieces: a
file prefix, a block prefix, and a suffix.  A k-mer needs 2*k bits to
represent it.  The file prefix is 6 bits wide, representing one of the 64
possible files in a database.  Inside a file, kmers are stored in blocks,
each kmer in a block will have the same block prefix.  The suffix data for a
block is stored using Elias-Fano encoding (CITE) where each suffix is split
into two pieces.  The first piece is encoded as a unary offset from the last
piece, and the second piece is a N-log2(N)+1 binary value.  At present,
values are stored as plain 32-bit integers, stored immediately after the kmer
data.


