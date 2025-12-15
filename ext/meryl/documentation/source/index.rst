Meryl
=====

.. toctree::
   :hidden:

   quick-start
   grammar
   reference
   usage
   history
   faq

`Meryl <http://github.com/marbl/meryl>`_ is a tool for working with sets of
kmers.  A set of kmers, when annotated with an integer value (typically the
number of times the kmer occurs in a set of sequences) and optionally a
label, is termed a database.  Meryl comprises three (modules), one for
generating kmer databases, one for filtering and combining databases, and one
for searching databases.  A simple, but yet to be documented, C++ API to
access kmer databases exists.

Publications
============

Rhie A, Walenz BP, Koren S, Phillippy AM.
`Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies <https://doi.org/10.1186/s13059-020-02134-9>`_. Genome Biology 21, 245 (2020).

Install
=======

The easiest way to get started is to download a
`release <https://github.com/marbl/meryl/releases>`_.
If you encounter any issues, please report them using the
`github issues <http://github.com/marbl/meryl/issues>`_ page.

Alternatively, you can also build the latest unreleased from github:

::

  git clone https://github.com/marbl/meryl.git
  cd meryl/src
  make -j <number of threads>

Learn
=====

NOTE!  This documentation applies to meryl2.  The original meryl is quite
similar, but doesn't support filters or labels, and the specification of
actions is much simpler.

*  :ref:`Quick Start             <quick-start>` - no experience or data required, download and analyze *Escherichia coli* repeats today!
*  :ref:`Usage                   <usage>`       - command line usage cheat sheet
*  :ref:`Grammar                 <grammar>`     - action grammar
*  :ref:`Reference               <reference>`   - all the details
*  :ref:`History                 <history>`     - the history of meryl
*  :ref:`FAQ                     <faq>`         - Frequently asked questions
