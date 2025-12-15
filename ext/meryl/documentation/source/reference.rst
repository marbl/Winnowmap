.. _reference:

==============
What is Meryl?
==============

Meryl is a tool for creating and working with DNA sequence k-mers.  Each
k-mer is annotated with an integer **value** and a binary **label**.  Sets of
k-mers are stored in binary **databases**, which can be combined using
**actions**.  An action **assigns** a new value and label to each output
k-mer and also **selects** which k-mers are to be output to a new databae or
passed as input to another action.

Initial values are determined by **counting** how many times each k-mer
sequence occurs in some collection of sequences.  Assigned values 
are computed from input values.
Values are usually intepreted as a quantity -- the number of times a k-mer
occurs in a sequence, the number of sequences that contain the k-mer -- and
this is reflected in the operations that can be performed on them; add,
subtract, minimum, maximum, etc.

Initial labels can set to a user-supplied constant counting k-mers; like values, they are
computed from input labels.
Labels can be interpreted as a set of single-bit yes/no flags or multi-bit
values.  For example, 7 bits could be used to indicate which of 7 samples the
k-mer is present in, or a set of 10 bits could be used to `unary encode
<https://en.wikipedia.org/wiki/Unary_coding>`_ the `melting temperature
<https://www.sigmaaldrich.com/US/en/technical-documents/protocol/genomics/pcr/oligos-melting-temp>`_
of the k-mer.

.. include:: reference-kmer.rst

.. include:: reference-acgt-order.rst

Given at least one sequence file, meryl will find the list of k-mers present
in it and ``count`` how many times each one occurs.  This number becomes the
``value`` of the k-mer.  These are stored in a meryl database.  The example
in the sidebar would store:

.. code-block:: none

   ATC 2
   AGA 1
   CTC 1
   TCA 1

Values are typically interpreted as the frequency of the k-mer in some set of
sequences, however, they are simply unsigned 32-bit integers (a maximum value
of 4,294,967,295) and can be used to store any arithmetic data.  Note that a
k-mer **cannot** have a value of zero - this is interpreted as meaning the
k-mer does not exist.

The (optional) label of a k-mer can contain up to 64 bits worth of
non-arithmetic data.  The label can, for example, be used to assign a
``type`` or ``class`` to certain k-mers, or to mark k-mers as coming from a
specific source, etc.  Labels are operated on by the binary logic operations
(AND, OR, XOR, NOT) and can also be shifted to the left or right.

The primary purpose of meryl is to combine multiple k-mer databases into a
single output database by computing new values and labels and filtering
k-mers based on their value, label, base composition and presence or absence
from specific databases.  It does this by passing each k-mer through a tree
of ``actions``.  A leaf node of the tree reads k-mers from input databases
(or by counting k-mers in input sequence files), filtering k-mers via an
action, and emitting k-mers to other nodes or output databases.

=========
Databases
=========

.. include:: reference-database-details.rst

A set of k-mers, each k-mer with a value and a label, is stored in a
**database**.  The database is a directory with 129 binary files in it -- 64
data files, 64 index files and one master index.  This division lets meryl
easily process each of these files independently, making effective use of up
to 64 compute threads.

Databases also store the k-mer size (**-k** option), label size (**-l**
option), and any simple sequence reductions (**--compress** and **--ssr**
options) applied.  It is not possible to combine databases with different
parameters.

Each k-mer is stored at most once per database - thus a k-mer cannot have
multiple values of labels associated with it (though we did envision doing
this at one time).


===============
Counting K-mers
===============

The **count** action reads sequence from any number of input files and counts
the number of times each (canonical) k-mer occurs.  Actions **count-forward**
and **count-reverse** will instead count k-mers as they are oriented in the
input sequence or the reverse-complement sequence, respectively.

Input sequences can be in either FASTA, FASTQ, raw bases, or if compiled with
Canu support, in a Canu seqStore database.  Sequence files can be gzip, bzip2
or xz compressed.

An output database must be supplied to all count actions.  K-mers are both
written to the output database and provided as input to destination actions.

Count actions, unless accompanied by an action that reads input from an
existing database, MUST specify the k-mer size on the command line with the
**-k** option.

Count actions can include a value or label assignment, but cannot include any
selectors.  A value assignment could be used to assign each k-mer a constant value
instead of the count; a label assignment could be used to assign each k-mer
a constant representing the input file.

Counting is resource intense.  Meryl will use memory and threads up to a
limit supplied by: the operating system (usually physical memory and the
number of CPUs), a grid manager (such as Slurm, PBS or SGE) or a command line
option (**-m** and **-t**).

Two algorithms are used for counting k-mers.  The algorithm that is expected
to use the least memory is used.  The choice depends on the size of the input
sequences and the k-mer size.

Counting Small k-mers (k < 17)
------------------------------

.. warning::
  does count really only use one thread here?

.. warning::
  is this method used even for small amount of input sequence?

.. sidebar:: Small k-mer Counting Implementation Details

  Each integer counter is initially a 16-bit value.  Once any count exceeds
  2\ :sup:`16` = 65,535 another bit is added to all value, resulting in
  17-bit values for every k-mer.  Once any count then exceeds 2\ :sup:`17` =
  131,072, another bit is added, and so on.  Thus, memory usage is 512 MB *
  log\ :sub:`2` maximum_count_value

For k at most 16, meryl counts k-mers directly, that is, by associating an
integer count with each possible k-mer.  This has the benefit of being simple
and uses a constant amount of memory regardless of the size of the input, but
quickly exhausts memory for even moderate k-mer sizes.

There are 4\ :sup:`k` k-mers of size k; for k=16, there are 4,294,967,296
possible k-mers.  Counting 16-mers with this method will use at least 8
GB of memory, independent of input size: counting 16-mers in an E.coli genome
will use 8 GB of memory, despite there being only 5 million or so k-mers.  Further,
memory usage can increase depending on the maximum count value.

This method uses only a single thread to read the input sequence and
increment counters in the array, but multiple threads can be used to generate
the output database.

Counting Large k-mers (k > 15)
------------------------------

.. sidebar:: Large k-mer Counting Implementation Details

  Each k-mer is split into a prefix and a suffix.  The prefix is used to
  select a list to which the suffix is added.  When the (approximate) size of
  all lists exceeds a user-supplied threshold, each list is sorted, the
  suffixes are counted, and this subset of counted k-mers is output to an
  intermediate database.  After all k-mers are processed, the intermediate
  databases are merged into one.

For k larger than 15, or for small amounts of input sequence, meryl counts
k-mers by first converting the sequence to a list of k-mers, duplicates
included, then sorts the list, then scans the list to count the number of
times each k-mer is present.

If all k-mers in an input sequence do not fit in memory, a partial result is
written to disk.  After all input sequences have been processed, the partial
results are combined into a single output database.  In practice, this method
requires several additional gigabytes of memory to minimize the overhead of
writing and merging partial results.

This method can use multiple threads for every stage.

=======
Actions
=======

Meryl processing is built around **actions**.  An action loads a k-mer from
one or multiple databases (or, for counting actions, computes the k-mer from
a sequence file) computes value and label for it (based on the input kmer
values and labels), decides if it should be output or discarded (e.g., "if
the new value is more than 100, output the k-mer"), and saves it to an output
database or text file or displays it on the screen or passes it to
another action for further processing.

An action is specified as an alias (listed below) or by explicitly stating
all parameters.  The parameters describe:

 * how to compute the value of the k-mer
 * how to compute the label of the k-mer
 * conditions when a k-mer should be output or discarded:
    - based on which input databases it came from
    - based on the input and/or output values of the k-mer
    - based on the input and/or output labels of the k-mer
    - based on the sequence of the k-mer
 * what to do with output k-mers
    - output them to a new database
    - print them to ASCII output files
    - display them on the terminal
    - pass them to other actions

.. note::
  K-mers are read "in order" from the inputs.  If an input does not contain
  the "next" k-mer, it does not participate in the action processing.  For example,
  suppose we have three input databases with the following 4-mers and their counts:

  .. code-block:: none
    :caption: Sample databases.
    :linenos:

    input-1  input-2  input-3
    AAAA/1   AAAA/2   AAAA/3
    AAAC/1   CAAT/2   CCCC/3
    CAAT/1            GGGG/3
    GGGG/1

  A ``union-sum`` action with these three databases as input will output:

  .. code-block:: none
    :caption: Sample output from union-sum action.
    :linenos:

    AAAA/6 (using the k-mer from input-1, input-2 and input-3)
    AAAC/1 (... from input-1)
    CAAT/3 (... from input-1 and input-2)
    CCCC/3 (... from input-3)
    GGGG/4 (... from input-1 and input-3)

An **assignment** computes the output value (label) for each k-mer from among
the input values (labels).  At most one assignment can be supplied for the
value and one for the label.

A **selector** decides if the k-mer should be output or discarded.  Selectors can
use input values (labels), the computed output value (label), the base composition
of the k-mer and which specific inputs the k-mer was present in.
Any number of selectors can be supplied, linked with **and**, **or** and
**not** operators.  See SELECTORS.

Though it is possible to specify all those choices explicitly, **aliases** are
provided for most common operations.

Aliases exist to support common operations.  An alias sets the ``value``,
``label`` and ``input`` options and so these are not allowed to be used with
aliases.  Examples of aliases and their explicit configuration:

.. warning::
  OBSOLETE!

.. table:: Action Aliases
  :widths: 25 30 45

  +--------------------+--------------------+----------------------------------------------+
  | Alias              | Output k-mer if... | Sets value to the...                         |
  +====================+====================+==============================================+
  | union              | ...k-mer is in any | ...number of databases the k-mer is in.      |
  +--------------------+ input database.    +----------------------------------------------+
  | union-min          |                    | ...smallest input value.                     |
  +--------------------+                    +----------------------------------------------+
  | union-max          |                    | ...largest input value.                      |
  +--------------------+                    +----------------------------------------------+
  | union-sum          |                    | ...sum of the input values.                  |
  +--------------------+--------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+
  | intersect          | ...k-mer is in all | ...value of the k-mer in the first database. |
  +--------------------+ input databases.   +----------------------------------------------+
  | intersect-min      |                    | ...smallest input value.                     |
  +--------------------+                    +----------------------------------------------+
  | intersect-max      |                    | ...largest input value.                      |
  +--------------------+                    +----------------------------------------------+
  | intersect-sum      |                    | ...sum of the input values.                  |
  +--------------------+--------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+
  | subtract           | ...k-mer is in the | ...value of the k-mer in the first database  |
  |                    | first database.    | minus all other values.                      |
  +--------------------+                    +----------------------------------------------+
  | difference         |                    | ...value of the k-mer in the first database. |
  +--------------------+--------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+
  | less-than X        | ...k-mer is in the | ...value of the k-mer.                       |
  +--------------------+ first and only     |                                              |
  | greater-than X     | database and the   |                                              |
  +--------------------+ value meets the    |                                              |
  | at-least X         | speficied          |                                              |
  +--------------------+ condition.         |                                              |
  | at-most X          |                    |                                              |
  +--------------------+                    |                                              |
  | equal-to X         |                    |                                              |
  +--------------------+                    |                                              |
  | not-equal-to X     |                    |                                              |
  +--------------------+--------------------+----------------------------------------------+
  +--------------------+--------------------+----------------------------------------------+
  | increase X         | ...k-mer is in the | ...value of the k-mer modified by            |
  +--------------------+ first and only     | the specified operation.                     |
  | decrease X         | database.          |                                              |
  +--------------------+                    | (divide-round rounds 0 up to 1)              |
  | multiple X         |                    |                                              |
  +--------------------+                    |                                              |
  | divide X           |                    |                                              |
  +--------------------+                    |                                              |
  | divide-round X     |                    |                                              |
  +--------------------+                    |                                              |
  | modulo X           |                    |                                              |
  +--------------------+--------------------+----------------------------------------------+

.. warning::
  OBSOLETE!  This table has not been verified!

.. table:: Action Aliases
  :widths: 19 19 19 16 14 13

  +----------------+---------------------------------------------------------------------------------+
  |                |                                    Action                                       |
  | Alias          +------------------------------------+--------------------------------------------+
  |                + Assignment                         | Selector                                   |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | union          | value=count       | label=or       | input:any    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | union-min      | value=min         | label=min-value| input:any    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | union-max      | value=max         | label=max-value| input:any    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | union-sum      | value=sum         | label=or       | input:any    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | intersect      | value=first       | label=and      | input:all    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | intersect-min  | value=min         | label=min-value| input:all    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | intersect-max  | value=max         | label=max-value| input:all    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | intersext-sum  | value=sum         | label=and      | input:all    | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | subtract       | value=sub         | label=first    | input:first  | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | difference     | value=sub         | label=first    | input:first  | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | less-than X    | value=first       | label=first    | input:only   | value:<X     | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | greater-than X | value=first       | label=first    | input:only   | value:>X     | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | at-least X     | value=first       | label=first    | input:only   | value:>=X    | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | at-most X      | value=first       | label=first    | input:only   | value:<=X    | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | equal-to X     | value=first       | label=first    | input:only   | value:==X    | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | not-equal-to X | value=first       | label=first    | input:only   | value:!=X    | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | increase X     | value=\@1+X       | label=first    | input:only   | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | decrease X     | value=\@1-X       | label=first    | input:only   | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | multiply X     | value=\@1*X       | label=first    | input:only   | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | divide X       | value=\@1/X       | label=first    | input:only   | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | divide-round X | value=\@1/X [#a]_ | label=first    | input:only   | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+
  | modulo X       | value=\@1%X       | label=first    | input:only   | value:       | label:       |
  +----------------+-------------------+----------------+--------------+--------------+--------------+

.. [#a] The ``divide-round`` alias rounds values of 0 up to 1.

A full action is:

.. code-block:: none
  :caption: Fully general action template.
  :linenos:

  [ action-name
      output:database=<database.meryl>            # to a binary database
      output:list=<files.##.mers>                 # to ascii files
      output:show                                 # to stdout
      output:pipe=<label>                         # to a 'named pipe' action input
      output:histogram=<hist-file>                # compute a histogram
      output:statistics=<stats-file>              # compute statistics
      assign:value=<rule-to-create-output-value>
      assign:label=<rule-to-create-output-label>
      [not] select:value:<rule-to-select-k-mer-for-output> [or | and]  # if none specified
      [not] select:label:<rule-to-select-k-mer-for-output> [or | and]  # and is assumed
      [not] select:bases:<rule-to-select-k-mer-for-output> [or | and]
      [not] select:input:<rule-to-select-k-mer-for-output>
      [input:database] input-1                           # from a database or
      [input:list]     input-2                           # list file
      [input:pipe]     input-3                           # 'named pipe' label
      [input:action]   [ action-name ... ]
      ...
  ]

The 'action-name' is used to either set default values from an alias, or to
define a new alias (when no inputs are supplied, see ALIASES).

Notice that each parameter is described by a parameter-class ('output') and a
parameter-name ('database').  There are four parameter-classes: 'output',
'assign', 'select' and 'input'.  Parameter-class 'input' is optional as meryl
can query the input source to determine what type it is.  Parameter-classes
and parameter-names may be shortened to any prefix of the class or name, all
the way to a single letter.  In some cases, a different word can be
substituted, however, these may not be shortened and must be used as shown:

.. table:: Parameter class abbreviations and aliases
  :widths: 34 33 33

  +-----------------+---------------------------+------------------------------+
  | Parameter Class | Suggested Abbreviations   | Aliases                      |
  +-----------------+---------------------------+------------------------------+
  | output          | o, out,                   |                              |
  +-----------------+---------------------------+------------------------------+
  | assign          | a                         | set                          |
  +-----------------+---------------------------+------------------------------+
  | select          | s, sel                    | get                          |
  +-----------------+---------------------------+------------------------------+
  | input           | i, inp, in                |                              |
  +-----------------+---------------------------+------------------------------+

.. table:: Parameter name abbreviations and aliases
  :widths: 34 33 33

  +-----------------+---------------------------+------------------------------+
  | Parameter Name  | Suggested Abbreviations   | Aliases                      |
  +-----------------+---------------------------+------------------------------+
  | database        | d, dat                    | db                           |
  +-----------------+---------------------------+------------------------------+
  | list            | l, lis                    | t, txt, text                 |
  +-----------------+---------------------------+------------------------------+
  | show            | s, sho                    | display, dis, print, stdout  |
  +-----------------+---------------------------+------------------------------+
  | pipe            | p, pip                    |                              |
  +-----------------+---------------------------+------------------------------+
  | histogram       | h, hist, histo            |                              |
  +-----------------+---------------------------+------------------------------+
  | statistics      | st, stats                 | NOTE: 's' is NOT an abbrev.  |
  +-----------------+---------------------------+------------------------------+
  | action          | a, act                    |                              |
  +-----------------+---------------------------+------------------------------+
  | value           | v, val                    |                              |
  +-----------------+---------------------------+------------------------------+
  | label           | l, lab                    |                              |
  +-----------------+---------------------------+------------------------------+
  | bases           | b, bas                    | acgt, bp                     |
  +-----------------+---------------------------+------------------------------+
  | input           | i, inp, in                |                              |
  +-----------------+---------------------------+------------------------------+


.. note::

  Selectors are 'positive-sense', meaning they initially assume every k-mer
  is to be discarded, then 'select' only those k-mers that (positively) match
  some set of conditions.  A 'negative-sense' scheme would initially assume
  all k-mers are selected, and then de-select those that match the set of
  conditions.  One can always convert selectors between senses, but the two
  styles cannot be mixed.  Here's why.

  We'll use 'select' for a positive-sense selector and 'reject' for a
  negative-sense selector.  These following two sets are equivalent:

  .. code-block:: none
    :linenos:

    select:value:>=3 and select:value:<=10
    reject:value:<3  or  reject:value:>10

  either will unambiguously discard items with value less than three, accept
  items with value between 3 and 10, and discard items with value more
  than 10.

  De Morgan rests comfortably; these are all equivalent:

  .. code-block:: none
    :linenos:

    select:value:>=3 and select:value:<=10
    not (not (select:value:>=3 and select:value:<=10))
    not (not select:value:>=3) or (not select:value:<=10)
    not (select:value:<3) or (select:value:>10)

  But something odd happens if we replace a single 'select' term with a
  seemingly equivalent 'reject' term.  Replcing the first term
  (`select:value:>=3`) in the positive-sense selector (line 1 in both examples
  above) with an 'equivalent' reject term (`reject:value:<3`) results in a
  contradiction:

  .. code-block:: none
    :linenos:

    reject:value:<3 or select:value:<=10

  For values less than three, we're told to both 'reject' and 'select' the k-mer!
  We could probably devise rules of precedence to resolve this problem, but
  that would quickly lead to confusion.

  A more sinister contradiction is that the 'reject' term implicitly accepts
  items with value 3 or larger unless some other term says to discard it; the
  'select' term does the opposite - it implicitly discards items with value
  more than ten.  Values greater than ten are both (implicitly) discarded and
  accepted!  discarding and implicitly accepting them.

  For simplicity, we chose to implement only one sense of select term.  It
  might be possible to later implement the other sense, or even add
  'exceptions' of the form `select:value>=3 except reject:label=2`.




Square brackets MUST surround every action (exception: the first action in a
command tree can omit the brackets).

``output:database`` (``o:d``) is optional, but may occur at most once per
action.  If present, the k-mers generated by this action will
be written to the specified meryl database.  If an existing database is
supplied, it will be overwritten.

``output:list`` (``o:l``) is optional.  If present, the k-mers generated by this action will be
written to ASCII file(s) in the format ``<k-mer><tab><value><tab><label>``,
one k-mer per line.  The k-mers will be in sorted order: A, C, T, G.  If the
file name includes the string ``##``, the data will be written to 64 files,
in parallel, using up to 64 threads.  Appending suffix ``.gz``, ``.bz2`` or
``.xz`` will cause the output file to be compressed.

``output:show`` (``o:s``) is optional.  It behaves like ``output:list``,
except the k-mers are written to the ``stdout``, the terminal, unless
redirected.

``output:pipe`` (``o:p``) is optional.  If present, the k-mers generated by
this acction will be supplied to other meryl actions that read input from the
same pipe name.

.. note::

  ``listACGT`` is the same as ``list``, but modifies the ordering of k-mers
  from ``A < C < T < G`` to ``A < C < G < T`` when forming canonical k-mers.
  While this generates correct canonical k-mers, the output k-mers are not
  sorted.

  Consider 3-mers from string ``GGAGAGCT``:

  .. table:: ACTG order vs ACGT order
    :widths: 20 20 20 20 20

    +------------+---------+---------+---------------------------+
    |            |         |         |     canonical k-mer in    |
    |            |         |         +-------------+-------------+
    | ``GGAGCT`` | forward | reverse |  ACTG order |  ACGT order |
    +------------+---------+---------+-------------+-------------+
    | ``GGA...`` | ``GGA`` | ``TCC`` |   ``TCC``   |   ``GGA``   |
    +------------+---------+---------+-------------+-------------+
    | ``.GAG..`` | ``GAG`` | ``CTC`` |   ``CTC``   |   ``CTC``   |
    +------------+---------+---------+-------------+-------------+
    | ``..AGC.`` | ``AGC`` | ``GCT`` |   ``AGC``   |   ``AGC``   |
    +------------+---------+---------+-------------+-------------+
    | ``...GCT`` | ``GCT`` | ``AGC`` |   ``AGC``   |   ``AGC``   |
    +------------+---------+---------+-------------+-------------+

  When meryl builds the datase, it uses the ``A < C < T < G`` order.  These
  k-mers will be stored in the database in order: ``AGC``, ``CTC``, ``TCC``,
  ``GCT``.  But when output using ``listACGT``, the k-mers will be reported as
  ``AGC``, ``CTC``, ``GGA``, ``GCT``.  Notice that because of the change in
  canonical k-mer from ``TCC`` to ``GGA`` the last k-mer is not in sorted order.

``assign:value=`` (``a:v=``) and ``assign:label=`` (``a:l=``) describe how to
combine the input values and labels into a single output value and label.

``select:value:`` (``s:v``), ``select:label:`` (``s:l``), ``select:bases:``
(``s:b``) and ``select:input:`` (``s:i``) describe the conditions required
for a k-mer to output.  Any number of these may be supplied.

An ``input`` is either a meryl database, a list of k-mers, or another meryl
action.

Some actions require exactly one input, others require more than
one - this is specified in the ``select:input:`` rule.

===========
Assignments
===========

.. warning::
  HOW IS THIS IMPLEMENTED?

  When value:#c, value:first, value:min or value:max are used, the label
  operation acts ONLY on the k-mers matching the value selection.  For example,
  if value:min finds value=5 is the minimu, label=or would combine the labels
  of all k-mers with value=5.  Contrast this with value:add (which would set the
  output value to the sum of the k-mer values in all databases) and label:and
  (which would set each bit in the output label to true if the corresponding
  bit is true in all inputs).

  Likewise for label:#c, label:first, label:minweight and label:maxweight.  For
  example, when label:#c is used, value:add would sum the values of all labels
  that are the same as constant c.



Value Assignment
----------------

A **value assignment** computes the output value of the k-mer
based on the input values and possibly a single integer constant.

.. note::
  The optional parameter ``(#X)`` means to also include constant ``X`` in the
  computation.

.. note::
  Constants can be decimal integers (``123`` or ``123d``), hexadecimal (``abch``),
  octal (``147o``) or binary (``0101010b``).  SI suffixes can be used on plain
  decimal integers (``123k`` == 123,000; ``1mi`` == 1,048,576).  For example,
  ``value=add#10`` would set the output value to the sum of the input values
  plus ten; ``value=min#10`` would set the output value to the smallest input
  value or 10 if all input values are larger than 10.

.. warning::
  How to form complex expressions?

.. warning::
  Things like value=@1-@2 are NOT supported.  Even the potentially useful
  value=@1 isn't supported (though it is listed below).

.. warning::
  value=selected isn't implemented.

.. table:: Value Assignments
  :widths: 20 80

  +--------------------+-------------------------------------------------+
  | Assignment         | Set value to ...                                |
  +====================+=================================================+
  | value=#X           | ...constant X.                                  |
  +--------------------+-------------------------------------------------+
  | value=@X           | ...that of the k-mer in the Xth input           |
  +--------------------+-------------------------------------------------+
  | value=first        | ...that of the k-mer in the first input.        |
  +--------------------+-------------------------------------------------+
  | value=selected     | ...that of the k-mer selected by the label=     |
  |                    | selector.  When multiple k-mers are selected,   |
  |                    | the value of the first is used.                 |
  +--------------------+-------------------------------------------------+
  | value=min(#X)      | ...the minimum of all input values.             |
  +--------------------+-------------------------------------------------+
  | value=max(#X)      | ...the maximum of all input values.             |
  +--------------------+-------------------------------------------------+
  | value=add(#X)      | ...the sum of all input values.                 |
  +--------------------+-------------------------------------------------+
  | value=sum(#X)      | ...the sum of all input values.                 |
  +--------------------+-------------------------------------------------+
  | value=sub(#X)      | ...the value of the k-mer in the first input    |
  |                    | minus all other values.                         |
  +--------------------+-------------------------------------------------+
  | value=dif(#X)      | ...the value of the k-mer in the first input    |
  |                    | minus all other values.                         |
  +--------------------+-------------------------------------------------+
  | value=mul(#X)      | ...the product of all input values.             |
  +--------------------+-------------------------------------------------+
  | value=div(#X)      | ...the value of the k-mer in the first input    |
  |                    | divided by all other values.                    |
  +--------------------+-------------------------------------------------+
  | value=divzero(#X)  | ...the value of the k-mer in the first input    |
  |                    | divided by all other values, rounding zero up   |
  |                    | to one.                                         |
  +--------------------+-------------------------------------------------+
  | value=mod(#X)      | ...the remainder after the value of the k-mer in|
  |                    | the first input is divided by all other values. |
  +--------------------+-------------------------------------------------+
  | value=rem(#X)      | ...the remainder after the value of the k-mer in|
  |                    | the first input is divided by all other values. |
  +--------------------+-------------------------------------------------+
  | value=count        | ...the number of inputs the k-mer is present in.|
  +--------------------+-------------------------------------------------+

Label Assignment
----------------

A **label assignment** computes the output label of the k-mer
based on the input label and possibly a single 64-bit constant.

.. warning::
  How to form complex expressions?

.. table:: Value Assignments
  :widths: 20 80

  +------------------------+-------------------------------------------------+
  | Assignment             | Set label to ...                                |
  +========================+=================================================+
  | label=#X               | ...constant X.                                  |
  +------------------------+-------------------------------------------------+
  | label=@X               | ...that of the k-mer in the Xth input           |
  +------------------------+-------------------------------------------------+
  | label=first            | ...that of the k-mer in the first input.        |
  +------------------------+-------------------------------------------------+
  | label=selected         | ...that of the k-mer selected by the value=     |
  |                        | selector.  When multiple k-mers are selected,   |
  |                        | the label of the first is used.                 |
  +------------------------+-------------------------------------------------+
  | label=min(#X)          | ...the minimum of all input labels.             |
  +------------------------+-------------------------------------------------+
  | label=max(#X)          | ...the maximum of all input labels.             |
  +------------------------+-------------------------------------------------+
  | label=and(#X)          | ...the bitwise AND of all input labels.         |
  +------------------------+-------------------------------------------------+
  | label=or(#X)           | ...the bitwise OR of all input labels.          |
  +------------------------+-------------------------------------------------+
  | label=xor(#X)          | ...the bitwise XOR of all input labels.         |
  +------------------------+-------------------------------------------------+
  | label=difference(#X)   | ...that of the k-mer in the first input,        |
  |                        | with all bits set in other inputs turned off.   |
  +------------------------+-------------------------------------------------+
  | label=lightest(#X)     | ...the label with the fewest bit set.           |
  +------------------------+-------------------------------------------------+
  | label=heaviest(#X)     | ...the label with the most bits set.            |
  +------------------------+-------------------------------------------------+
  | label=invert(#X)       | ...the bitwise invert of the first input.       |
  +------------------------+-------------------------------------------------+
  | label=shift-left(#X)   | ...the first input shifted left by X places.    |
  +------------------------+-------------------------------------------------+
  | label=shift-right(#X)  | ...the first input shifted right by X places.   |
  +------------------------+-------------------------------------------------+
  | label=rotate-left(#X)  | ...the first input rotated left by X places.    |
  +------------------------+-------------------------------------------------+
  | label=rotate-right(#X) | ...the first input rotated right by X places.   |
  +------------------------+-------------------------------------------------+

=========
Selectors
=========

A selector decide if the k-mer should be output.  They can use the values and
labels of the input k-mers, the computed value and label of the k-mer to be
output, the number and location of inputs that supplied an input k-mer, and
the base composition of the k-mer.  A single select term tests one condition,
e.g., ``value:>3``, and multiple terms are connected together in a
sum-of-products form (e.g., 'and' has higher precedence than 'or'):

.. code-block:: none
  :caption: Sum-of-Products filters.
  :linenos:

  value:@1>=20 or value:@2>=20 or value:>30 and input:#2

will output a k-mer if it has a value of at least 20 in either input
database, or the output value is more than 30 and the k-mer occurs in both
inputs, or both conditions are met.

The 'not' keyword has highest precedence and can be used to invert the sense
of the next term, and only the next term.  While this seems restrictive,
`De Morgan's laws <https://en.wikipedia.org/wiki/De_Morgan's_laws>`_ are useful:

.. code-block:: none
  :caption: De Morgan's laws
  :linenos:

    not (A and B) = (not A) or  (not B)
    not (A or  B) = (not A) and (not B)

Do not confuse selectors (which use a ``:`` - 'select:value:', 'select:label:', 'select:input:', 'select:bases:') with
assignments (which use an ``=`` - 'assign:value=' and 'assign:label=').

Value Selectors
---------------

A value selector discards the k-mer from output if the input or output values
are undesired.  When the selector is TRUE the k-mer is output.  The syntax and
options are similar to **label selectors**, but value selectors are typically
integer functions.

.. code-block:: none

  value:<ARG1><OP><ARG2>

ARG1 and ARG2 can be an input file (``@3``), a constant (``#4`` or ``4``), a
special function (ARG2 only) or empty (ARG1 only).

.. table::
  :widths: 20 10 20 50

  +------------------+------+-------------------+--------------------------------------------------------+
  | ARG1             |  OP  | ARG2              | Meaning                                                |
  +==================+======+===================+========================================================+
  | ``@n``           |      | ``@n``            | Use the value from the k-mer in the ``n``\th input.    |
  +------------------+------+-------------------+--------------------------------------------------------+
  | ``#n`` or ``n``  |      | ``#n`` or ``n``   | Use the constant ``n``.                                |
  +------------------+------+-------------------+--------------------------------------------------------+
  | <not-present>    |      |                   | Use the value of the selected output k-mer.            |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |      | ``distinct=f``    | Use the value such that ``f`` fraction of the distinct |
  |                  |      |                   | k-mers have at most this value.  That is,              |
  |                  |      |                   | ``value:ge:distinct=0.9`` will output the 10% most     |
  |                  |      |                   | repetitive k-mers in the database.                     |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |      | ``word-freq=f``   | Compute the word-frequency of a k-mer as its value     |
  |                  |      |                   | divided by the total number of k-mers in the database. |
  |                  |      |                   |                                                        |
  |                  |      |                   |                                                        |
  |                  |      |                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |      | ``threshold=n``   | Use the constant ``n``.                                |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``==``|                   | TRUE if ARG1 equals ARG2.                              |
  |                  |``=`` |                   |                                                        |
  |                  |``eq``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``!=``|                   | TRUE if ARG1 does not equal ARG2.                      |
  |                  |``<>``|                   |                                                        |
  |                  |``ne``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``<=``|                   | TRUE if ARG1 is less than or equal to ARG2.            |
  |                  |``le``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``>=``|                   | TRUE if ARG1 is greater than or equal to ARG2.         |
  |                  |``ge``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``<`` |                   | TRUE if ARG1 is less than ARG2.                        |
  |                  |``lt``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``>`` |                   | TRUE if ARG1 is greater than ARG2.                     |
  |                  |``gt``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+

Note that ``@1`` is not necessarily the first file supplied to the action.  If
the k-mer occurs only in the last file, ``@1`` will be the value of the k-mer in
that file.

Examples:

.. table::
  :widths: 25 75

  +--------------------+------------------------------------------------------------------------------+
  |                    | TRUE if ...                                                                  |
  +--------------------+------------------------------------------------------------------------------+
  | ``value:>5``       | ...the output value is more than 5.                                          |
  +--------------------+------------------------------------------------------------------------------+
  | ``value:@2<=#52o`` | ...the value of the second input is at most 52\ :sub:`8` (or 42\ :sub:`10`). |
  +--------------------+------------------------------------------------------------------------------+
  | ``value:4>@2``     | ...4 is larger than the value of the second input.                           |
  +--------------------+------------------------------------------------------------------------------+
  | ``value:@1>@2``    | ...the value of the first input is more than the second input.               |
  +--------------------+------------------------------------------------------------------------------+



Label Selectors
---------------

A label selector discards the k-mer from output if the input or output labels
are undesired.  When the selector is TRUE the k-mer is output.  The syntax and
options are similar to **value selectors**, but label selectors are typically
binary functions.

.. code-block:: none

  label:<ARG1><OP><ARG2>

ARG1 and ARG2 can be an input file (``@3``), a constant (``#0100b`` or ``4h``), a
special function (ARG2 only) or empty (ARG1 only).

.. table::
  :widths: 20 10 20 50

  +------------------+------+-------------------+--------------------------------------------------------+
  | ARG1             |  OP  | ARG2              | Meaning                                                |
  +==================+======+===================+========================================================+
  | ``@n``           |      | ``@n``            | Use the label from the k-mer in the ``n``\th input.    |
  +------------------+------+-------------------+--------------------------------------------------------+
  | ``#n`` or ``n``  |      | ``#n`` or ``n``   | Use the constant ``n``.                                |
  +------------------+------+-------------------+--------------------------------------------------------+
  | <not-present>    |      |                   | Use the label of the selected output k-mer.            |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |      | ``distinct=f``    | Use the label such that ``f`` fraction of the distinct |
  |                  |      |                   | k-mers have at most this label.                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |      | ``word-freq=f``   | (same, but for total k-mers?)                          |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |      | ``threshold=n``   | Use the constant ``n``.                                |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``==``|                   | TRUE if ARG1 equals ARG2.                              |
  |                  |``=`` |                   |                                                        |
  |                  |``eq``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``!=``|                   | TRUE if ARG1 does not equal ARG2.                      |
  |                  |``<>``|                   |                                                        |
  |                  |``ne``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``<=``|                   | TRUE if ARG1 is less than or equal to ARG2.            |
  |                  |``le``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``>=``|                   | TRUE if ARG1 is greater than or equal to ARG2.         |
  |                  |``ge``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``<`` |                   | TRUE if ARG1 is less than ARG2.                        |
  |                  |``lt``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``>`` |                   | TRUE if ARG1 is greater than ARG2.                     |
  |                  |``gt``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+

Note that ``@1`` is not necessarily the first file supplied to the action.  If
the k-mer occurs only in the last file, ``@1`` will be the label of the k-mer in
that file.

.. table:: Proposed Flters
  :widths: 20 80

  +--------------------+------------------------------------------------------------+
  | Selector           | TRUE if...                                                 |
  +====================+============================================================+
  | label:all#c        | ...all bits set in c are also set in the label             |
  |                    |                                                            |
  | label:and#c        | equivalent to ``l & c == c`` or ``l | c == l``             |
  |                    |                                                            |
  |                    | equivalent to ``~l & c == 0`` (not expressible in meryl)   |
  +--------------------+------------------------------------------------------------+
  | label:any#c        | ...any bits set in c are also set in the label             |
  |                    |                                                            |
  | label:or#c         | equivalent to ``l & c != 0``                               |
  +--------------------+------------------------------------------------------------+
  | label:none#c       | ...no bits set in c are set in the label                   |
  |                    |                                                            |
  | label:not#c        | equivalent to ``l & c == 0``                               |
  +--------------------+------------------------------------------------------------+
  | label:only#c       | ...only the bits set in c are set in the label             |
  |                    |                                                            |
  | label:xor#c ??     | equivalent to ``l | c == c`` or ``l && c == l``            |
  |                    |                                                            |
  |                    | equivalent to ``none#~c``                                  |
  +--------------------+------------------------------------------------------------+
  | label:and#c=d      | ...not expressible in meryl                                |
  |                    |                                                            |
  +--------------------+------------------------------------------------------------+
  | label:or#c=d       | ...not expressible in meryl                                |
  +--------------------+------------------------------------------------------------+

Examples:

We want to find k-mers in an that are in none, one or two different read
sets.  We'll assign distinct indicator bits to each input, ``union`` everything
together, then pick out k-mers that have the `assembly` indicator set.

.. code-block:: none
  :caption: Finding unsupported k-mers, version 1.
  :linenos:

  meryl \
    union \
      output=labeled.meryl \
      value=@3 \
      label=or \
      label:and#0b100 \
      [ label=0b001 db1.meryl ] \
      [ label=0b010 db2.meryl ] \
      [ label=0b100 asm.meryl ]

The result will by k-mers with the value from the assembly and labeled with:

.. table::
  :widths: 20 80

  +-------+---------------------------------------------------+
  | 0b0.. | (label never occurs)                              |
  +-------+---------------------------------------------------+
  | 0b100 | appears only in asm                               |
  +-------+---------------------------------------------------+
  | 0b101 | appears in asm and db1                            |
  +-------+---------------------------------------------------+
  | 0b110 | appears in asm and db2                            |
  +-------+---------------------------------------------------+
  | 0b111 | appears in asm and both db1 and db1               |
  +-------+---------------------------------------------------+

An alternate method is to first intersect the read k-mers with the assembly,
then merge those two sets:

.. code-block:: none
  :caption: Finding unsupported k-mers, version 2.
  :linenos:

  meryl \
    union \
    output=labeled.meryl \
      value=sum \
      label=or \
      [ intersect value=@2 label=0b01 asm.meryl db1.meryl ] \
      [ intersect value=@2 label=0b10 asm.meryl db2.meryl ]

The result is slightly different.  We no longer output k-mers that exist only
in the assembly; the value of k-mers will by the sum of the values in the
read databases, and the label will be:

.. table::
  :widths: 20 80

  +-------+---------------------------------------------------+
  |  0b00 | (label never occurs)                              |
  +-------+---------------------------------------------------+
  |  0b01 | appears in the assembly and db1                   |
  +-------+---------------------------------------------------+
  |  0b10 | appears in the assembly and db2                   |
  +-------+---------------------------------------------------+
  |  0b11 | appears in the assembly and both db1 and db1      |
  +-------+---------------------------------------------------+

Find k-mers with a value of at least ten that exist in two or more databases,
report which and how many databases contained the k-mer.

.. code-block:: none
  :caption: Finding supported k-mers from multiple databases.
  :linenos:

  meryl \
    display \
    value=count \
    label=or \
    inputs:2-4 \
    [ value:>=10 label:0001b a.meryl ] \
    [ value:>=10 label:0010b b.meryl ] \
    [ value:>=10 label:0100b c.meryl ] \
    [ value:>=10 label:1000b d.meryl ]

Output a k-mer if it exists in at least three databases with count greater
than 100, but output the minimum count the k-mer has in any input database.

.. code-block:: none
  :caption: Minimum value of well-supported k-mers.
  :linenos:

  meryl \
  intersect \
    value:@2 \
    [ input:3-5 \
      [ value:>100 a.meryl ] \
      [ value:>100 b.meryl ] \
      [ value:>100 c.meryl ] \
      [ value:>100 d.meryl ] \
      [ value:>100 e.meryl ] \
    ]
    [ union-min \
      a.meryl \
      b.meryl \
      c.meryl \
      d.meryl \
      e.meryl \
    ] \

The first sub-action generates a list of k-mers that are well-supported in at
least three inputs.  Its sub-actions return lists of k-mers with value
greater than 100.  The second sub-action returns a list of all k-mers with
their actual value.  Finally, ``intersect`` returns a list of k-mers that
are both "well-supported in at least three inputs" and "in any input" and
sets the output value to whatever was in the second input.

A Generalized Label Selector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A general label selector can be devised by supplying a function that converts
each bit in the label to some testable output bit, then testing those output
bits.

An example will follow the tables.

The four functions are:

.. table:: 

  +-----------+------+--------------+------------------------+
  | function  | code | output value |                        |
  +===========+======+==============+========================+
  | zero(bit) | 0    | 0            | always false           |
  +-----------+------+--------------+------------------------+
  | one(bit)  | 1    | 1            | always true            |
  +-----------+------+--------------+------------------------+
  | pass(bit) | \+   | bit          | true if label is set   |
  +-----------+------+--------------+------------------------+
  | flip(bit) | \-   | !bit         | true if label is unset |
  +-----------+------+--------------+------------------------+

And the five tests are:

.. table:: 

  +-------------------+------+---------------------------------------+
  | test              | code |                                       |
  +===================+======+=======================================+
  | all-must-be-true  | T    | All 'T' bits must be true.            |
  +-------------------+------+---------------------------------------+
  | any-must-be-true  | t    | At least one 't' bit must be true.    |
  +-------------------+------+---------------------------------------+
  | any-must-be-false | f    | At least one 'f' bit must be false.   |
  +-------------------+------+---------------------------------------+
  | all-must-be-false | F    | All 'F' bits must be false.           |
  +-------------------+------+---------------------------------------+
  | don't care        | x    | Bit is not tested.                    |
  +-------------------+------+---------------------------------------+

Example: With a five-bit label, the selector ``label:+--++:FttTT`` will output
the k-mer if its label begins with ``0``, has at least one ``0`` in the next
two bits, and ends with ``11``.

Aliases ``all`` (all tests are ``T``), any (all tests are ``t``), and none
(all tests are ``F``) exist.

The default function is ``+`` and the default test is ``T``.

A selector ``label:0101011`` needs to be special case interpreted to
mean "the label equals 0101011".

A selector ``label:....011`` likewise should be special cased to mean "the
label ends in 011".

Examples on 2-bit labels:

.. table:: 

  +----------------+-----------------------+-------------------+
  | Selector       | Meaning               | Label Example     |
  |                |                       +----+----+----+----+
  |                |                       | 00 | 01 | 10 | 11 |
  +================+=======================+====+====+====+====+
  | label:00:all   | always false          |  F |  F |  F |  F |
  +----------------+-----------------------+----+----+----+----+
  | label:11:all   | always true           |  T |  T |  T |  T |
  +----------------+-----------------------+----+----+----+----+
  | label:-+:all   | label must be 01      |  F |  T |  F |  F |
  +----------------+-----------------------+----+----+----+----+
  | label:1+:all   | label must be x1      |  F |  T |  F |  T |
  +----------------+-----------------------+----+----+----+----+
  +----------------+-----------------------+----+----+----+----+
  | label:00:any   | always false          |  F |  F |  F |  F |
  +----------------+-----------------------+----+----+----+----+
  | label:11:any   | always true           |  T |  T |  T |  T |
  +----------------+-----------------------+----+----+----+----+
  | label:-+:any   | label cannot be 10    |  T |  T |  F |  T |
  +----------------+-----------------------+----+----+----+----+
  | label:0+:any   | label cannot be x0    |  F |  T |  F |  T |
  +----------------+-----------------------+----+----+----+----+
  +----------------+-----------------------+----+----+----+----+
  | label:00:none  | always true           |  T |  T |  T |  T |
  +----------------+-----------------------+----+----+----+----+
  | label:11:none  | always false          |  F |  F |  F |  F |
  +----------------+-----------------------+----+----+----+----+
  | label:-+:none  | label must be 10      |  F |  F |  T |  F |
  +----------------+-----------------------+----+----+----+----+
  | label:0+:none  | label cannot be x1    |  T |  F |  T |  F |
  +----------------+-----------------------+----+----+----+----+

Examples:

.. code:: none

  'equal' label:+-+-:all   --> true if the label is exactly 1010
  'all'   label:+1+1:any   --> true if none of the '+' bits are set
  'any'   label:+0+0:any   --> true if any of the '+' bits are set
  'none'  label:-1-1:all   --> true if none of the '-' bits are set
  'only'  label:1-1-:all   --> true if the '-' bits are all zero

Examples:

.. code:: none

  label:101011

    SPECIAL CASE, outputs k-mer if the label is exactly 101011.

  label:+-+-++:all

    Outputs k-mer if the label is exactly 101011.

  label=and label:+11+:all label:@2:1++1:all

    Compute the output label as the AND of all input labels.
    Require that the output label have the first and last bits set, AND
    require that the label on the second input have the middle two bits set.

Some tests that can be implemented:

A. ``L == C`` or ``L != C``

   For testing equality, convert the constant into a string of +'s (for 1
   bits) and -'s (for 0 bits) then check that all bits are set.

   For testing non-equality, invert the conversion so that 1's are set to -
   and 0's to +, then check that any bit is set.

   For testing (non-)equality of only a portion of the label, set the bits
   that should not be tested to 1 (for equality) or 0 (for non-equality).

   Proof: **not verified recently** if the label is the same as the constant,
   1 bits in the label will be inverted (to 0) and 0 bits in the label will
   be output true (so also 0) resulting in a modified label of all 0's.  If a
   1 bit in the label corresponds to a 0 bit in the constant, it will be
   passed true (a 1 in the modified label), similarly, a if a 0 bit in the
   label corresponds to a 1 bit in the constant it will be passed inverted (a
   1 in the modified label), either of which will make the modified label not
   equal to zero.

B. ``L & C == L`` or ``L | C == C``

   These test that C dominates L: that every bit set in L is also set in C, equivalently,
   that there is no bit set in L that is not set in C.

   Convert the constant c into a string of 0's (for 1 bits) and +'s (for 0
   bits), then check that no bit is set.

   Where C is a 1, we don't care what L is; 'l&c == l' is true for both l=0
   and l=1.  By forcing these bits to 0 in the modified string, they will
   never result in the check failing.

   Where C is a 0, howeveer, L must be 0.  Hence, passing the L bit true will
   result in a 1 bit output when L=1, which will cause the check to properly
   fail.  When L=0, a 0 is outpout, and the check passes.

   If instead of testing that no bit is set, we test that any bit is set, the
   sense of the selector is inverted; we test that C does not dominate L; that
   there is a bit set in L that is not set in C.

C. ``L & C == C`` or ``L | C == L``

   This is the dual of case B: L dominates C; that every bit set in C is also set in L.

   Convert the constant c into a string of +'s (for 1 bits) and 1's (for 0
   bits), then check that all bits are set.

   The selector is inverted if we test that some testable bit is false.
   
D. ``L & C == 0`` or ``L & C != 0``

   These test that L and C have no bits in common (or at least one bit in
   common).

   Convert the constant to +'s for 1's and 0's for 0's, then check that no
   bit is set (for no bits in common) or that some bit is set (for some bits
   in common).




Base Composition Selectors
--------------------------

The base composition selector selects k-mers for output based on the number of
A's, C's, G's and T's in the k-mer sequence.

.. code-block:: none

  bases:<BASES>:<OP><NUMBER>

Where ``<BASES>`` is a string containing ``A``, ``C``, ``G`` and ``T``
letters; case, order and quantity are unimportant.  The selector will count the
number of the specified letters in the k-mer and compare aginst ``<NUMBER>``
using the specified numeric comparison operator ``<OP>``.

.. table::
  :widths: 20 10 20 50

  +------------------+------+-------------------+--------------------------------------------------------+
  | BASES            |  OP  | NUMBER            | Meaning                                                |
  +==================+======+===================+========================================================+
  | ``A``            |      |                   | Count the number of ``A``'s in the k-mer.              |
  +------------------+------+-------------------+--------------------------------------------------------+
  | ``AC``           |      |                   | Count the number of ``A``'s and ``C``'s in the k-mer.  |
  +------------------+------+-------------------+--------------------------------------------------------+
  | ``GAAGAA``       |      |                   | Count the number of ``A``'s and ``G``'s in the k-mer.  |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |      | ``#n`` or ``n``   | Use the constant ``n``.                                |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``==``|                   | TRUE if the number of bases                            |
  |                  |``=`` |                   | is equal to the constant.                              |
  |                  |``eq``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``!=``|                   | TRUE if the number of bases                            |
  |                  |``<>``|                   | is not equal to the constant.                          |
  |                  |``ne``|                   |                                                        |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``<=``|                   | TRUE if the number of bases                            |
  |                  |``le``|                   | is less than or equal to the constant.                 |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``>=``|                   | TRUE if the number of bases                            |
  |                  |``ge``|                   | is greater than or equal to the constant.              |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``<`` |                   | TRUE if the number of bases                            |
  |                  |``lt``|                   | is less than the constant.                             |
  +------------------+------+-------------------+--------------------------------------------------------+
  |                  |``>`` |                   | TRUE if the number of bases                            |
  |                  |``gt``|                   | is greater than the constant.                          |
  +------------------+------+-------------------+--------------------------------------------------------+



Input Selectors
---------------

The input selector selects k-mers for output based on which inputs
supplied the k-mer.

.. code-block:: none

  input:<CONDITION>[:<CONDITION>[...]]

A ``<CONDITION>`` is either an input number (``@n``) or input count (``n`` or ``n-m``).
For the selector to be TRUE, all the CONDITIONS must be met.

.. warning::
  Do input-counts require ``#n`` or just integers ``n``?

Note that a k-mer is always present in at least one input.

Assuming 9 input files, some examples are:

.. table::
  :widths: 30 70

  +----------------+------------------------------------------------------+
  | Selector       | Output k-mer if it is present in...                  |
  +----------------+------------------------------------------------------+
  | input:@1       | ...the first input file.                             |
  +----------------+------------------------------------------------------+
  | input:@1-@3    | ...the first three input files.                      |
  +----------------+------------------------------------------------------+
  | input:#4:#5:@1 | ...4 or 5 input files, including the first           |
  +----------------+------------------------------------------------------+
  | input:#4-#6:#8 | ...4 or 5 or 6 or 8 input files.                     |
  +----------------+------------------------------------------------------+
  | input:#3-#9    | ...3 or more input files.                            |
  +----------------+------------------------------------------------------+
  | input:#1-#6    | ...at most 6 input files.                            |
  +----------------+------------------------------------------------------+

A few aliases exist:

.. table::
  :widths: 25 25 50

  +-------------+-------------+---------------------------------------------------------+
  | Alias       | Selector    | Meaning                                                 |
  +=============+=============+=========================================================+
  | input:any   | input:#1-#9 | k-mer is in any number of inputs.                       |
  +-------------+-------------+---------------------------------------------------------+
  | input:all   | input:#9    | k-mer is in all inputs.                                 |
  +-------------+-------------+---------------------------------------------------------+
  | input:only  | input:@1:#1 | k-mer is in the first input, and in exactly one input.  |
  +-------------+-------------+---------------------------------------------------------+
  | input:first | input:@1    | k-mer is in the first input, and maybe other inputs.    |
  +-------------+-------------+---------------------------------------------------------+

The difference between 'only' and 'first' is subtle: 'only' is true if the
k-mer is present exactly only in the firt file, while 'first' is true if the
k-mer exists in the first file and any other files.  'only' will effect a set
difference action, while 'first' is more akin to a set intersection.

================
Processing Trees
================

Meryl processes k-mers using a tree of actions.  An action reads k-mers from
multiple inputs, computes a function on the values and labels of all inputs
with the same k-mer, and outputs a single k-mer with a single value and a
single label.

(An action can also read sequence files and count the k-mers.)

Each action in the tree is enclosed in square brackets.  Square brackets
around the top-level / outermost action are optional.

The input to an action can be either a meryl database on disk or the output of
a different action.

The 'union' action below reads input from meryl databases 'input-1.meryl' and
'input-2.meryl'.  All three forms below are equivalent.

.. code-block:: none
  :caption: A simple union action reading from two inputs.
  :linenos:

  [ union input-1.meryl input-2.meryl ]

.. code-block:: none
  :caption: A simple union action reading from two inputs, but formatted.
  :linenos:

    union
      input-1.meryl
      input-2.meryl

.. code-block:: none
  :caption: A simple union action reading from two inputs, as sub-actions.
  :linenos:

    union
      [ input-1.meryl ]   //  This form technically makes input-1 and input-2 into
      [ input-2.meryl ]   //  sub-actions instead of direct inputs to 'union'.

Sub-actions can pre-process inputs.  The 'intersect' action below reads input
from two counting actions, and the one after computes a `union` before the
`intersection`.

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  intersect 
    [ count input-1.fasta output=input-1.meryl ]
    [ count input-2.fasta output=input-2.meryl ]

Each action will automatically pass its output k-mers to the parent action,
and can optionally write them to an output database.

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  intersect output=abINT12.meryl
    [ union input-a.meryl input-b.meryl output=ab.meryl ]
    [ union input-1.meryl input-2.meryl output=12.meryl ]

The original meryl allowed sub-actions to be supplied without surrounding
square brackets, but this led to great ambiguity in which action the output
modifier was associated with.  Without brackets, the following is ambiguous:

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  meryl
    union
      intersect
        a.meryl
        b.meryl
      intersect
        c.meryl
        d.meryl

As written, the intent is clear, but meryl interprets the second 'intersect' action
as an input to the first:

.. code-block:: none
  :caption: Sample databases.
  :linenos:

  meryl
    union
      intersect
        a.meryl
        b.meryl
        intersect
          c.meryl
          d.meryl

Therefore, meryl2 **requires** actions (except the very first) to be
surrounded by square brackets.
