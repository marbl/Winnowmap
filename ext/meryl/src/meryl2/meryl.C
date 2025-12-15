
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

//
//  Process each word on the command line:
//
//    Debug options and requests for help.
//     - dumpIndex
//     - dumpFile
//     - -h, -help, --help, help
//
//    Global options.  These are easier to do outside the command builder
//    (and besides, they're global options).
//      -k <kmer-size>
//      -l <label-size>
//      -m <mem-in-gb>    and -memory, --memory
//      -t <thread-limit> and -threads, --threads
//      -V, -VV, etc.
//  
//    Legacy options are k=<kmer-size>, memory=<mem-in-gb> and
//    threads=<thread-limit>.
//
//    If none of the above match, the usual case, toss the word to the
//    command builder and let it take care of it.
//
//  Once the command line has been scanned, finish building the command
//  trees, linking outputs to inputs, computing computable parameters, etc.
//  Then check for errors and display/fail if any are found.
//
//  Counting operations are a big headache.  They don't fit into the
//  tree nicely:
//   - they do their own threading, so only one thread can start the operation
//   - when done, they transform themselves into a pass-through operation that
//     simply reads the (just created) database and passes kmers through.
//
//  So, we special case them here.  This steps through each action in the tree,
//  counting kmers, writing to a new database, and finally converting the action
//  to a null pass-through.
//
//  Once counting is done, the action tree is expanded into 64 copies, one
//  for each database slice, that are run in parallel.  After the run,
//  deleting the commandBuilder will (recursively) delete all the action
//  nodes which will close and open files, etc.
//



//  Parse text from a file into our pTxt string.
//
//  Single and double quotes behave as expected, joining
//  multiple words into a single token and allowing the other quote
//  mark.  These are all single words and the outer-most quotes are
//  removed:
//    "one word"      -> |one word|
//    "it's here"     -> |it's here|
//    'sample "A"'    -> |sample "A"|
//
//  Escapes can be used for the same effect:
//     one\ word         |one word|
//     it\'s here        |it's here|
//     sample\ \"A\"  -> |sample "A"|
//
//  Note that the escape symbol '\' is treated as plain text
//  in a quoted string:
//    "one\ word"     -> |one\ word"
//  
//  An escape at the end of a line 

void
merylCommandBuilder::loadProgramText(char const *f) {
  bool        esc = false;   //  Escape mode
  bool        sgl = false;   //  In a single-quoted word
  bool        dbl = false;   //  In a double-quoted word

  compressedFileReader  *cft = new compressedFileReader(f);

  while (cft->readLine() == true) {
    char const *line  = cft->line();
    uint32      lineL = cft->lineLen();

    //  Grow the array by 16 KB if appending the new line would exceed allocated space.
    increaseArray(_pTxt, _pTxtLen+lineL+2, _pTxtMax, 16384);

    for (uint32 ll=0; ll < lineL; ll++) {
      char  ch = line[ll];

      bool  nesc = ((esc == false) && (sgl == false) && (dbl == false));
      bool  osgl = ((esc == false) && (sgl == true)  && (dbl == false));
      bool  odbl = ((esc == false) && (sgl == false) && (dbl == true ));

      bool  com1 = ((ll == 0) &&                        (line[ll] == '#'));
      bool  com2 = ((ll >  1) && (line[ll-1] == ' ') && (line[ll] == '#'));

      //  Handle mode switches.
      //
      //   - If not in an escape mode and a backslash, single or double quote
      //     is encountered, enter the appropriate mode.
      //   - If in a quote mode and a close quote is seen, exit the mode.
      //
      if      (nesc && (ch == '\\'))  esc = true;
      else if (nesc && (ch == '\''))  sgl = true;
      else if (nesc && (ch == '"'))   dbl = true;
      else if (osgl && (ch == '\''))  sgl = false;
      else if (odbl && (ch == '"'))   dbl = false;

      //  Handle comments.  Just bail on the rest of this line.
      //
      else if (nesc && (com1 || com2))
        break;

      //  Add a letter.
      //
      //   - If not in an escape/quote mode, replace spaces with word-separators.
      //   - Otherwise, add the letter verbatim and exit escape mode.
      //
      else if (nesc && (ch == ' '))  _pTxt[_pTxtLen++] = '\v';
      else {
        _pTxt[_pTxtLen++] = ch;
        esc = false;
      }
    }

    if ((sgl == true) || (dbl == true))
      fprintf(stderr, "WARNING: end-of-line encountered in quoted string.\n");

    if (esc == false)             //  Add a word-sep at the end of each line,
      _pTxt[_pTxtLen++] = '\v';   //  unless it is escaped.
    esc = false;                  //  Exit escape mode.

    _pTxt[_pTxtLen]   = '\0';
  }

  delete cft;
}


//  Words are separated by a vertical-tab, VT, \v.
//
//  Words are appended to the program text verbatim.  This will
//  pass filenames correctly, but will NOT pass quoted command
//  line options correctly:
//    meryl2 "output : database " <file>
//  but this is only a problem for malicious users, and we don't
//  support those.
//
void
merylCommandBuilder::appendProgramWord(char const *w) {

  if (w == nullptr)
    return;

  for (uint32 l=strlen(w)+2; _pTxtLen+l >= _pTxtMax; ) //  Get space for the word
    increaseArray(_pTxt, _pTxtLen+l, _pTxtMax, 16384); //  and terminating bytes.

  while (*w)                                           //  Copy word to program text.
    _pTxt[_pTxtLen++] = *w++;                          //

  _pTxt[_pTxtLen++] = '\v';                            //  Separate the word from any next
  _pTxt[_pTxtLen]   = '\0';                            //  word and terminate the string.
}


int
main(int argc, char **argv) {
  merylCommandBuilder  *B = new merylCommandBuilder;
  int                   r = 0;

  argc = globals.initialize(argc, argv);  //  Handles --version, initializes max memory and threads.

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if ((globals.processDebugOption (arg, argv, err) == true) ||   //  The process() function
        (globals.processGlobalOption(arg, argv, err) == true))     //  handles the option;
      ;                                                            //  nothing for us to do.
    else if (strcmp(argv[arg], "-f") == 0)
      B->loadProgramText(argv[++arg]);
    else
      B->appendProgramWord(argv[arg]);
  }

  B->processProgramText();

  B->finalizeTrees();   //  Finalize the processing tree and check for errors.

  //
  //  Report command line usage if any errors so far.
  //
  if ((argc == 1) ||                //  No commands
      (B->numTrees() == 0) ||       //  No actions
      (err.size() > 0)) {           //  Command line errors
    fprintf(stderr, "usage: %s ...\n", argv[0]);

#include "meryl-usage.H"

    for (uint32 ii=0; ii<err.size(); ii++)
      if (err[ii] != NULL)
        fprintf(stderr, "%s\n", err[ii]);

    r=1;  //  Exit fail.
  }

  //
  //  Or report the tree and fail if there are errors.
  //
  else if (B->displayTreesAndErrors() == true) {
    r=1;
  }

  //
  //  Or just stop after showing a successful tree.
  //
  else if (globals.stopAfterConfigure()) {
    fprintf(stderr, "Stopping after configuration.\n");
  }

  //
  //  Or actually run!
  //
  else {
    B->performCounting(globals.allowedMemory(), globals.allowedThreads());
    B->finalizeParameters();
    B->spawnThreads(globals.allowedThreads());
    B->runThreads(globals.allowedThreads());

    if (globals.showStandard() == true)
      fprintf(stderr, "\n"
                      "Cleaning up.\n");
  }

  delete B;

  if (globals.showStandard() == true)
    fprintf(stderr, "\n"
                    "Bye.\n");
  return r;
}
