
/******************************************************************************
 *
 *  This file is part of meryl-utility, a collection of miscellaneous code
 *  used by Meryl, Canu and others.
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

#include "system.H"
#include "arrays.H"
#include "files.H"

namespace merylutil::inline files::inline v1 {


compressedFileReader::compressedFileReader(const char *filename, int32 nThreads, cftType type) {
  _filename = duplicateString(filename);

  if (type == cftType::cftNONE)                //  If no type is supplied, guess it
    _type   = compressedFileType(_filename);   //  based on the extension (the default),
  else                                         //  otherwise, trust what the user told
    _type   = type;                            //  us (mostly for compressed pipe input).

  reopen(nThreads);
}


compressedFileReader::~compressedFileReader() {
  close();
  delete [] _filename;   _filename = nullptr;
  delete [] _line;       _line     = nullptr;
}


void
compressedFileReader::close(void) {

  errno = 0;

  if ((_file) && (_stdi == false) && (_pipe ==  true))   pclose(_file);
  if ((_file) && (_stdi == false) && (_pipe == false))   closeFile(_file, _filename);

  if (errno)
    fprintf(stderr, "WARNING:  Failed to cleanly close input file '%s': %s\n", _filename, strerror(errno));

  _file = nullptr;
}


void
compressedFileReader::reopen(int32 nThreads) {
  char   cmd[FILENAME_MAX];

  if (nThreads > 0)        //  Reset number of reading threads if a valid value
    _nThreads = nThreads;  //  is supplied; default to 2, set in .H.

  //  If and existing input is from stdin, do nothing.  reopen() on this
  //  makes no sense, and doing nothing is _possibly_ more correct than
  //  failing.
  if (_stdi)
    return;

  //  Close any existing file.
  close();

  //  Blow up if the file doesn't exist.
  if ((_type != cftSTDIN) && (fileExists(_filename) == false))
    fprintf(stderr, "ERROR:  Failed to open input file '%s': %s\n", _filename, strerror(ENOENT)), exit(1);

  //  Open the file!
  errno = 0;

  //  We used to allow pigz here, but some codes want to open many many input
  //  files and read them in parallel.  This uses far too many threads, and
  //  we've hit TOO MANY PROCESSES more than once.  So now we don't allow
  //  pigz.

  switch (_type) {
    case cftGZ:
      snprintf(cmd, FILENAME_MAX, "gzip -dc '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;
      break;

    case cftLZIP:
      snprintf(cmd, FILENAME_MAX, "lzip -dc '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;
      break;

    case cftBZ2:
      snprintf(cmd, FILENAME_MAX, "bzip2 -dc '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;
      break;

    case cftXZ:
      snprintf(cmd, FILENAME_MAX, "xz -dc '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;
      break;

    case cftZSTD:
      snprintf(cmd, FILENAME_MAX, "zstd -q -c -d -f -k '%s'", _filename);
      _file = popen(cmd, "r");
      _pipe = true;
      break;

    case cftSTDIN:
      _file = stdin;
      _stdi = true;
      break;

    default:
      _file = fopen(_filename, "r");
      _pipe = false;
      break;
  }

  //  Catch errors.
  //   - popen() does not set errno, so all we can do is fail.
  //   - otherwise, we can say something intelligent.

  if (_file == nullptr) {
    if (_pipe)
      fprintf(stderr, "ERROR:  Failed to open file with command '%s'\n", cmd);
    else
      fprintf(stderr, "ERROR:  Failed to open input file '%s': %s\n", _filename, strerror(errno));

    exit(1);
  }
}

}  //  merylutil::files::v1
