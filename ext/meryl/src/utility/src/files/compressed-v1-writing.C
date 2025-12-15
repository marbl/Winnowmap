
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


compressedFileWriter::compressedFileWriter(const char *filename, uint32 cLevel, uint32 nThreads) {
  char   cmd[FILENAME_MAX];

  _file     = NULL;
  _filename = duplicateString(filename);
  _pipe     = false;
  _stdi     = false;

  cftType   ft = compressedFileType(_filename);

  if (nThreads == 0)             //  Use max threads if an invalid
    nThreads = getNumThreads();  //  number is requested.

  errno = 0;

  switch (ft) {
    case cftGZ:                          //  If 'pigz' looks like it works, use
      if (commandAvailable("pigz -h"))   //  that with a few threads.
        snprintf(cmd, FILENAME_MAX, "pigz -%dc -p %d > '%s'", cLevel, nThreads, _filename);
      else
        snprintf(cmd, FILENAME_MAX, "gzip -%dc > '%s'", cLevel, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftLZIP:
      if (commandAvailable("plzip -h"))
        snprintf(cmd, FILENAME_MAX, "plzip -%dc -n %d > '%s'", cLevel, nThreads, _filename);
      else
        snprintf(cmd, FILENAME_MAX, "lzip -%dc > '%s'", cLevel, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftBZ2:
      snprintf(cmd, FILENAME_MAX, "bzip2 -%dc > '%s'", cLevel, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftXZ:
      snprintf(cmd, FILENAME_MAX, "xz -%dc > '%s'", cLevel, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftZSTD:
      snprintf(cmd, FILENAME_MAX, "zstd -q -c -%d -T%d > '%s'", cLevel, nThreads, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftSTDIN:
      _file = stdout;
      _stdi = 1;
      break;

    default:
      _file = fopen(_filename, "w");
      _pipe = false;
      break;
  }

  if (errno)
    fprintf(stderr, "ERROR:  Failed to open output file '%s': %s\n", _filename, strerror(errno)), exit(1);
}


compressedFileWriter::~compressedFileWriter() {
  close();
  delete [] _filename;   _filename = nullptr;
}


void
compressedFileWriter::close(void) {

  errno = 0;

  if ((_file) && (_stdi == false) && (_pipe ==  true))   pclose(_file);
  if ((_file) && (_stdi == false) && (_pipe == false))   closeFile(_file, _filename);

  if (errno)
    fprintf(stderr, "ERROR:  Failed to cleanly close output file '%s': %s\n", _filename, strerror(errno)), exit(1);

  _file = nullptr;
}


}  //  merylutil::files::v1
