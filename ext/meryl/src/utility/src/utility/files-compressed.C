
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

#include "files.H"
#include "arrays.H"


cftType
compressedFileType(char const *filename) {

  if ((filename == NULL) || (filename[0] == 0) || (strcmp(filename, "-") == 0))
    return(cftSTDIN);

  int32  len = strlen(filename);

  if      ((len > 3) && (strcasecmp(filename + len - 3, ".gz") == 0))
    return(cftGZ);

  else if ((len > 4) && (strcasecmp(filename + len - 4, ".bz2") == 0))
    return(cftBZ2);

  else if ((len > 3) && (strcasecmp(filename + len - 3, ".xz") == 0))
    return(cftXZ);

  else
    return(cftNONE);
}



compressedFileReader::compressedFileReader(const char *filename) {
  char    cmd[FILENAME_MAX];
  int32   len = 0;

  _file     = NULL;
  _filename = duplicateString(filename);
  _pipe     = false;
  _stdi     = false;

  cftType   ft = compressedFileType(_filename);

  if ((ft != cftSTDIN) && (fileExists(_filename) == false))
    fprintf(stderr, "ERROR:  Failed to open input file '%s': %s\n", _filename, strerror(errno)), exit(1);

  errno = 0;

  switch (ft) {
    case cftGZ:
      snprintf(cmd, FILENAME_MAX, "gzip -dc '%s'", _filename);
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

      if (_file == NULL)    //  popen() returns NULL on error.  It does not reliably set errno.
        fprintf(stderr, "ERROR:  Failed to open input file '%s': popen() returned NULL\n", _filename), exit(1);

      errno = 0;
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

  if (errno)
    fprintf(stderr, "ERROR:  Failed to open input file '%s': %s\n", _filename, strerror(errno)), exit(1);
}



compressedFileReader::~compressedFileReader() {

  if (_file == NULL)
    return;

  if (_stdi)
    return;

  if (_pipe)
    pclose(_file);
  else
    AS_UTL_closeFile(_file);

  delete [] _filename;
}



compressedFileWriter::compressedFileWriter(const char *filename, int32 level) {
  char   cmd[FILENAME_MAX];

  int32  nThreads      = omp_get_max_threads();
  bool   pigzAvailable = false;

  _file     = NULL;
  _filename = duplicateString(filename);
  _pipe     = false;
  _stdi     = false;

  cftType   ft = compressedFileType(_filename);

  //  Decide if we have pigz or gzip available.

  if (ft == cftGZ) {
    snprintf(cmd, FILENAME_MAX, "pigz -h > /dev/null 2>&1");

    FILE *F = popen(cmd, "r");
    int32 e = pclose(F);

    if (e == 0)
      pigzAvailable = true;
  }

#if 0
  if (pigzAvailable)
    fprintf(stderr, "Using pigz for compression.\n");
  else
    fprintf(stderr, "Using gzip for compression.\n");
#endif

  //  Open the output processor for input.

  errno = 0;

  switch (ft) {
    case cftGZ:
      if (pigzAvailable)
        snprintf(cmd, FILENAME_MAX, "pigz -%dc -p %d > '%s'", level, nThreads, _filename);
      else
        snprintf(cmd, FILENAME_MAX, "gzip -%dc > '%s'", level, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftBZ2:
      snprintf(cmd, FILENAME_MAX, "bzip2 -%dc > '%s'", level, _filename);
      _file = popen(cmd, "w");
      _pipe = true;
      break;

    case cftXZ:
      snprintf(cmd, FILENAME_MAX, "xz -%dc > '%s'", level, _filename);
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

  if (_file == NULL)
    return;

  if (_stdi)
    return;

  errno = 0;

  if (_pipe)
    pclose(_file);
  else
    AS_UTL_closeFile(_file);

  if (errno)
    fprintf(stderr, "ERROR:  Failed to cleanly close output file '%s': %s\n", _filename, strerror(errno)), exit(1);

  delete [] _filename;
}
