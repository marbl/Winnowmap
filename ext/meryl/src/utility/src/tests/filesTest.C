
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

using merylutil::compressedFileReader;
using merylutil::compressedFileWriter;

char tempname[64] = { 0 };
char tempnagz[64] = { 0 };

bool
testMkdirRmdir(void) {

  fprintf(stderr, "\n");
  fprintf(stderr, "Testing mkdir/rmdir.\n");

  if (merylutil::directoryExists(tempname) == true) {
    fprintf(stderr, " - Directory exists, rmdir it.\n");
    merylutil::rmdir(tempname);
  }
  if (merylutil::fileExists(tempname) == true) {
    fprintf(stderr, " - File exists, unlink it.\n");
    merylutil::unlink(tempname);
  }
  if (merylutil::pathExists(tempname) == true) {
    fprintf(stderr, " - Path exists, FAIL!\n");
    return false;
  }

  fprintf(stderr, " - mkdir.\n");
  if (merylutil::mkdir(tempname, false) == false) {
    fprintf(stderr, " - Failed to mkdir '%s'.\n", tempname);
    return false;
  }

  fprintf(stderr, " - rmdir.\n");
  if (merylutil::rmdir(tempname, false) == false) {
    fprintf(stderr, " - Failed to rmdir '%s'.\n", tempname);
    return false;
  }

  fprintf(stderr, " - mkdir over file.\n");
  merylutil::createEmptyFile(tempname);

  if (merylutil::mkdir(tempname, false) == false)
    fprintf(stderr, "   - correctly failed.\n");
  else
    fprintf(stderr, "   - INCORRECTLY SUCCEEDED.\n");

  fprintf(stderr, " - Cleanup.\n");
  merylutil::unlink(tempname);

  if (merylutil::pathExists(tempname) == true) {
    fprintf(stderr, "   - Failed to cleanup.\n");
    return false;
  }

  fprintf(stderr, " - Pass!\n");

  return true;
}


bool
testFileIO(uint16 *array, uint64 nObj) {

  fprintf(stderr, "\n");
  fprintf(stderr, "Testing compressedFileWriter/compressedFileReader using '%s'.\n", tempnagz);

  {
    compressedFileWriter *out = new compressedFileWriter(tempnagz);

    merylutil::writeToFile(array, "array", nObj, out->file());

    for (uint64 ii=0; ii<nObj; ii++)
      merylutil::writeToFile(array[ii], "value", out->file());

    delete out;
  }
  fprintf(stderr, " - data written.\n");

  if (merylutil::fileExists(tempnagz) == false) {
    fprintf(stderr, " - file doesn't exist.\n");
    return false;
  }

  {
    compressedFileReader *in = new compressedFileReader(tempnagz);

    merylutil::loadFromFile(array, "array", nObj, in->file());

    for (uint64 ii=0; ii<nObj; ii++)
      assert(array[ii] == (uint16)ii);

    for (uint64 ii=0; ii<nObj; ii++)
      merylutil::loadFromFile(array[ii], "value", in->file());

    for (uint64 ii=0; ii<nObj; ii++)
      assert(array[ii] == (uint16)ii);

    if (merylutil::loadFromFile(array[0], "value", in->file(), false) == 0)
      fprintf(stderr, " - eof detected.\n");
    else
      fprintf(stderr, " - eof NOT detected.\n");

    delete in;
  }
  fprintf(stderr, " - data read.\n");

  merylutil::unlink(tempnagz);

  fprintf(stderr, " - Pass!\n");

  return true;
}


bool
testUnlink(void) {

  fprintf(stderr, "\n");
  fprintf(stderr, "Testing unlink.\n");
  fprintf(stderr, " - of non-existent file\n");
  merylutil::unlink(tempname);

  fprintf(stderr, " - of existent file\n");
  merylutil::createEmptyFile(tempname);
  merylutil::unlink(tempname);

  fprintf(stderr, " - of directory\n");
  merylutil::mkdir(tempname);
  merylutil::unlink(tempname);
  merylutil::rmdir(tempname);

  if (merylutil::pathExists(nullptr) == true) {
    fprintf(stderr, "   - nullptr file exists?\n");
    return(false);
  }

  fprintf(stderr, " - Pass!\n");

  return true;
}


bool
testPermissions(void) {

  fprintf(stderr, "\n");
  fprintf(stderr, "Testing permissions.\n");

  if (merylutil::makeReadOnly(tempname) == true) {
    fprintf(stderr, "   - made non-existent file read only?\n");
    return false;
  }
  if (merylutil::makeWritable(tempname) == true) {
    fprintf(stderr, "   - made non-existent file write only?\n");
    return false;
  }

  merylutil::createEmptyFile(tempname);

  if (merylutil::fileExists(tempname, 0, nullptr, true) == false) {  //  TEST: should be writable
    fprintf(stderr, "   - initial file isn't writable.\n");
    return false;
  }

  if (merylutil::makeReadOnly(tempname) == false) {
    fprintf(stderr, "   - failed to make read only.\n");
    return false;
  }
  if (merylutil::fileExists(tempname, 0, nullptr, true) == true) {   //  TEST: should be read only
    fprintf(stderr, "   - file claims to be writable, but should be read-only.\n");
    return false;
  }

  if (merylutil::makeWritable(tempname) == false) {
    fprintf(stderr, "   - failed to make writable.\n");
    return false;
  }
  if (merylutil::fileExists(tempname, 0, nullptr, true) == false) {  //  TEST: should be writable
    fprintf(stderr, "   - file claims to be unwritable, but should be writable.\n");
    return false;
  }

  merylutil::unlink(tempname);

  fprintf(stderr, " - Pass!\n");

  return true;
}


int32
main(int32 argc, char **argv) {
  uint64     nObj      = (uint64)16 * 1024 * 1024;
  uint16    *array     = new uint16 [nObj];
  uint32     tests     = 0;
  bool       success   = true;

  strcpy(tempname, "./filesTest-XXXXXXXXXXXXXXXX");
  mktemp(tempname);
  strcpy(tempnagz, tempname);
  strcat(tempnagz, ".gz");

  std::vector<char const *>  err;
  for (int32 arg=1; arg < argc; arg++) {
    if      (strcmp(argv[arg], "-mkdir") == 0)        tests = 1;
    else if (strcmp(argv[arg], "-io") == 0)           tests = 2;
    else if (strcmp(argv[arg], "-unlink") == 0)       tests = 3;
    else if (strcmp(argv[arg], "-permissions") == 0)  tests = 4;
    else if (strcmp(argv[arg], "-suffix") == 0) {
      if (strlen(argv[++arg]) < 32) {
        strcpy(tempnagz, tempname);
        strcat(tempnagz, ".");
        strcat(tempnagz, argv[arg]);
      }
      else {
        fprintf(stderr, "error: suffix '%s' too long.\n", argv[arg]);
        return 1;
      }
    }
    else if (strcmp(argv[arg], "-h") == 0)            tests = 5;
    else                                              tests = 5;
  }
  if (tests == 5) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -mkdir        run just mkdir/rmdir tests.\n");
    fprintf(stderr, "  -io           run just compressed file create/read/write tests.\n");
    fprintf(stderr, "  -unlink       run just unlink tests,\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  -suffix suf   use suffix 'suf' for compressed files\n");
    fprintf(stderr, "                ('gz', 'lz', 'bz2', 'xz', 'zstd')\n");
    fprintf(stderr, "  \n");
    fprintf(stderr, "  by default, all tests are run.\n");
    fprintf(stderr, "  \n");
    return 0;
  }

  for (uint64 ii=0; ii<nObj; ii++)
    array[ii] = ii;

  fprintf(stderr, "Testing using temporary file/directory '%s' and '%s'.\n", tempname, tempnagz);

  if ((tests == 0) || (tests == 1))   success &= testMkdirRmdir();
  if ((tests == 0) || (tests == 2))   success &= testFileIO(array, nObj);
  if ((tests == 0) || (tests == 3))   success &= testUnlink();
  if ((tests == 0) || (tests == 4))   success &= testPermissions();

  delete [] array;

  if (success)
    fprintf(stderr, "\nAll tests passed!\n");
  else
    fprintf(stderr, "\nSome test failed.\n");

  return (success == false);
}

