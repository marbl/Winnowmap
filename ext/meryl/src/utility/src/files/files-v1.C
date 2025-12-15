
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

//#include "arrays.H"
//#include "strings.H"
#include "files.H"

using namespace merylutil::files::v1;



//  Searches for a file in the 'share/' directory.
//
//  Checks for $CANU_INSTALL_PATH/relpath/filename
//             $MERYL_INSTALL_PATH/relpath/filename
//             $PATH/../relpath/filename
//             ./filename
//  and returns the first one that exists.  If no file is found,
//  NULL is returned.
//
//  'relpath' should be something like 'share/sequence'.
//
char const *
findSharedFile(char const *binname, char const *relpath, char const *filename) {
  static
  char     fp[FILENAME_MAX + 1] = {0};
  char    *env;

  //  Does the file exist as is?

  if (fileExists(filename))
    return(filename);

  //  Does the file exist in any Canu installation?

  env = getenv("CANU_INSTALL_PATH");
  if (env != NULL) {
    snprintf(fp, FILENAME_MAX, "%s/%s/%s", env, relpath, filename);

    if (fileExists(fp))
      return(fp);
  }

  //  Does the file exist in any Meryl installation?

  env = getenv("MERYL_INSTALL_PATH");
  if (env != NULL) {
    snprintf(fp, FILENAME_MAX, "%s/%s/%s", env, relpath, filename);

    if (fileExists(fp))
      return(fp);
  }

  //  Does the file exist near our binary?

  strcpy(fp, binname);
  env = strrchr(fp, '/');
  if (env != nullptr) {
    *++env = 0;
    strcat(env, "../");
    strcat(env, relpath);
    strcat(env, "/");
    strcat(env, filename);

    if (fileExists(fp))
      return(fp);
  }

  //  Does the file exist in any component of the path?

  env = getenv("PATH");
  if (env != NULL) {
    while ((*env != ':') && (*env != 0)) {   //  Until we exhaust the PATH,
      int32  fpp = 0;

      while ((*env != ':') && (*env != 0))   //  Copy up to the first ':'
        fp[fpp++] = *env++;                  //  or the terminating NUL.

      if (*env == ':')                       //  Skip over the delimiter.
        env++;

      fp[fpp] = 0;

      strcat(fp, "/../");                    //  Append the relative path.
      strcat(fp, relpath);                   //  and file we're searching
      strcat(fp, "/");                       //  for.
      strcat(fp, filename);

      if (fileExists(fp))
        return(fp);
    }
  }

  //  Nope, not found.

  return(NULL);
}



