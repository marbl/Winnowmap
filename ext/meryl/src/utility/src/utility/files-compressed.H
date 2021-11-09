
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

#ifndef FILES_COMPRESSED_H
#define FILES_COMPRESSED_H

//  Do not include directly.  Use 'files.H' instead.

enum cftType {
  cftNONE  = 0,
  cftGZ    = 1,
  cftBZ2   = 2,
  cftXZ    = 3,
  cftSTDIN = 4
};

cftType  compressedFileType(char const *filename);



class compressedFileReader {
public:
  compressedFileReader(char const *filename);
  ~compressedFileReader();

  void  reopen(void);

  FILE *operator*(void)     {  return(_file);              };
  FILE *file(void)          {  return(_file);              };

  char *filename(void)      {  return(_filename);          };

  bool  isCompressed(void)  {  return(_pipe == true);      };
  bool  isNormal(void)      {  return((_pipe == false) &&
                                      (_stdi == false));   };

private:
  FILE     *_file;
  char     *_filename;

  cftType   _type;

  bool      _pipe;
  bool      _stdi;
};



class compressedFileWriter {
public:
  compressedFileWriter(char const *filename, int32 level=1);
  ~compressedFileWriter();

  FILE *operator*(void)     {  return(_file);          };
  FILE *file(void)          {  return(_file);          };

  char *filename(void)      {  return(_filename);          };

  bool  isCompressed(void)  {  return(_pipe == true);  };

private:
  FILE  *_file;
  char  *_filename;
  bool   _pipe;
  bool   _stdi;
};



#endif  //  FILES_COMPRESSED_H
