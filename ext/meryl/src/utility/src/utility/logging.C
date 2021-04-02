
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

#include "logging.H"
#include "arrays.H"

#include <stdarg.h>





class logFileLevel {
public:
  logFileLevel(char const *levelName, bool levelEnabled=false) {
    setName(levelName);
    _lVerbosity = 0;
    _lEnabled   = levelEnabled;
  };

  ~logFileLevel() {
  };

  bool   isEnabled(uint32 messageVerbosity=0) {
    return((true             == _lEnabled) &&
           (messageVerbosity <= _lVerbosity));
  };

  uint32 verbosity(void) {
    return(_lVerbosity);
  }


private:
  void  setName(char const *levelName) {
    uint32 x=0;

    for (; (levelName != NULL) && (x < logFileLevelNameLen) && (levelName[x] != 0); x++)
      _lName[x] = levelName[x];

    for (; x < logFileLevelNameLen; x++)
      _lName[x] = 0;
  };

public:
  bool      operator==(char const *name) const {
    return(strncmp(_lName, name, logFileLevelNameLen) == 0);
  };

  bool      operator<(const logFileLevel &that) const {
    return(strncmp(_lName, that._lName, logFileLevelNameLen) < 0);
  };

private:
  char     _lName[logFileLevelNameLen];    //  Name of this logging class
  uint32   _lVerbosity;                    //  Strength of messages to show (incremented by logFile itself)
  bool     _lEnabled;                      //  Enabled if true

  friend class logFile;
};






class logFileInstance {
public:
  logFileInstance(char   const *prefix,
                  uint32 const  threadID,
                  uint32 const  bufferSize = 1024 * 1024) {

    memset(_prefix,     0, sizeof(char) * (FILENAME_MAX + 1));
    memset(_name,       0, sizeof(char) * (FILENAME_MAX + 1));
    memset(_filePrefix, 0, sizeof(char) * (FILENAME_MAX + 1));
    memset(_fileName,   0, sizeof(char) * (FILENAME_MAX + 1));

    if (prefix)
      strncpy(_prefix, prefix, FILENAME_MAX);

    _order      = 0;
    _threadID   = threadID;

    _part       = 0;

    _length     = 512 * 1024 * 1024;   //  Forces a rotate() on the first write.
    _lengthMax  = 512 * 1024 * 1024;

    _output     = NULL;
  };


  ~logFileInstance() {
    delete _output;
  };


  void    setPrefix(char const *prefix) {

    memset(_prefix, 0, sizeof(char) * (FILENAME_MAX + 1));

    if (prefix)
      strncpy(_prefix, prefix, FILENAME_MAX);
  };


  void    setMaxSize(uint64 maxsize) {
    _lengthMax = maxsize;
  };


  void    setName(char const *name) {

    memset(_name, 0, sizeof(char) * (FILENAME_MAX + 1));

    if (name)
      strncpy(_name, name, FILENAME_MAX);

    _order += 1;
    _part   = 0;

    rotate();
  };


  //  Rotate the log file please, HAL.
  void    rotate(void) {

    delete _output;

    _part++;

    if (_prefix[0] == 0)
      fprintf(stderr, "_prefix not set for thread %d\n", _threadID);
    assert(_prefix[0] != 0);

    if (_threadID < UINT32_MAX) {
      snprintf(_filePrefix, FILENAME_MAX, "%s.%03u.%s",         _prefix, _order, _name);
      snprintf(_fileName,   FILENAME_MAX, "%s.%03u.%s.thr%03d", _prefix, _order, _name, _threadID);
    } else {
      snprintf(_filePrefix, FILENAME_MAX, "%s.%03u.%s",         _prefix, _order, _name);
      snprintf(_fileName,   FILENAME_MAX, "%s.%03u.%s",         _prefix, _order, _name);
    }

    _length = 0;
    _output = new writeBuffer(_fileName, "w", _lengthMax);
  };


  void    rotateMessage(void) {
    char  rotmes[256] = {0};

    if (_output) {
      snprintf(rotmes, 1024, "logFile()--  next message will exceed file size limit of " F_U64 "; rotate to new file.\n",
               _lengthMax);

      _output->write(rotmes, strlen(rotmes));
    }
  };


  void    writeLog(char const *fmt, va_list ap) {

    errno = 0;

    char   *mes = NULL;
    int32   len = vasprintf(&mes, fmt, ap);

    if (len == -1) {
      fprintf(stderr, "writeLog()-- error writing log with fmt '%s': %s\n", fmt, strerror(errno));
      return;
    }

    if (_length + len > _lengthMax) {
      rotateMessage();
      rotate();
    }

    _output->write(mes, len);

    free(mes);
  };


  void    flush(void) {
    if (_output)
      _output->flush();
  }


private:
  char          _prefix[FILENAME_MAX + 1];
  char          _name  [FILENAME_MAX + 1];

  char          _filePrefix[FILENAME_MAX + 1];   //  e.g., 'prefix.###.name'
  char          _fileName  [FILENAME_MAX + 1];   //  e.g., 'prefix.###.name.thr###.part###.log'

  uint32        _order;
  uint32        _threadID;

  uint32        _part;

  writeBuffer  *_output;
  uint64        _length;
  uint64        _lengthMax;

  friend class logFile;
};






logFile::logFile(char const *prefix, uint64 maxSize) {

  _threadMax = 1024;
  _threadNum = omp_get_max_threads();

  _maxSize   = maxSize;

  _mainI     = new logFileInstance(prefix, UINT32_MAX, maxSize);
  _threadI   = new logFileInstance * [_threadMax];

  for (uint32 ii=0; ii<_threadMax; ii++)
    _threadI[ii] = nullptr;

  for (uint32 ii=0; ii<_threadNum; ii++)
    _threadI[ii] = new logFileInstance(prefix, ii, maxSize);

  _levelsLen = 0;
  _levelsMax = 16;
  _levels    = new logFileLevel * [_levelsMax];

  _verbosity = 0;
}


logFile::~logFile() {

  delete    _mainI;

  for (uint32 ii=0; ii<_threadMax; ii++)
    delete _threadI[ii];
  delete [] _threadI;

  for (uint32 ii=0; ii<_levelsLen; ii++)
    delete _levels[ii];
  delete [] _levels;
}


void
logFile::setPrefix(char const *prefix) {

  _mainI->setPrefix(prefix);

  for (uint32 ii=0; ii<_threadMax; ii++) {
    if (_threadI[ii])
      _threadI[ii]->setPrefix(prefix);
  }
}


char const *
logFile::getPrefix(void) {
  return(_mainI->_prefix);
}


char const *
logFile::getLogName(void) {
  return(_mainI->_filePrefix);
}



void
logFile::setName(char const *name) {

  _mainI->setName(name);

  for (uint32 ii=0; ii<_threadMax; ii++)
    if (_threadI[ii])
      _threadI[ii]->setName(name);
}


void
logFile::setMaxSize(uint64 size) {

  _mainI->setMaxSize(size);

  for (uint32 ii=0; ii<_threadMax; ii++)
    if (_threadI[ii])
      _threadI[ii]->setMaxSize(size);
}


logFileHandle
logFile::addLevel(char const *levelName,
                  bool        enabled) {
  logFileLevel   key(levelName);

  increaseArray(_levels, _levelsLen, _levelsMax, 8);

  //_levelsIndex[key]   = _levelsLen;
  _levels[_levelsLen] = new logFileLevel(levelName, enabled);

  _levelsLen++;

  return(logFileHandle(_levelsLen-1));
}


logFileHandle
logFile::addLevel(char const *levelName,
                  uint32      verbosity,
                  bool        enabled) {
  logFileLevel   key(levelName);

  increaseArray(_levels, _levelsLen, _levelsMax, 8);

  //_levelsIndex[key]   = _levelsLen;
  _levels[_levelsLen] = new logFileLevel(levelName, enabled);

  _levelsLen++;

  return(logFileHandle(_levelsLen-1));
}




//  We're expecting a few tens of levels maximum, and so just do a linear
//  search on the levels we know to find the one that matches the supplied
//  name.  (Plus, using a map<> makes the client code below much much more
//  complicated.)
//
uint32
logFile::findLevelIndex(char const *levelName) {
  uint32  ret = UINT32_MAX;

  for (uint32 ii=0; ii<_levelsLen; ii++)
    if (strcmp(levelName, _levels[ii]->_lName) == 0)
      ret = ii;

  return(ret);
}


int32
logFile::enable(char const *optionString, char const *levelName) {
  uint32    verbosity = 0;

  while ((*optionString != 0) && (*optionString == '-')) {
    optionString++;
  }

  while ((*optionString != 0)) {
    verbosity++;
    optionString++;
  }

  if (levelName == NULL) {
    _verbosity = verbosity;
    return(0);
  }

  else {
    enable(levelName, verbosity);
    return(1);
  }
}


void
logFile::enable(char const *levelName, uint32 verbosity) {
  uint32   idx = findLevelIndex(levelName);

  if (idx < _levelsLen) {
    _levels[idx]->_lEnabled   = true;
    _levels[idx]->_lVerbosity = verbosity;
  }
}

void
logFile::disable(char const *levelName) {
  uint32   idx = findLevelIndex(levelName);

  if (idx < _levelsLen) {
    _levels[idx]->_lEnabled = false;
  }
}

void
logFile::increment(char const *levelName) {
  uint32   idx = findLevelIndex(levelName);

  if (idx < _levelsLen) {
    _levels[idx]->_lVerbosity++;
  }
}





void
logFile::enable(logFileHandle levelName, uint32 verbosity) {
  uint32   idx = levelName._index;

  if (idx < _levelsLen) {
    _levels[idx]->_lEnabled   = true;
    _levels[idx]->_lVerbosity = verbosity;
  }
}

void
logFile::disable(logFileHandle levelName) {
  uint32   idx = levelName._index;

  if (idx < _levelsLen) {
    _levels[idx]->_lEnabled = false;
  }
}

void
logFile::increment(logFileHandle levelName) {
  uint32   idx = levelName._index;

  if (idx < _levelsLen) {
    _levels[idx]->_lVerbosity++;
  }
}





void
logFile::enable(uint32 verbosity) {
  _verbosity = verbosity;
}

void
logFile::increment(void) {
  _verbosity++;
}






void
logFile::writeStatus(char const *fmt, va_list ap) {
  vfprintf(stderr, fmt, ap);
}


void
logFile::writeLog(char const *fmt, va_list ap) {
  int32   nt = omp_get_num_threads();   //  Number of threads currently active
  int32   tn = omp_get_thread_num();    //  ID of this thread

  //  If tn is more than we have space for we need to allocate a new
  //  _threadI array.  But this is hard.  So just blow up.

  if (tn >= _threadMax) {
    fprintf(stderr, "TOO MANY THREADS!\n");
    assert(0);
  }

  //  If we're only running a single thread, or we have already allocated an
  //  output for this thread, we can immediately write the log.

  if (nt == 1) {
    _mainI->writeLog(fmt, ap);
  }

  else if (_threadI[tn]) {
    _threadI[tn]->writeLog(fmt, ap);
  }

  //  Otherwise, we need to allocate a new thread output and set it up before
  //  we can write.

  else {
    _threadI[tn] = new logFileInstance(getPrefix(), tn, _maxSize);
    _threadI[tn]->writeLog(fmt, ap);
  }
}




void
logFile::writeStatus(char const *fmt, ...) {
  va_list           ap;

  va_start(ap, fmt);
  writeStatus(fmt, ap);
  va_end(ap);
}

void
logFile::writeLog(char const *fmt, ...) {
  va_list           ap;

  va_start(ap, fmt);
  writeLog(fmt, ap);
  va_end(ap);
}

void
logFile::writeStatus(logFileHandle levelName, char const *fmt, ...) {
  va_list           ap;

  if (levelEnabled(levelName)) {
    va_start(ap, fmt);
    writeStatus(fmt, ap);
    va_end(ap);
  }
}

void
logFile::writeLog(logFileHandle levelName, char const *fmt, ...) {
  va_list           ap;

  if (levelEnabled(levelName)) {
    va_start(ap, fmt);
    writeLog(fmt, ap);
    va_end(ap);
  }
}


inline
bool
logFile::verbosityEnabled(uint32 verbosity) {   //  True if message 'verbosity' is
  return(verbosity <= _verbosity);              //  below our max '_verbosity' allowed.
}

inline
bool
logFile::levelEnabled(logFileHandle level, uint32 verbosity) {
  uint32  idx = level._index;

  if (idx < _levelsLen)
    return(_levels[idx]->isEnabled(verbosity));

  return(false);
}



void
logFile::writeStatus(uint32 verbosity, char const *fmt, ...) {
  va_list           ap;

  if (verbosityEnabled(verbosity)) {
    va_start(ap, fmt);
    writeStatus(fmt, ap);
    va_end(ap);
  }
}

void
logFile::writeLog(uint32 verbosity, char const *fmt, ...) {
  va_list           ap;

  if (verbosityEnabled(verbosity)) {
    va_start(ap, fmt);
    writeLog(fmt, ap);
    va_end(ap);
  }
}

void
logFile::writeStatus(logFileHandle levelName, uint32 verbosity, char const *fmt, ...) {
  va_list           ap;

  if (levelEnabled(levelName, verbosity)) {
    va_start(ap, fmt);
    writeStatus(fmt, ap);
    va_end(ap);
  }
}

void
logFile::writeLog(logFileHandle levelName, uint32 verbosity, char const *fmt, ...) {
  va_list           ap;

  if (levelEnabled(levelName, verbosity)) {
    va_start(ap, fmt);
    writeLog(fmt, ap);
    va_end(ap);
  }
}

void
logFile::flush(void) {
  _mainI->flush();

  for (uint32 ii=0; ii<_threadMax; ii++)
    if (_threadI[ii])
      _threadI[ii]->flush();
}

