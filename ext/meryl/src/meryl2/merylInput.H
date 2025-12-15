
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

#ifndef MERYLINPUT_H
#define MERYLINPUT_H

#include "types.H"
#include "files.H"
#include "sequence.H"
#include "kmers.H"

#ifdef CANU                   //  If we're building under
#include "sqStore.H"          //  Canu, include support
#else                         //  to read the seqStore
typedef  uint32  sqStore;     //  directly, otherwise,
typedef  uint32  sqRead;      //  add fake types to reduce
#endif                        //  the number of ifdef blocks.


using namespace merylutil;
using merylutil::kmers::v2::merylFileReader;


class merylOpTemplate;
class merylOpCompute;


enum class merylInputType {                               // he's a real nowhere man
  inNowhere,    //  From nowhere, yet <------------------ \\ sitting in his nowhere land
  inTemplate,   //  From another template.                // making all his nowhere plans
  inAction,     //  From a compute slice.                 \\          . . .
  inPipe,       //  From a compute slice pipe.            // nowhere man please listen
  inDB,         //  From a database on disk.              \\ you don't know what you're missing
  inList,       //  From a meryl list output file.        // nowhere man, the world is at your command
  inSequence,   //  From a sequence file, for counting.   \\  -- John Lennon, 1965
  inCanu        //  From a sequence file, for counting.   //
};


//
//  Test if path p looks like a meryl database, a canu seqStore or just a
//  sequence file.
//
inline bool isMerylDatabase(char const *p) { return (fileExists(p, '/', "merylIndex")); }
inline bool isCanuSeqStore(char const *p)  { return (fileExists(p, '/', "info") &&
                                                     fileExists(p, '/', "reads"));      }
inline bool isSequenceFile(char const *p)  { return (fileExists(p));                    }



class merylInput {
public:
  merylInput() {};
  ~merylInput();

  //  The 'main' function - load the next kmer from the input source
  //  and save it internally so it can be returned in _kmer.

  void   nextMer(void);

  kmer  _kmer;
  bool  _valid     = false;

  //  Description of the input.

  merylInputType   _type = merylInputType::inNowhere;

  bool   isFromNowhere(void)      { return _type == merylInputType::inNowhere;  }
  bool   isFromTemplate(void)     { return _type == merylInputType::inTemplate; }
  bool   isFromAction(void)       { return _type == merylInputType::inAction;   }
  bool   isFromPipe(void)         { return _type == merylInputType::inPipe;     }
  bool   isFromDatabase(void)     { return _type == merylInputType::inDB;       }
  bool   isFromList(void)         { return _type == merylInputType::inList;     }
  bool   isFromSequence(void)     { return _type == merylInputType::inSequence; }
  bool   isFromStore(void)        { return _type == merylInputType::inCanu;     }

  bool   isCompressedFile(void)   { return isFromSequence() && _sequence->isCompressed(); }

  char const *inputName(void) {
    char const *n = nullptr;
    switch (_type) {
      case merylInputType::inNowhere:   n = "nowhere";         break;
      case merylInputType::inTemplate:  n = _templateName;     break;
      case merylInputType::inAction:    n = _actionName;       break;
      case merylInputType::inPipe:      n = _pipeName;         break;
      case merylInputType::inDB:        n = _dbName;           break;
      case merylInputType::inList:      n = _listName;         break;
      case merylInputType::inSequence:  n = _sequenceName;     break;
      case merylInputType::inCanu:      n = _storeName;        break;
      default:                          n = "unnamed-input";   break;
    }
    return n;
  }

  char const *inputType(void) {
    char const *t = nullptr;
    switch (_type) {
      case merylInputType::inNowhere:   t = "undefined-input";  break;
      case merylInputType::inTemplate:  t = "meryl-template";   break;
      case merylInputType::inAction:    t = "meryl-action";     break;
      case merylInputType::inPipe:      t = "meryl-pipe";       break;
      case merylInputType::inDB:        t = "meryl-database";   break;
      case merylInputType::inList:      t = "meryl-list-file";  break;
      case merylInputType::inSequence:  t = "sequence-file";    break;
      case merylInputType::inCanu:      t = "canu-seqStore";    break;
      default:                          t = "invalid-input";    break;
    }
    return t;
  }

  //
  //  Register an input.  Do checking only if the file(s) exist.
  //  Returns false if there is no way this can be opened.
  //

  bool        registerTemplate  (merylOpTemplate *t,                            std::vector<char const *> &err, bool displayErrors=true);
  bool        registerAction    (merylOpCompute *a,                             std::vector<char const *> &err, bool displayErrors=true);
  bool        registerActionPipe(char const *pipeName,                          std::vector<char const *> &err, bool displayErrors=true);
  bool        registerMerylDB   (char const *dbName,                            std::vector<char const *> &err, bool displayErrors=true);
  bool        registerMerylList (char const *listName,                          std::vector<char const *> &err, bool displayErrors=true);
  bool        registerSeqFile   (char const *seqName, bool doCompression,       std::vector<char const *> &err, bool displayErrors=true);
  bool        registerSeqStore  (char const *sqName, uint32 seg, uint32 segMax, std::vector<char const *> &err, bool displayErrors=true);

  merylInput *registerMerylDB   (char const *dbName);

  //
  //  Add per-slice inputs.
  //

  merylInput(merylOpTemplate *t)       { _type = merylInputType::inTemplate;  _template = t; }
  merylInput(merylOpCompute *a)        { _type = merylInputType::inAction;    _action   = a; }
  merylInput(merylFileReader *d)       { _type = merylInputType::inDB;        _db       = d; }
  merylInput(compressedFileReader *r)  { _type = merylInputType::inList;      _list     = r; }
  //rylInput(merylPipe *p)             { _type = merylInputType::inPipe;      _pipe     = p; }

  //id addActionSlice(merylOpCompute *a,     uint32 s)  { _type = merylInputType::inAction;  _action = a;  _sliceID = s; }
  //id addDbaseSlice(merylFileReader *d,     uint32 s)  { _type = merylInputType::inDB;      _db     = d;  _sliceID = s; }
  //id addListSlice(compressedFileReader *r, uint32 s)  { _type = merylInputType::inList;    _list   = r;  _sliceID = s; }
  //id addPipeSlice()   { _type = merylInputType::inPipe;    _pipe   = p; }

  //
  //  Once we're ready to start, this will actually open the inputs.
  //
private:
  void openinputSeqStore(void);
public:
  void openInput(std::vector<char const *> &err);

  //
  //  All the data!
  //

  //nt32                    _sliceID        = uint32max;

  //  Template input
  char const               *_templateName   = nullptr;
  merylOpTemplate          *_template       = nullptr;

  //  Action input
  char const               *_actionName     = nullptr;
  merylOpCompute           *_action         = nullptr;

  //  Action pipe input
  char const               *_pipeName       = nullptr;
  //rylOpCompute           *_pipe           = nullptr;

  //  Meryl database input
  char const               *_dbName         = nullptr;
  merylFileReader          *_db             = nullptr;

  //  Meryl list file input
  char const               *_listName       = nullptr;
  compressedFileReader     *_list           = nullptr;

  //  Sequence file input
  char const               *_sequenceName   = nullptr;
  dnaSeqFile               *_sequence       = nullptr;

  bool                      _squish         = false;      //  Sequence file inputs can
  char                      _lastByte       = 0;          //  be homopoly compressed.

  //  Canu seqStore input
  char const               *_storeName      = nullptr;
  sqStore                  *_store          = nullptr;
  uint32                    _storeSeg       = 0;
  uint32                    _storeSegMax    = 0;

  uint32                    _sqBgn          = 0;          //  First read to load.
  uint32                    _sqEnd          = 0;          //  Last read to load (not inclusive!)

  sqRead                   *_read           = nullptr;    //  Space to load the read.
  uint32                    _readID         = 0;
  uint32                    _readPos        = uint32max;



private:
  void   openInputSeqStore(void);
  bool   loadBasesFromCanu(char    *seq,
                           uint64   maxLength,
                           uint64  &seqLength,
                           bool    &endOfSequence);

public:
  bool   loadBases(char    *seq,
                   uint64   maxLength,
                   uint64  &seqLength,
                   bool    &endOfSequence);
};

#endif  //  MERYLINPUT_H
