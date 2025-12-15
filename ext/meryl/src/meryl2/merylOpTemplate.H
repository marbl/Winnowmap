
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

#ifndef MERYLOPTEMPLATE_H
#define MERYLOPTEMPLATE_H

#ifndef MERYLINCLUDE
#error "Do not use merylOpTemplate.H, use meryl.H instead."
#endif


//  This class is used for setting up the processing tree, checking args,
//  keeping track of parameters, etc.  Objects are made/modified while
//  the command line is being parsed.  Once parsing is done, they are
//  used to make merylOpCompute objects for each thread.
//
//  Counting operations are also done here.

//  The selectors are a sum-of-products expression returning TRUE if the kmer should
//  be output and FALSE if the kmer should be ignored.
//
//  isKmerSelectoredOut() computes the sum of those product selectors applied
//  to the kmer we are thinking of outputting (_kmer).
//
//  Ideally, this should be captured in its own object, but it needs access
//  to the _act lists.
//
//  isTrue(k) returns true if the output kmer 'k' meets the conditions of
//  the selector and should be passed to the output.  It has access to all
//  the other input kmers through the _act inputs.
//
//  addNewSelectorProduct() will add a new product term to the sum list.
//  addSelectorToProduct()  will add a new selector to the most recent product term.

class merylOpTemplate {
public:
public:
  merylOpTemplate(uint32 ident);
  ~merylOpTemplate();

public:
  bool    isEmpty(void) {
    return ((_isCounting == false) && (_counting == nullptr) &&
            (_isSelector == false) &&
            (_inputs.size() == 0) &&
            (_select.size() == 0) &&
            (_outDbse == nullptr) &&
            (_outList == nullptr) &&
            (_outShow == nullptr) &&
            (_outStats == nullptr) &&
            (_outHisto == nullptr));
  }

public:
  void    addInput(merylInput *in) {;
    _inputs.push_back(in);
  }

  void    addOutputToDB  (char const *wrName,                 std::vector<char const *> &err);
  void    addOutputToList(char const *prName, bool ACGTorder, std::vector<char const *> &err);
  void    addOutputToPipe(char const *pipeName,               std::vector<char const *> &err) {};

  void    addStatsOutput (char const *hiName,                 std::vector<char const *> &err);
  void    addHistoOutput (char const *hiName,                 std::vector<char const *> &err);

public:
  bool    acceptsSequenceInputs(void)  { return _isCounting; };

public:
  void    finalizeTemplateInputs(uint32 oo, std::vector<char const *> &err);
  void    finalizeTemplateParameters(void);

  void    finishAction(void);

  //
  //  Counting operations.
  //

  void    doCounting(uint64 allowedMemory,
                     uint32 allowedThreads);

  merylOpCounting               *_counting = nullptr;

  //
  //  Inputs.  The minimum and maximum are set by the action or selector:
  //    Histograms must have exactly one input
  //    Subtract needs two
  //

  std::vector<merylInput *>      _inputs;

  uint32                         _inputsMin = 0;
  uint32                         _inputsMax = uint32max;

  //
  //  Action parameters, a description of how to generate output values/labels.
  //

  uint32                         _ident         = uint32max;
  bool                           _isSelector    = false;
  bool                           _isCounting    = false;

  bool                           _onlyConfig    = false;

  merylOpCompute                *_computes[64]  = { nullptr };

  void                            setValueConstant(char const *c, merylAssignValue a, kmvalu v);
  char                          *_valueString   = nullptr;
  merylAssignValue               _valueAssign   = merylAssignValue::valueNOP;
  kmvalu                         _valueConstant = 0;

  void                            setLabelConstant(char const *c, merylAssignLabel a, kmlabl l);
  char                          *_labelString   = nullptr;
  merylAssignLabel               _labelAssign   = merylAssignLabel::labelNOP;
  kmlabl                         _labelConstant = 0;

  char const     *displayValueAssignment(char *ts);
  char const     *displayLabelAssignment(char *ts);

  //
  //  Selectors.
  //

  void            addNewSelectorProduct(void);
  void            addSelectorToProduct(merylSelector const f);
  bool            isSelectorProductEmpty(void);
  merylSelector  &getSelector(uint32 p, uint32 t);
  merylSelector  &getLastSelector(void);

private:
  std::vector< std::vector<merylSelector> >  _select;


  //
  //  OUTPUTS.
  //
public:
  char const                    *_outDbseName    = nullptr;   //  Pointers to command
  char const                    *_outListName    = nullptr;   //  line memory.
  char const                    *_outShowName    = nullptr;
  char const                    *_outPipeName    = nullptr;

  char const                    *_outStatsName   = nullptr;
  char const                    *_outHistoName   = nullptr;

  merylFileWriter               *_outDbse        = nullptr;   //  Master outputs.
  compressedFileWriter          *_outList        = nullptr;
  compressedFileWriter          *_outShow        = nullptr;
  //ar const                    *_outPipe        = nullptr;

  compressedFileWriter          *_outStats       = nullptr;
  compressedFileWriter          *_outHisto       = nullptr;

  //  deprecated
  bool _printACGTorder = false;



  //  Only so it can access _select when displaying the command tree.
  friend class merylCommandBuilder;
  friend class merylOpCompute;
};






inline
void
merylOpTemplate::addNewSelectorProduct(void) {
  if ((_select.size() > 0) &&        //  Add a new product if the existing one
      (_select.back().size() > 0))   //  is not empty.
    _select.emplace_back();
}

inline
void
merylOpTemplate::addSelectorToProduct(merylSelector const f) {

  _isSelector = true;

  if (_select.size() == 0)          //  If there are no selectors yet,
    _select.emplace_back();         //  add a new product term.

  _select.back().emplace_back(f);   //  Add 'f' to the last product term.
}

inline
bool
merylOpTemplate::isSelectorProductEmpty(void) {
  return ((_select.size() == 0) ||          //  Current selector product is empty if
          (_select.back().size() == 0));    //  it either doesn't exist or...is empty.
}

inline
merylSelector &
merylOpTemplate::getSelector(uint32 p, uint32 t) {
  return(_select[p][t]);
}

inline
merylSelector &
merylOpTemplate::getLastSelector(void) {
  return(_select.back().back());
}




#endif  //  MERYLOPTEMPLATE_H
