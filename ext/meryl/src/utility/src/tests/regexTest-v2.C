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

#include "types.H"
#include "strings.H"
#include "math.H"


class rtNode {                     //  A node in a tree representing a regex.
public:                            //  Used for building random regular expressions
  rtNode() {}                      //  and finding random strings that match.
  ~rtNode() {
    for (uint32 ii=0; ii<_lLen; ii++)
      delete _l[ii];
  }

  void     recursiveInit(rtNode  *parent,
                         merylutil::mtRandom &rn);
  char    *display(char *&str, uint64 &strLen, uint64 &strMax);

public:
  static
  uint32   _identC;                //  Global counter of ident.
  uint32   _ident     = 1;         //  Mostly unused ident of the node.
  rtNode  *_parent    = nullptr;   //  Pointer to parent node, also mostly unused.
  uint32   _depth     = 0;         //  Depth this node is at in the tree.

  bool     _isString  = false;     //  Mutually exclusive flags for node type.
  bool     _isConcat  = false;
  bool     _isAltern  = false;

  bool     _isCapture = false;     //  Annotations; is the operator a capture
  bool     _isClosure = false;     //  group?  A closure group with iteration
  uint32   _cMin  = 1;             //  limits {min,max}?
  uint32   _cMax  = 1;

  char     _chars[8] = { 0 };      //  String in a leaf node, or
  uint32   _lLen  = 0;             //  children in operator nodes.
  rtNode  *_l[16] = { nullptr };
};

uint32 rtNode::_identC = 1;


class rtLeaf {                     //  A leaf found during random walk
public:                            //  of the rtNode tree.
  char     _chars[8] = { 0 };      //    String the node emits
  uint32   _cgLen    = 0;          //    Capture groups the string is in
  uint32   _cg[32]   = {0};        //
};



void
rtNode::recursiveInit(rtNode  *parent,
                      merylutil::mtRandom &rn) {

  _ident  = _identC++;
  _parent = parent;
  _depth  = (parent == nullptr) ? (0) : (parent->_depth + 1);

  double isStringRN  = rn.mtRandomRealOpen();   //  Compute probabilistic events.
  double isConcatRN  = rn.mtRandomRealOpen();   //  Not the most concise method,
  double isClosureRN = rn.mtRandomRealOpen();   //  but easy to modify.
  double isCaptureRN = rn.mtRandomRealOpen();

  //  If a closure exists in a parent, decrase the probability.
  //
  uint32  pCmax = 1;
  for (rtNode *p=parent; p; p = p->_parent) {
    uint32  sib = 1;

    for (uint32 ii=0; ii<p->_lLen; ii++)   //  Sibling choices.
      sib += (p->_l[ii] != nullptr) ? (p->_l[ii]->_cMax - 1) : 0;

    pCmax += sib * p->_cMax;               //  Times parent.
  }

  if      (_depth == 0) {
    _isString  =   false;   //  Real boring regex if only a string at the root.
    _isConcat  =  (isConcatRN  < 0.9) && (_isString == false);
    _isAltern  = !(isConcatRN  < 0.9) && (_isString == false);
    _isClosure =  (isClosureRN < 0.2) && (_isString == false);
    _isCapture =  (isCaptureRN < 0.2) && (_isString == false);

    _cMin = (_isClosure) ? (1 + 2 * rn.mtRandomRealOpen()) : 1;  //  Small outer closure
    _cMax = (_isClosure) ? (2 + 3 * rn.mtRandomRealOpen()) : 1;  //  range; 1-4.
  }
  else if (_depth == 1) {
    _isString  =  (isStringRN  < 0.1);
    _isConcat  =  (isConcatRN  < 0.8) && (_isString == false);
    _isAltern  = !(isConcatRN  < 0.8) && (_isString == false);
    _isClosure =  (isClosureRN < 0.2) && (_isString == false);
    _isCapture =  (isCaptureRN < 0.2) && (_isString == false);

    _cMin = (_isClosure) ? (2 + 2 * rn.mtRandomRealOpen()) : 1;  //  Slightly larger;
    _cMax = (_isClosure) ? (3 + 3 * rn.mtRandomRealOpen()) : 1;  //  2-5.
  }
  else if (_depth == 2) {
    double pcc = (pCmax < 15) ? 0.2 : 0.0;

    _isString  =  (isStringRN  < 0.6);
    _isConcat  =  (isConcatRN  < 0.7) && (_isString == false);
    _isAltern  = !(isConcatRN  < 0.7) && (_isString == false);
    _isClosure =  (isClosureRN < pcc) && (_isString == false);
    _isCapture =  (isCaptureRN < 0.5) && (_isString == false);

    _cMin = (_isClosure) ? (1 + 4 * rn.mtRandomRealOpen()) : 1;  //  1-6
    _cMax = (_isClosure) ? (4 + 3 * rn.mtRandomRealOpen()) : 1;  //
  }
  else if (_depth == 3) {
    _isString  =   true;
    _isConcat  =   false;
    _isAltern  =   false;
    _isClosure =  (isClosureRN < 0.0) && (_isString == false);
    _isCapture =  (isCaptureRN < 0.1) && (_isString == false);

    _cMin = (_isClosure) ? (1 + 4 * rn.mtRandomRealOpen()) : 1;
    _cMax = (_isClosure) ? (4 + 4 * rn.mtRandomRealOpen()) : 1;
  }
  else {
    _isString  =   true;
    _isConcat  =   false;
    _isAltern  =   false;
    _isClosure =  (isClosureRN < 0.0) && (_isString == false);
    _isCapture =  (isCaptureRN < 0.1) && (_isString == false);

    _cMin = (_isClosure) ? (1 + 4 * rn.mtRandomRealOpen()) : 1;
    _cMax = (_isClosure) ? (4 + 4 * rn.mtRandomRealOpen()) : 1;
  }



  if      (_isString == true) {                                 //  Create the node:
    _chars[0] = _chars[1] = _chars[2] = _chars[3] = 0;          //  a short random string,
    _chars[4] = _chars[5] = _chars[6] = _chars[7] = 0;

    uint32 nc = (parent->_isAltern) ? (3 + rn.mtRandomRealOpen() * 4) :
                                      (1 + rn.mtRandomRealOpen() * 16 / parent->_lLen);

    nc = std::min(nc, 8u);

    for (uint32 ii=0; ii<nc; ii++)
      _chars[ii] = "acgt"[(uint32)floor(4 * rn.mtRandomRealOpen())];
    _chars[nc] = 0;
  }
  else {                                                        //  or 8-12 concat nodes,
    _lLen = ((_isConcat) ? 8 : 2) + 5 * rn.mtRandomRealOpen();  //  or 2-6 altern nodes.
    _lLen = std::min(_lLen, 16u);
  }

  //if (_isString)   fprintf(stderr, "[%04d]d=%d chars  - '%s'\n", _ident, _depth, _chars);
  //if (_isConcat)   fprintf(stderr, "[%04d]d=%d concat - %u children%s\n", _ident, _depth, _lLen, _isCapture ? " capture" : "");
  //if (_isAltern)   fprintf(stderr, "[%04d]d=%d altern - %u children%s\n", _ident, _depth, _lLen, _isCapture ? " capture" : "");

  for (uint32 ll=0; ll<_lLen; ll++) {                           //  Log the node we created,
    _l[ll] = new rtNode;                                        //  then create children.
    _l[ll]->recursiveInit(this, rn);
  }
}


char *
rtNode::display(char *&str, uint64 &strLen, uint64 &strMax) {
  merylutil::increaseArray(str, strLen + 4 + strlen(_chars) + _lLen + 1 + 1, strMax, 512, merylutil::_raAct::copyDataClearNew);
  //                                  '({c}' + <string> + '|' + ')' + NUL

  if      (_isString) {
    for (uint32 ii=0; _chars[ii]; )                //  Append the string
      str[strLen++] = _chars[ii++];
  }
  else {
    str[strLen++] = '(';                           //  Always append '('.
    str[strLen]   = '{';   strLen += _isCapture;   //  Append '{c}' and step
    str[strLen]   = 'c';   strLen += _isCapture;   //  over it if capture
    str[strLen]   = '}';   strLen += _isCapture;   //  is enabled.

    for (uint32 nn=0; nn<_lLen; nn++) {            //  Display each child,
      _l[nn]->display(str, strLen, strMax);        //

      if ((nn < _lLen-1) && (_isAltern))
        str[strLen++] = '|';
    }

    str[strLen++] = ')';                           //  Always append ')'.
  }

  if (_isClosure) {
    str[strLen++] = '{';
    toDec(_cMin, str, strLen, strMax);
    str[strLen++] = ',';
    toDec(_cMax, str, strLen, strMax);
    str[strLen++] = '}';
  }

  str[strLen] = 0;

  return str;
}


uint32  cgIdent=1;

void
dfs(rtNode *cur,
    merylutil::stack<rtLeaf>  &leafs,
    merylutil::stack<uint32>  &cg,
    merylutil::mtRandom       &rn) {
  uint32  iters  = cur->_cMin + (cur->_cMax - cur->_cMin) * rn.mtRandomRealOpen();

  if (cur->_isCapture)
    cg.push(cgIdent++);

  for (uint32 mm=0; mm<iters; mm++) {

    if      (cur->_isString) {
      rtLeaf  leaf;

      for (uint32 ii=0; ii<8; ii++)
        leaf._chars[ii] = cur->_chars[ii];

      for (uint32 ii=0; ii<cg.depth(); ii++)
        leaf._cg[leaf._cgLen++] = cg[ii];

      leafs.push(leaf);
    }

    else if (cur->_isConcat) {
      for (uint32 c=0; c < cur->_lLen; c++)
        dfs(cur->_l[c], leafs, cg, rn);
    }

    else if (cur->_isAltern) {
      uint32  a = cur->_lLen * rn.mtRandomRealOpen();

      dfs(cur->_l[a], leafs, cg, rn);
    }
    else {
      assert(0);
    }
  }

  if (cur->_isCapture)
    cg.pop();
}

class regexResult {
public:
  void     release(void) {
    for (uint32 ii=0; ii<cgLen; ii++)
      delete [] cg[ii];
    delete [] cg;
  }
  uint32   cgLen;
  char   **cg;
};

regexResult
iterate(rtNode *root, merylutil::mtRandom &rn) {
  merylutil::stack<rtLeaf>  leafs;
  merylutil::stack<uint32>  cg;

  cg.push(0);                             //  Add an initial capture group for the
  dfs(root, leafs.clear(), cg, rn);       //  whole string, then pick a random path.

  char **caps = new char * [cgIdent];  //  Allocate and initialize space for
  char **capi = new char * [cgIdent];  //  each capture.

  uint32 maxlen = 0;
  for (uint32 ii=0; ii<leafs.depth(); ii++)
    maxlen += strlen(leafs[ii]._chars);

  for (uint32 ii=0; ii<cgIdent; ii++) {
    caps[ii] = capi[ii] = new char [maxlen + 1];
    capi[ii][0] = 0;
  }

  //  Iterate over each leaf, copying the string to the capture groups
  //  it is included in.

  for (uint32 ii=0; ii<leafs.depth(); ii++) {
    for (uint32 jj=0; jj<leafs[ii]._cgLen; jj++) {
      char * l = leafs[ii]._chars;           //  Pointer to chars to copy.
      char *&c = capi[ leafs[ii]._cg[jj] ];  //  REFERENCE to array to copy to.

      while (*l)
        *c++ = *l++;
      *c = 0;
    }
  }

  //for (uint32 jj=0; jj<cgIdent; jj++)
  //  if (caps[jj][0])
  //    fprintf(stdout, "[%02d] %4u %s\n", jj, (uint32)(capi[jj] - caps[jj]), caps[jj]);

  //for (uint32 ii=0; ii<cgIdent; ii++)
  //  delete [] caps[ii];
  delete [] capi;
  //delete [] caps;

  regexResult r;
  r.cgLen = cgIdent;
  r.cg    = caps;

  return r;
}



int
main(int argc, char const **argv) {
  merylutil::mtRandom  rn;

  uint32             regex    = 0;
  uint32             strings  = 0;

  uint32             tests    = 0;
  uint32             iters    = 1;

  merylutil::regEx  re;

  if (argc < 3) {
    fprintf(stderr, "usage: %s [-v] <regEx-string> <text-string> ...\n", argv[0]);
    exit(1);
  }

  int               arg = 1;
  int               err = 0;
  for (arg=1; arg<argc; arg++) {
    if      (strcmp(argv[arg], "-v")    == 0)   re.enableVerbose();
    //se if (strcmp(argv[arg], "-r")    == 0)   tests = strtouint32(argv[++arg]);

    else if (strcmp(argv[arg], "-test") == 0)   tests = strtouint32(argv[++arg]);

    else if (strcmp(argv[arg], "-iter") == 0)   iters = strtouint32(argv[++arg]);
    else if (strcmp(argv[arg], "-seed") == 0)   rn.mtSetSeed(strtouint32(argv[++arg]));
    else
      break;
  }

  if (regex > 0) {
    fprintf(stdout, "COMPILE: '%s'\n", argv[regex]);
    re.compile(argv[regex]);

    for (uint32 ss=strings; ss<argc; ss++) {
      fprintf(stdout, "EVAL: '%s'\n", argv[ss]);
      if (re.match(argv[ss]))
        for (uint64 ii=0; ii<re.numCaptures(); ii++)
          fprintf(stdout, "  %2lu %s %3lu-%3lu '%s'\n", ii,
                  re.isValid(ii) ? "valid" : "inval", re.getBgn(ii), re.getEnd(ii), re.get(ii));
    }
  }

  for (uint32 test=0; test<tests; test++) {
    uint64  strMax = 0;
    uint64  strLen = 0;
    char   *str    = nullptr;

    rtNode *root   = new rtNode;

    fprintf(stdout, "INITIALIZE: seed %u iter %lu\n", rn.mtGetSeed(), rn.mtGetIterations());

    root->recursiveInit(nullptr, rn);
    root->display(str, strLen, strMax);

    fprintf(stdout, "COMPILE: %u '%s'\n", test, str);
    re.compile(str);

    for (uint32 ii=0; ii<iters; ii++) {
      fprintf(stdout, "EVALUTE: %u\n", iters);

      regexResult r = iterate(root, rn);

      fprintf(stdout, "EVALUTE: %s\n", r.cg[0]);

      if (re.match(r.cg[0])) {
        for (uint64 ii=0; ii<re.numCaptures(); ii++)
          fprintf(stdout, "  %2lu %s %3lu-%3lu '%s'\n", ii,
                  re.isValid(ii) ? "valid" : "inval", re.getBgn(ii), re.getEnd(ii), re.get(ii));
      }
      else {
        fprintf(stdout, "Fail.\n");
        exit(1);
      }

      r.release();
    }

    delete    root;
    delete [] str;
  }

  return 0;
}
