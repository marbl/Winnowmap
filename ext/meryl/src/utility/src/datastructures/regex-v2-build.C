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

#include "regex-v2.H"
#include "arrays.H"

#include <set>
#include <map>
#include <algorithm>

namespace merylutil::inline regex::inline v2 {



void
regExState::addMatch(uint64 tid, regExToken tok, bool verbose) {

  if (_mid    != uint64max)  fprintf(stderr, "ERROR: can't add second match rule to state id=%lu.\n", _id);
  if (_lid[0] != uint64max)  fprintf(stderr, "ERROR: can't add match rule to state id=%lu; has lambda rules.\n", _id);

  assert(_mid    == uint64max);
  assert(_lid[0] == uint64max);

  if (verbose)
    fprintf(stderr, "  addMatch()   %-3lu match  -> %lu\n", _id, tid);

  _tok = tok;
  _mid = tid;
}

void
regExState::addEpsilon(uint64 tid, bool verbose) {

  if (_mid    != uint64max)  fprintf(stderr, "ERROR: can't add epsilon to state id=%lu; has existing _match rule.\n", _id);
  if (_lid[1] != uint64max)  fprintf(stderr, "ERROR: can't add third epsilon rule to state id=%lu.\n", _id);

  assert(_mid    == uint64max);
  assert(_lid[1] == uint64max);

  if (_lid[0] == uint64max)  _lid[0] = tid;
  else                       _lid[1] = tid;

  if (verbose)
    fprintf(stderr, "  addEpsilon() %-3lu lambda -> %lu\n", _id, tid);
}





//  Make sure there is enough space in the rs array to hold at least 'more'
//  more nodes.  Reinitialize the node id's if we grow space - we could
//  probably get away with only initializing from rsLen to (the new) rsMax.
//
void
expandrs(regExState *&rs, uint64 &rsLen, uint64 &rsMax, uint64 more) {
  if (rsLen + more < rsMax)
    return;

  increaseArray(rs, rsLen + more, rsMax, 1024);

  for (uint32 ii=0; ii<rsMax; ii++)
    rs[ii]._id = ii;
}


regExExpr
regEx::concat(regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExExpr a, regExExpr b) {

  if (vBuild)
    fprintf(stderr, "concat() %lu..%lu -> %lu..%lu:\n",
            a._bid, a._eid,
            b._bid, b._eid);

  a.end(rs)->addEpsilon(b._bid, vBuild);

  return regExExpr{._bid=a._bid, ._eid=b._eid};
}


regExExpr
regEx::alternate(regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExExpr a, regExExpr b) {

  expandrs(rs, rsLen, rsMax, 2);

  //regExState  *bgn = &rs[rsLen];   uint64 bid = rsLen++;
  //regExState  *end = &rs[rsLen];   uint64 eid = rsLen++;
  uint64 bid = rsLen++;
  uint64 eid = rsLen++;

  if (vBuild)
    fprintf(stderr, "alternate() %lu -> %lu..%lu | %lu..%lu -> %lu\n",
            bid,
            a._bid, a._eid,
            b._bid, b._eid,
            eid);

  rs[bid].addEpsilon(a._bid, vBuild);
  rs[bid].addEpsilon(b._bid, vBuild);

  a.end(rs)->addEpsilon(eid, vBuild);
  b.end(rs)->addEpsilon(eid, vBuild);

  return regExExpr{._bid=bid, ._eid=eid};
}



void
findReachable(uint64 aid, regExState *rs, std::set<uint64> &reachable) {

  if ((aid == uint64max) || (reachable.contains(aid) == true))
    return;

  reachable.insert(aid);

  if (rs[aid]._mid    != uint64max)   findReachable(rs[aid]._mid,    rs, reachable);
  if (rs[aid]._lid[0] != uint64max)   findReachable(rs[aid]._lid[0], rs, reachable);
  if (rs[aid]._lid[1] != uint64max)   findReachable(rs[aid]._lid[1], rs, reachable);
}


regExExpr
duplicate(regExExpr re,
          regExState *&rs, uint64 &rsLen, uint64 &rsMax,
          bool verbose) {

  std::set<uint64>         reachable;
  std::map<uint64,uint64>  idMap;

  if (verbose)
    fprintf(stderr, "  duplicate() %lu..%lu\n", re._bid, re._eid);

  //  Find the reachable nodes from the beginning of this expression
  findReachable(re._bid, rs, reachable);

  //  Build a map from original node ident to new node ident.
  //  (but do not otherwise initialize nodes)
  for (auto a : reachable) {
    idMap[a] = rsLen++;
    //if (verbose)
    //  fprintf(stderr, "duplicate()    map %lu -> %lu\n", a, idMap[a]);
  }

  expandrs(rs, rsLen, rsMax, 0);

  //  Duplicate nodes.
  for (auto a : reachable) {
    regExState *o = &rs[       a  ];
    regExState *d = &rs[ idMap[a] ];

    assert(      a  < rsLen);
    assert(idMap[a] < rsLen);

    assert(o != nullptr);
    assert(d != nullptr);

    //>_id        = o->_id
    d->_accepting = o->_accepting;
    d->_tok       = o->_tok;

    if (o->_mid    != uint64max)   d->_mid    = idMap[o->_mid   ];
    if (o->_lid[0] != uint64max)   d->_lid[0] = idMap[o->_lid[0]];
    if (o->_lid[1] != uint64max)   d->_lid[1] = idMap[o->_lid[1]];
  }

  regExExpr ne;
  ne._bid = idMap[re._bid];
  ne._eid = idMap[re._eid];

  if (verbose)
    fprintf(stderr, "duplicate() %lu..%lu into %lu..%lu\n",
            re._bid, re._eid,
            ne._bid, ne._eid);

  return ne;
}


//  Cases:
//    1) min=0, max=uint64max - *     - the classic closure
//    2) min=0, max=1         - ?     - zero or one
//    3) min=1, max=uint64max - +     - one or more
//    4) min>0, max<uint64max - {x,y} - general case
//
regExExpr
regEx::closure(regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExExpr a, uint64 min, uint64 max) {

  expandrs(rs, rsLen, rsMax, 2);

  /* regExState  *bgn = &rs[rsLen]; */  uint64 bid = rsLen++;
  /* regExState  *end = &rs[rsLen]; */  uint64 eid = rsLen++;

  uint64       nnn =       min;
  uint64       mmm = max - min;
  regExExpr   *nnndups = nullptr;
  regExExpr   *mmmdups = nullptr;

  if (vBuild)
    fprintf(stderr, "closure() %lu -> %lu..%lu -> %lu {%s,%s}\n",
            bid, a._bid, a._eid, eid,
            (min == uint64max) ? "inf" : toDec(min),
            (max == uint64max) ? "inf" : toDec(max));

  if ((min == 0) && (max == uint64max)) {   //  Basic Kleene star, no minimum, no maximum.
    rs[bid].addEpsilon(eid, vBuild);           //   - no matches, jump straight to the exi.
    rs[bid].addEpsilon(a._bid, vBuild);        //   - or enter 'a'

    a.end(rs)->addEpsilon(eid, vBuild);        //   - from end of 'a', jump to the exit.
    a.end(rs)->addEpsilon(a._bid, vBuild);     //   - or go back to the start of 'a'.
    goto finishclosure;
  }

  if ((min == 0) && (max == 1)) {           //  For zero or one, it's just like star,
    rs[bid].addEpsilon(eid, vBuild);           //  but without the 'go back' lambda.
    rs[bid].addEpsilon(a._bid, vBuild);

    a.end(rs)->addEpsilon(eid, vBuild);
    //end(rs)->addEpsilon(a._bid, vBuild);     //  Do not go back to the start!
    goto finishclosure;
  }

  //  Add 'nnn' copies of 'a', chained together.

  if (nnn > 0) {
    nnndups = new regExExpr [nnn];
    for (uint64 ii=0; ii<nnn; ii++)
      nnndups[ii] = duplicate(a, rs, rsLen, rsMax, vBuild);

    rs[bid].addEpsilon(nnndups[0]._bid, vBuild);                    //  Enter the chain of 'a' copies

    for (uint64 ii=0; ii<nnn-1; ii++) {
      nnndups[ii].end(rs)->addEpsilon(nnndups[ii+1]._bid, vBuild);  //  Leave 'a' to next copy.
      //ndups[ii].end(rs)->addEpsilon(eid, vBuild);                 //  But NO short-cut to accepting yet!
    }

    nnndups[nnn-1].end(rs)->addEpsilon(eid, vBuild);                //  Leave (last copy of) 'a' to accepting.
  }

  //  If 'max' is infinite, add 'a' itself as a normal closure.

  if (max == uint64max) {
    if (nnn > 0) {                                               //  Enter 'a' from either the last
      //ndups[nnn-1].end(rs)->addEpsilon(eid, vBuild);              //  in the chain of min's or from
      nnndups[nnn-1].end(rs)->addEpsilon(a._bid, vBuild);           //  the newly created bgn state.
    } else {                                                     //
      rs[bid].addEpsilon(eid, vBuild);                              //  If from the chain, we've already
      rs[bid].addEpsilon(a._bid, vBuild);                           //  got the short-cut to the end.
    }

    a.end(rs)->addEpsilon(eid, vBuild);                             //   Leave 'a' to accepting.
    a.end(rs)->addEpsilon(a._bid, vBuild);                          //   Leave 'a' back to the start of it.
  }

  //  Otherwise, add 'mmm' copies of a, chained together, but allowing any one
  //  to jump to the accepting state.

  else if (mmm > 0) {
    mmmdups = new regExExpr [mmm];
    for (uint64 ii=0; ii<mmm; ii++)
      mmmdups[ii] = duplicate(a, rs, rsLen, rsMax, vBuild);

    if (nnn > 0) {                                               //  Enter the chain of 'a' copies from
      //ndups[nnn-1].end(rs)->addEpsilon(eid, vBuild);              //  either the last in the chain of
      nnndups[nnn-1].end(rs)->addEpsilon(mmmdups[0]._bid, vBuild);  //  min's or from the newly created
    } else {                                                     //  bgn state.
      rs[bid].addEpsilon(eid, vBuild);                              //
      rs[bid].addEpsilon(mmmdups[0]._bid, vBuild);                  //  Like above, we've already got
    }                                                            //  the short-cut to the end.

    for (uint64 ii=0; ii<mmm-1; ii++) {
      mmmdups[ii].end(rs)->addEpsilon(mmmdups[ii+1]._bid, vBuild);  //  Leave 'a' to next copy.
      mmmdups[ii].end(rs)->addEpsilon(eid, vBuild);                 //  Leave 'a' to accepting.
    }

    mmmdups[mmm-1].end(rs)->addEpsilon(eid, vBuild);                //  Leave (last copy of) 'a' to accepting.
  }

  delete [] nnndups;   //  We're done with these expressions; the nodes
  delete [] mmmdups;   //  are saved elsewhere.

 finishclosure:
  return regExExpr{._bid=bid, ._eid=eid};
}


regExExpr
regEx::symbol(regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExToken tok) {

  expandrs(rs, rsLen, rsMax, 3);

  //regExState  *bgn = &rs[rsLen];   uint64 bid = rsLen++;
  //regExState  *mat = &rs[rsLen];   uint64 mid = rsLen++;
  //regExState  *end = &rs[rsLen];   uint64 eid = rsLen++;
  uint64 bid = rsLen++;
  uint64 mid = rsLen++;
  uint64 eid = rsLen++;

  if (vBuild)
    fprintf(stderr, "symbol() %lu -> %lu -> %lu for %s\n", bid, mid, eid, tok.display());

  rs[bid].addEpsilon(mid, vBuild);
  rs[mid].addMatch(eid, tok, vBuild);

  return regExExpr{._bid=bid, ._eid=eid};
}


regExExpr
regEx::epsilon(regExState *&rs, uint64 &rsLen, uint64 &rsMax, regExToken tok) {

  expandrs(rs, rsLen, rsMax, 2);

  uint64 bid = rsLen++;
  uint64 eid = rsLen++;

  if (vBuild)
    fprintf(stderr, "epsilon() %lu -> %lu for %s\n", bid, eid, tok.display());

  rs[bid].addEpsilon(eid, vBuild);

  return regExExpr{._bid=bid, ._eid=eid};
}



bool
regEx::build(void) {

  if (vBuild) {
    fprintf(stderr, "\n");
    fprintf(stderr, "BUILDING\n");
    fprintf(stderr, "\n");
  }

  //  Iterate over all states.  Push charater-matching states onto the stack
  //  as a simple expression, then pop off the required number of expressions
  //  when an operator shows up, merge them according to the operator, and
  //  push the new expression onto the stack.

  merylutil::stack<regExExpr>  st;

  for (uint64 tt=0; tt<tlLen; tt++) {
    regExExpr a;
    regExExpr b;

    switch (tl[tt]._type) {
      case regExTokenType::rtAlternation:
        a = st.pop();
        b = st.pop();
        st.push(alternate(rs, rsLen, rsMax, b, a));
        break;

      case regExTokenType::rtClosure:
        a = st.pop();
        st.push(closure(rs, rsLen, rsMax, a, tl[tt]._min, tl[tt]._max));
        break;

      case regExTokenType::rtConcat:
        a = st.pop();
        b = st.pop();  //  if this doesn't exist we need to add a lambda
        st.push(concat(rs, rsLen, rsMax, b, a));
        break;

      case regExTokenType::rtCharClass:
        st.push(symbol(rs, rsLen, rsMax, tl[tt]));
        break;

      case regExTokenType::rtEpsilon:
        st.push(epsilon(rs, rsLen, rsMax, tl[tt]));
        break;

      default:
        fprintf(stderr, "ERROR: unknown type: %s\n", tl[tt].display());
        assert(0);
        break;
    }
  }

  if (vBuild)
    fprintf(stderr, "stLen %lu rsLen %lu\n", st.depth(), rsLen);

  assert(st.depth() == 1);
  re = st.pop();

  re.end(rs)->_accepting = true;

#if 0
  for (uint64 ii=0; ii<rsLen; ii++)
    fprintf(stderr, "rs[%02lu] - id %02lu acc %d match %-4s lambda %4s:%-4s\n",
            ii,
            rs[ii]._id,
            rs[ii]._accepting,
            rs[ii]._mid    != uint64max ? toDec(rs[ii]._mid)    : "none",
            rs[ii]._lid[0] != uint64max ? toDec(rs[ii]._lid[0]) : "none",
            rs[ii]._lid[1] != uint64max ? toDec(rs[ii]._lid[1]) : "none");
#endif

  return true;
}



}  //  merylutil::regex::v2
