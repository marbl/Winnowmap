
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

#include "meryl.H"

//
//  Convert a template operation tree into a set of 64 compute operation
//  trees.
//


void
merylOpCompute::addSliceOutput(merylFileWriter *outDbse, uint32 slice) {
  if (outDbse) {
    _outDbseSlice = outDbse->getStreamWriter(slice);
    _outDbse      = _outDbseSlice;
  }
}


void
merylOpCompute::addSliceLister(opPname which, char const *name, compressedFileWriter *list, uint32 slice) {

  if (name == nullptr)                    //  If no name defined, no list output requested.
    return;

  if (list != nullptr) {                  //  If we're given an input cFW to use, use it.
    if (which == opPname::pnList)         //  This is either 'show' to stdout or a 'list'
      _outList = list;                    //  to a non-parallel single file.
    if (which == opPname::pnShow)         //
      _outShow = list;                    //  If list is nullptr, set up for a list output
    return;                               //  to parallel files below.
  }

  uint32      nm = strlen(name) + 2;      //  +1 for at most one extra digit (# -> ##), +1 for a NUL.
  uint32      ni = 0;                     //  We are the Knights Who Say "Ni!"
  uint32      nh = 0;
  char       *N  = new char [nm];
  char const *T  = name;

  for (; (*T != 0) && (*T != '#'); )      //  Copy the prefix to N.
    N[ni++] = *T++;                       //  Ni!

  for (; (*T == '#'); T++, nh++)          //  Count and skip the hashes.
    ;

  N = toDec(slice, N, std::max(nh, 2u));  //  Add the slice id, at least two digits.

  for (; (*T != 0); )                     //  Copy the suffix to N.
    N[ni++] = *T++;                       //  Ni!
  N[ni] = 0;                              //  Ni!  Ni!

  assert(ni < nm);                        //  Stop saying that!
  assert(which == opPname::pnList);

  _outListSlice = new compressedFileWriter(N);   //  This gets deleted in destructor.
  _outList      = _outListSlice;                 //  This does not get deleted.

  delete [] N;                            //  Ekke ekke ekke pitang zoo boing!
}


void
merylOpCompute::addSliceStatistics(merylOpTemplate *ot, uint32 slice) {
  if ((ot->_outStats != nullptr) ||
      (ot->_outHisto != nullptr))
    _statsAcc = new merylHistogram(1048576);
}



//
//  Clone the command tree(s) into thread-specific copies, one tree per thread.
//

void
merylCommandBuilder::spawnThreads(uint32 allowedThreads) {

  //  Allocate compute objects for each of our 64 file slices, then copy the
  //  list of templates into each slice.  These need to exist before we start
  //  creating inputs.

  for (uint32 ss=0; ss<64; ss++) {
    _thList[ss] = new merylOpCompute * [_opList.size()];

    for (uint32 oo=0; oo<_opList.size(); oo++)
      _thList[ss][oo] = new merylOpCompute(_opList[oo], ss, _opList[oo]->_inputs.size());
  }

  //  Save pointers to all the compute objects in each template.  This is
  //  used to collect statistics when we're done.  The arrays are annoyingly
  //  perpendicular to each other and so we need to save each of the 64
  //  pointers, instead of just pointing to an array of 64 elements.

  for (uint32 oo=0; oo<_opList.size(); oo++)
    for (uint32 ss=0; ss<64; ss++)
      _opList[oo]->_computes[ss] = _thList[ss][oo];

  //  Update all the input/output objects to be per-thread.

  for (uint32 ss=0; ss<64; ss++) {
    for (uint32 oo=0; oo<_opList.size(); oo++) {
      merylOpTemplate  *tpl = _opList[oo];       //  The template operation
      merylOpCompute   *cpu = _thList[ss][oo];   //  The per-thread operation we're creating.

      //  For each input to the template, add a new input to the compute.
      //
      for (uint32 ii=0; ii<tpl->_inputs.size(); ii++) {
        merylInput  *in = tpl->_inputs[ii];

        //  If the template input is from an action, we need to search for
        //  that action in the master list of templates, then set the
        //  per-thread input to be the corresponding object in the per-thread
        //  list.
        //

        if (in->isFromTemplate() == true) {
          uint32  inidx = UINT32_MAX;

          for (uint32 xx=0; xx<_opList.size(); xx++)
            if (in->_template == _opList[xx])
              inidx = xx;

          if (inidx == UINT32_MAX)
            fprintf(stderr, "Failed to find corresponding operation.\n"), exit(1);

          cpu->addSliceInput(new merylInput(_thList[ss][inidx]));
        }

        //  If the template input is from a database, make a new input for
        //  just the piece we're processing in this thread.
        //
        if (in->isFromDatabase() == true) {
          merylFileReader  *r = new merylFileReader(in->inputName(), ss);
          merylInput       *i = new merylInput(r);

          cpu->addSliceInput(i);
        }

        //  If the template input is from anything else, it's an error.
        //
        assert(in->isFromSequence() == false);
        assert(in->isFromStore()    == false);
      }

      //  Add a stats accumulator for this slice.
      //
      cpu->addSliceStatistics(tpl, ss);

      //  Add outputs for this slice.
      //
      cpu->addSliceOutput(tpl->_outDbse, ss);
      cpu->addSliceLister(opPname::pnList, tpl->_outListName, tpl->_outList, ss);
      cpu->addSliceLister(opPname::pnShow, tpl->_outShowName, tpl->_outShow, ss);
    }
  }
}
