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

namespace merylutil::inline regex::inline v2 {


char const *
regExToken::display(char *str) {
  static
  char  dis[1024];
  char *out = dis;
  char *app = dis;

  if (str)
    out = app = str;

  switch (_type) {
    case regExTokenType::rtGroupBegin:
      app += sprintf(app, "group bgn  id:%2lu cap:%03lu grp:%03lu%s%s", _id, _capIdent, _grpIdent, _cap ? " CAP" : "", _pfx ? " PFX" : "");
      break;
    case regExTokenType::rtGroupEnd:
      app += sprintf(app, "group end  id:%2lu cap:%03lu grp:%03lu", _id, _capIdent, _grpIdent);
      break;
    case regExTokenType::rtAlternation:
      app += sprintf(app, "alternat   id:%2lu", _id);
      break;
    case regExTokenType::rtClosure:
      app += sprintf(app, "closure    id:%2lu {%s,%s}", _id, toDec(_min), (_max == uint64max) ? "inf" : toDec(_max));
      break;
    case regExTokenType::rtConcat:
      app += sprintf(app, "concat     id:%2lu", _id);
      break;
    case regExTokenType::rtLineStart:
      app += sprintf(app, "line-start id:%2lu", _id);
      break;
    case regExTokenType::rtLineEnd:
      app += sprintf(app, "line-end   id:%2lu", _id);
      break;
    case regExTokenType::rtCharClass:
      app += sprintf(app, "charact    id:%2lu", _id);

      if (_capGrp != nullptr) {
        app += sprintf(app, " (");
        for (uint32 ii=0; ii<_capGrp->size(); ii++)
          if (*(app-1) == '(')
            app += sprintf(app, "%u", _capGrp->get(ii));
          else
            app += sprintf(app, " %u", _capGrp->get(ii));
        app += sprintf(app, ")");
      }

      app += sprintf(app, " [");
      for (uint32 cc=0, ee=0; cc<256; cc=++ee) {   //  uint8 unsafe!
        if (isMatch(cc) == false)                    //  Skip cc if not valid.
          continue;

        while ((ee < 256) && (isMatch(ee) == true))  //  Find last valid ee.
          ee++;
        ee--;

        auto displ = [&](uint32 x) {  if (isVisible(x)) return sprintf(app, "%c",      x);
                                      else              return sprintf(app, "\\x%02x", x);  };

        //fprintf(stderr, "DISPL cc=0x%02x ee=0x%02x\n", cc, ee);

        if      (cc   == ee)  { app += displ(cc);                                          }
        else if (cc+1 == ee)  { app += displ(cc);                      app += displ(ee);   }
        else                  { app += displ(cc);  app += displ('-');  app += displ(ee);   }
      }

      app += sprintf(app, "]");
      break;
    case regExTokenType::rtEpsilon:
      app += sprintf(app, "epsilon    id:%2lu cap:%03lu ", _id, _capIdent);
      break;
    case regExTokenType::rtNone:
      app += sprintf(app, "terminal   id:%2lu", _id);
      break;
    default:
      app += sprintf(app, "DEFAULT    id:%2lu", _id);
      break;
  }

  return out;
}



}  //    merylutil::regex::v2
