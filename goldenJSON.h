#ifndef goldenJSON_h
#define goldenJSON_h

bool isGoodMB(int MBPD, const UInt_t run)
{
  bool isGood = false;
  if(MBPD==2 && run>=263155 && run<=263797) isGood=true;  
  if((MBPD==3 || MBPD==4) && run>=263192 && run<=263797) isGood=true;  
  return isGood;
}

Bool_t isInGoldenJSON(const UInt_t run, const UInt_t lumi)
{
  if(run == 262620){
    if(lumi >= 99 && lumi <= 186) return true;
    if(lumi >= 192 && lumi <= 423) return true;
  }
  else if(run == 262640){
    if(lumi >= 87 && lumi <= 102) return true;
    if(lumi >= 105 && lumi <= 232) return true;
  }
  else if(run == 262656){
    if(lumi >= 1 && lumi <= 101) return true;
    if(lumi >= 104 && lumi <= 133) return true;
    if(lumi >= 135 && lumi <= 148) return true;
    if(lumi >= 151 && lumi <= 178) return true;
  }
  else if(run == 262694){
    if(lumi >= 72 && lumi <= 225) return true;
  }
  else if(run == 262695){
    if(lumi >= 3 && lumi <= 223) return true;
  }
  else if(run == 262697){
    if(lumi >= 1 && lumi <= 62) return true;
  }
  else if(run == 262698){
    if(lumi >= 1 && lumi <= 11) return true;
    if(lumi >= 15 && lumi <= 15) return true;
  }
  else if(run == 262699){
    if(lumi >= 1 && lumi <= 34) return true;
    if(lumi >= 36 && lumi <= 40) return true;
  }
  else if(run == 262701){
    if(lumi >= 1 && lumi <= 1) return true;
  }
  else if(run == 262702){
    if(lumi >= 1 && lumi <= 1) return true;
  }
  else if(run == 262703){
    if(lumi >= 1 && lumi <= 123) return true;
  }
  else if(run == 262726){
    if(lumi >= 64 && lumi <= 500) return true;
  }
  else if(run == 262733){
    if(lumi >= 1 && lumi <= 6) return true;
    if(lumi >= 60 && lumi <= 83) return true;
  }
  else if(run == 262735){
    if(lumi >= 1 && lumi <= 206) return true;
  }
  else if(run == 262768){
    if(lumi >= 59 && lumi <= 425) return true;
  }
  else if(run == 262777){
    if(lumi >= 1 && lumi <= 144) return true;
  }
  else if(run == 262783){
    if(lumi >= 2 && lumi <= 5) return true;
  }
  else if(run == 262784){
    if(lumi >= 1 && lumi <= 374) return true;
  }
  else if(run == 262811){
    if(lumi >= 7 && lumi <= 438) return true;
  }
  else if(run == 262813){
    if(lumi >= 1 && lumi <= 4) return true;
    if(lumi >= 6 && lumi <= 6) return true;
  }
  else if(run == 262814){
    if(lumi >= 1 && lumi <= 2) return true;
    if(lumi >= 5 && lumi <= 6) return true;
  }
  else if(run == 262815){
    if(lumi >= 1 && lumi <= 4) return true;
    if(lumi >= 6 && lumi <= 6) return true;
  }
  else if(run == 262816){
    if(lumi >= 1 && lumi <= 449) return true;
  }
  else if(run == 262817){
    if(lumi >= 1 && lumi <= 27) return true;
  }
  else if(run == 262818){
    if(lumi >= 1 && lumi <= 75) return true;
    if(lumi >= 77 && lumi <= 279) return true;
  }
  else if(run == 262819){
    if(lumi >= 1 && lumi <= 118) return true;
  }
  else if(run == 262834){
    if(lumi >= 73 && lumi <= 214) return true;
  }
  else if(run == 262836){
    if(lumi >= 1 && lumi <= 99) return true;
  }
  else if(run == 262837){
    if(lumi >= 1 && lumi <= 777) return true;
  }
  else if(run == 262893){
    if(lumi >= 33 && lumi <= 1054) return true;
  }
  else if(run == 262921){
    if(lumi >= 9 && lumi <= 887) return true;
    if(lumi >= 889 && lumi <= 1072) return true;
  }
  else if(run == 262987){
    if(lumi >= 71 && lumi <= 156) return true;
    if(lumi >= 161 && lumi <= 174) return true;
  }
  else if(run == 262988){
    if(lumi >= 1 && lumi <= 689) return true;
  }
  else if(run == 263005){
    if(lumi >= 63 && lumi <= 654) return true;
  }
  else if(run == 263007){
    if(lumi >= 1 && lumi <= 128) return true;
  }
  else if(run == 263022){
    if(lumi >= 3 && lumi <= 13) return true;
    if(lumi >= 27 && lumi <= 51) return true;
    if(lumi >= 66 && lumi <= 204) return true;
  }
  else if(run == 263031){
    if(lumi >= 1 && lumi <= 110) return true;
  }
  else if(run == 263033){
    if(lumi >= 1 && lumi <= 62) return true;
  }
  else if(run == 263035){
    if(lumi >= 1 && lumi <= 609) return true;
  }
  else if(run == 263233){
    if(lumi >= 43 && lumi <= 427) return true;
    if(lumi >= 429 && lumi <= 437) return true;
    if(lumi >= 439 && lumi <= 486) return true;
    if(lumi >= 488 && lumi <= 735) return true;
  }
  else if(run == 263234){
    if(lumi >= 1 && lumi <= 695) return true;
  }
  else if(run == 263261){
    if(lumi >= 82 && lumi <= 1357) return true;
  }
  else if(run == 263284){
    if(lumi >= 71 && lumi <= 223) return true;
  }
  else if(run == 263286){
    if(lumi >= 1 && lumi <= 602) return true;
  }
  else if(run == 263293){
    if(lumi >= 1 && lumi <= 373) return true;
  }
  else if(run == 263322){
    if(lumi >= 64 && lumi <= 1238) return true;
  }
  else if(run == 263333){
    if(lumi >= 78 && lumi <= 227) return true;
  }
  else if(run == 263349){
    if(lumi >= 60 && lumi <= 401) return true;
  }
  else if(run == 263355){
    if(lumi >= 1 && lumi <= 117) return true;
  }
  else if(run == 263357){
    if(lumi >= 1 && lumi <= 2) return true;
  }
  else if(run == 263358){
    if(lumi >= 1 && lumi <= 3) return true;
  }
  else if(run == 263359){
    if(lumi >= 1 && lumi <= 2) return true;
  }
  else if(run == 263362){
    if(lumi >= 1 && lumi <= 574) return true;
  }
  else if(run == 263379){
    if(lumi >= 82 && lumi <= 916) return true;
  }
  else if(run == 263400){
    if(lumi >= 63 && lumi <= 904) return true;
  }
  else if(run == 263406){
    if(lumi >= 82 && lumi <= 99) return true;
  }
  else if(run == 263410){
    if(lumi >= 78 && lumi <= 664) return true;
  }
  else if(run == 263412){
    if(lumi >= 1 && lumi <= 240) return true;
  }
  else if(run == 263491){
    if(lumi >= 52 && lumi <= 410) return true;
  }
  else if(run == 263502){
    if(lumi >= 74 && lumi <= 413) return true;
    if(lumi >= 415 && lumi <= 497) return true;
    if(lumi >= 499 && lumi <= 512) return true;
  }
  else if(run == 263511){
    if(lumi >= 1 && lumi <= 43) return true;
    if(lumi >= 45 && lumi <= 101) return true;
    if(lumi >= 104 && lumi <= 116) return true;
    if(lumi >= 118 && lumi <= 118) return true;
    if(lumi >= 120 && lumi <= 126) return true;
    if(lumi >= 128 && lumi <= 133) return true;
    if(lumi >= 137 && lumi <= 143) return true;
    if(lumi >= 147 && lumi <= 147) return true;
    if(lumi >= 149 && lumi <= 237) return true;
  }
  else if(run == 263584){
    if(lumi >= 78 && lumi <= 289) return true;
    if(lumi >= 316 && lumi <= 959) return true;
  }
  else if(run == 263604){
    if(lumi >= 81 && lumi <= 978) return true;
    if(lumi >= 316 && lumi <= 959) return true;
  }
  else if(run == 263614){
    if(lumi >= 80 && lumi <= 503) return true;
    if(lumi >= 505 && lumi <= 725) return true;
    if(lumi >= 753 && lumi <= 1015) return true;
    if(lumi >= 1017 && lumi <= 1050) return true;
    if(lumi >= 1053 && lumi <= 1092) return true;
  }

  return false;
}

#endif
