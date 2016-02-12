#ifndef INPUTSETTINGS
#define INPUTSETTINGS

#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

class Settings {
  public:
 
  static const int ntrkBins = 42;
  double xtrkbins[ntrkBins+1] = { 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 3.2 , 4.0 , 4.8 , 5.6 , 6.4 , 7.2 , 9.6 , 12.0, 14.4,19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8,73.6,86.4,103.6,120.8,140,165,190,220,250,280,310,350,400};
  
  static const int nTriggers = 4;
  double triggerBins[nTriggers+1] = {0,60,80,100,1200};
  double triggerOverlapBins[nTriggers] = {0,60,80,100};
  static const int HInTriggers = 5;
  double HItriggerBins[HInTriggers+1] = {0,60,80,100,120,1200};
  double HItriggerOverlapBins[HInTriggers] = {0,60,80,100,120};
  
  static const int njetBins = 240;
  static const int maxJetBin = 1200;
  
  static const int nTriggers_trk = 6;
  double triggerBins_trk[nTriggers_trk+1] = {0,20,26,36,47,55,1000};
  double triggerOverlapBins_trk[nTriggers_trk] = {0,20,26,36,47,55};
  static const int HInTriggers_trk = 6;
  double HItriggerBins_trk[HInTriggers_trk+1] = {0,14,20,26,36,47,1000};
  double HItriggerOverlapBins_trk[HInTriggers_trk] = {0,14,20,26,36,47};
  
  static const int nTrktriggerBins = 600;
  static const int maxTrktriggerBin = 600;
  

  static const int nCentBins = 28;
  int lowCentBin[nCentBins] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,0,6,10,2,6,10,14,16};
  int highCentBin[nCentBins] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,20,20,20,6,10,14,16,18};

  TH2D *spec[nTriggers],               *HIspec[HInTriggers][nCentBins];
  TH1D *evtCount[nTriggers],           *HIevtCount[HInTriggers][nCentBins];
  TH1D *nVtxMB,                        *HInVtxMB;
  TH2D *spec_trk[nTriggers_trk],       *HIspec_trk[HInTriggers_trk][nCentBins];
  TH1D *evtCount_trk[nTriggers_trk],   *HIevtCount_trk[HInTriggers_trk][nCentBins];
  TH1D *nVtxMB_trk,                    *HInVtxMB_trk;

  TH1D * pp,                           *HI[nCentBins];
  TH1D * ppByTrigger[nTriggers],       *HIByTrigger[HInTriggers][nCentBins];
  TH1D * ppUsedByTrigger[nTriggers],   *HIUsedByTrigger[HInTriggers][nCentBins];
  TH1D * ppJets,                       *HIJets[nCentBins];
  TH1D * ppJetsByTrigger[nTriggers],   *HIJetsByTrigger[HInTriggers][nCentBins];
  TH1D * pp_perMBTrigger,              *HI_perMBTrigger[nCentBins];               

  TH1D * pp_trk,                              *HI_trk[nCentBins];                                  
  TH1D * ppByTrigger_trk[nTriggers_trk],      *HIByTrigger_trk[HInTriggers_trk][nCentBins];            
  TH1D * ppUsedByTrigger_trk[nTriggers_trk],  *HIUsedByTrigger_trk[HInTriggers_trk][nCentBins];        
  TH1D * ppMaxtrk,                            *HIMaxtrk[nCentBins];                              
  TH1D * ppMaxtrkByTrigger[nTriggers_trk],    *HIMaxtrkByTrigger[HInTriggers_trk][nCentBins];        
  TH1D * pp_perMBTrigger_trk,                 *HI_perMBTrigger_trk[nCentBins];                     

  
  Settings();
};
 
Settings::Settings()
{
  std::cout << "Getting setting.." << std::endl;
  return;
}

#endif
