#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "getTrkCorr_simple.h"
#include "getTrkCorr_simple_trkTriggered.h"
#include "Settings.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void countTracks(std::vector<std::string> inputFiles, int jobNum, int isPP, bool isTest = false)
{
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
 
  Settings s; 
  TH2D *spec[s.nTriggers],             *HIspec[s.HInTriggers][20];
  TH1D *evtCount[s.nTriggers],         *HIevtCount[s.HInTriggers][20];
  TH1D *nVtxMB,                        *HInVtxMB;
  TH2D *spec_trk[s.nTriggers_trk],     *HIspec_trk[s.HInTriggers_trk][20];
  TH1D *evtCount_trk[s.nTriggers_trk], *HIevtCount_trk[s.HInTriggers_trk][20];
  TH1D *nVtxMB_trk,                    *HInVtxMB_trk;

  if(isPP){
    for(int i = 0; i<s.nTriggers; i++)
    {
      spec[i] = new TH2D(Form("spectrum_trigger%d",i),"",s.njetBins,0,s.maxJetBin,s.ntrkBins,s.xtrkbins);
      evtCount[i] = new TH1D(Form("evtCount%d",i),";max jet p_{T};N",s.njetBins,0,s.maxJetBin);
      evtCount[i]->SetMarkerColor(i);
    }
    nVtxMB = new TH1D("nVtxMB","nVtx;N Events",12,0,12);
    for(int i = 0; i<s.nTriggers_trk; i++)
    {
      spec_trk[i] = new TH2D(Form("spectrum_trigger%d_trk",i),"",s.nTrktriggerBins,0,s.maxTrktriggerBin,s.ntrkBins,s.xtrkbins);
      evtCount_trk[i] = new TH1D(Form("evtCount%d_trk",i),";max jet p_{T};N",s.nTrktriggerBins,0,s.maxTrktriggerBin);
      evtCount_trk[i]->SetMarkerColor(i);
    }
    nVtxMB_trk = new TH1D("nVtxMB_trk","nVtx;N Events",12,0,12);
  }else{  //end of pp loop, start of PbPb loop
    for(int i = 0; i<s.HInTriggers; i++)
    {
      for(int j = 0; j<20; j++){
        HIspec[i][j] = new TH2D(Form("HI_spectrum_trigger%d_cent%d",i,j),"",s.njetBins,0,s.maxJetBin,s.ntrkBins,s.xtrkbins);
        HIevtCount[i][j] = new TH1D(Form("HI_evtCount%d_cent%d",i,j),";max jet p_{T};N",s.njetBins,0,s.maxJetBin);
        HIevtCount[i][j]->SetMarkerColor(i);
      }
    }
    HInVtxMB = new TH1D("HI_nVtxMB","nVtx;N Events",12,0,12);
    for(int i = 0; i<s.HInTriggers_trk; i++)
    {
      for(int j = 0; j<20; j++){
        HIspec_trk[i][j] = new TH2D(Form("HI_spectrum_trigger%d_cent%d_trk",i,j),"",s.nTrktriggerBins,0,s.maxTrktriggerBin,s.ntrkBins,s.xtrkbins);
        HIevtCount_trk[i][j] = new TH1D(Form("HI_evtCount%d_cent%d_trk",i,j),";max jet p_{T};N",s.nTrktriggerBins,0,s.maxTrktriggerBin);
        HIevtCount_trk[i][j]->SetMarkerColor(i);
      }
    }
    HInVtxMB_trk = new TH1D("HI_nVtxMB_trk","nVtx;N Events",12,0,12);
  }   //end of PbPb loop
//******************************************************************************************************************************
//******************************************************************************************************************************
//******************************************************************************************************************************
  int nTrk;
  int nVtx;
  bool highPurity[60000];
  float trkPt[60000];
  float trkPtError[60000];
  float trkEta[60000];
  float trkPhi[60000];
  float trkMVA[60000];
  float trkDxy1[60000];
  float trkDxyError1[60000];
  float trkDz1[60000];
  float trkDzError1[60000];
  float pfEcal[60000];
  float pfHcal[60000];
  float trkChi2[60000];
  unsigned char trkNHit[60000];
  unsigned char trkNlayer[60000];
  unsigned char trkNdof[60000];
  unsigned char trkAlgo[60000];

  int pVtx;
  int pBeamScrape;
  //int NoiseFilter; 
  int pclusterCompatibilityFilter; 
  int pprimaryVertexFilter;  
  int phfCoincFilter3;
  int hiBin;

  int nref;
  float jtpt[200];
  float jteta[200];
  float rawpt[200];
  float chargedSum[200];
  float ecalSum[200];
  float hcalSum[200];

  int MB[20]={0};
  int j40=0;
  int j60=0;
  int j80=0;
  int t18=0;
  int t24=0;
  int t34=0;
  int t45=0;
  int t53=0;
  int HIMB[20]={0};
  int HIj40=0, HIj40_c30=0, HIj40_c50=0;
  int HIj60=0, HIj60_c30=0, HIj60_c50=0;
  int HIj80=0, HIj80_c30=0, HIj80_c50=0;
  int HIj100=0,HIj100_c30=0,HIj100_c50=0;
  int HIt12=0, HIt12_c10=0, HIt12_c30=0;
  int HIt18=0, HIt18_c10=0, HIt18_c30=0;
  int HIt24=0,              HIt24_c30=0;
  int HIt34=0,              HIt34_c30=0;
  int HIt45=0,              HIt45_c30=0;

  TrkCorr* trkCorr = new TrkCorr();
  TrkCorr_trkTriggered* trkCorr_trkTriggered = new TrkCorr_trkTriggered();
  TChain * trkCh;
  TChain * jetCh;
  TChain * evtCh;
  TChain * hltCh;
  TChain * hiCh;

  //for documenting which PD a file comes out of to avoid overlaps between PDs
  //0 is MB, 1 is jet40/60, 2 is jet80
  int PDindx[1000];
  for(unsigned int i = 0; i<inputFiles.size(); i++)
  {
    if(isPP){
      if(inputFiles.at(i).find("MinimumBias") != std::string::npos) PDindx[i]=0;
      if(inputFiles.at(i).find("HighPtLowerJets") != std::string::npos) PDindx[i]=1;
      if(inputFiles.at(i).find("HighPtJet80") != std::string::npos) PDindx[i]=2;
      if(inputFiles.at(i).find("FullTrack") != std::string::npos) PDindx[i]=3;
    }else{
      if(inputFiles.at(i).find("MinimumBias") != std::string::npos) PDindx[i]=0;
      if(inputFiles.at(i).find("HIHardProbes-") != std::string::npos) PDindx[i]=1;
      if(inputFiles.at(i).find("HIHardProbesPeripheral") != std::string::npos) PDindx[i]=2;
    }
  }

  if(isPP) trkCh = new TChain("ppTrack/trackTree");
  else     trkCh = new TChain("anaTrack/trackTree");
  
  for(unsigned int i = 0; i<inputFiles.size(); i++)  trkCh->Add(inputFiles.at(i).c_str());
  trkCh->SetBranchAddress("nTrk",&nTrk);
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  trkCh->SetBranchAddress("trkMVA",&trkMVA);
  trkCh->SetBranchAddress("trkNHit",&trkNHit);
  trkCh->SetBranchAddress("trkPtError",&trkPtError);
  trkCh->SetBranchAddress("pfHcal",&pfHcal);
  trkCh->SetBranchAddress("pfEcal",&pfEcal);
  trkCh->SetBranchAddress("trkDxy1",&trkDxy1);
  trkCh->SetBranchAddress("trkDxyError1",&trkDxyError1);
  trkCh->SetBranchAddress("trkDz1",&trkDz1);
  trkCh->SetBranchAddress("trkDzError1",&trkDzError1);
  trkCh->SetBranchAddress("nVtx",&nVtx);
  trkCh->SetBranchAddress("trkChi2",&trkChi2);
  trkCh->SetBranchAddress("trkNlayer",&trkNlayer);
  trkCh->SetBranchAddress("trkNdof",&trkNdof);
  trkCh->SetBranchAddress("trkAlgo",&trkAlgo);

  if(isPP) jetCh = new TChain("ak4CaloJetAnalyzer/t");
  else     jetCh = new TChain("akPu4CaloJetAnalyzer/t");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  jetCh->Add(inputFiles.at(i).c_str());
  jetCh->SetBranchAddress("nref",&nref);
  jetCh->SetBranchAddress("jtpt",&jtpt);
  jetCh->SetBranchAddress("jteta",&jteta);  
  jetCh->SetBranchAddress("rawpt",&rawpt);
  jetCh->SetBranchAddress("chargedSum",&chargedSum);  
  jetCh->SetBranchAddress("ecalSum",&ecalSum);
  jetCh->SetBranchAddress("hcalSum",&hcalSum);  
  trkCh->AddFriend(jetCh);

  evtCh = new TChain("skimanalysis/HltTree");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  evtCh->Add(inputFiles.at(i).c_str());
  if(isPP){
    evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
    evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
    //evtCh->SetBranchAddress("pHBHENoiseFilterResultProducer",&NoiseFilter);
  }else{
    evtCh->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);  
    evtCh->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);  
    evtCh->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3); 
  } 
  trkCh->AddFriend(evtCh);
 
  if(!isPP){
    hiCh = new TChain("hiEvtAnalyzer/HiTree");
    for(unsigned int i = 0; i<inputFiles.size(); i++)  hiCh->Add(inputFiles.at(i).c_str());
    hiCh->SetBranchAddress("hiBin",&hiBin);
    trkCh->AddFriend(hiCh);
  }
 
  hltCh = new TChain("hltanalysis/HltTree");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  hltCh->Add(inputFiles.at(i).c_str());
  if(isPP){
    for(int i = 0; i<20; i++) hltCh->SetBranchAddress(Form("HLT_L1MinimumBiasHF1OR_part%d_v1",i),&(MB[i]));
    hltCh->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&j40);
    hltCh->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&j60);
    hltCh->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&j80);
    hltCh->SetBranchAddress("HLT_FullTrack18ForPPRef_v3",&t18);
    hltCh->SetBranchAddress("HLT_FullTrack24ForPPRef_v3",&t24);
    hltCh->SetBranchAddress("HLT_FullTrack34ForPPRef_v4",&t34);
    hltCh->SetBranchAddress("HLT_FullTrack45ForPPRef_v3",&t45);
    hltCh->SetBranchAddress("HLT_FullTrack53ForPPRef_v3",&t53);
  }else{
    for(int i = 0; i<20; i++) hltCh->SetBranchAddress(Form("HLT_HIL1MinimumBiasHF2AND_part%d_v1",i),&(HIMB[i]));
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_v2",&HIj40);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_v1",&HIj60);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1",&HIj80);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1",&HIj100);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1",&HIj40_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1",&HIj60_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1",&HIj80_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent30_100_v1",&HIj100_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1",&HIj40_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1",&HIj60_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1",&HIj80_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_Cent50_100_v1",&HIj100_c50);
    hltCh->SetBranchAddress("HLT_HIFullTrack12_L1MinimumBiasHF2_AND_v1",&HIt12);
    hltCh->SetBranchAddress("HLT_HIFullTrack18_L1MinimumBiasHF2_AND_v1",&HIt18);
    hltCh->SetBranchAddress("HLT_HIFullTrack24_v1",&HIt24);
    hltCh->SetBranchAddress("HLT_HIFullTrack34_v1",&HIt34);
    hltCh->SetBranchAddress("HLT_HIFullTrack45_v1",&HIt45);
    //hltCh->SetBranchAddress("HLT_HIFullTrack12_L1Centrality010_v1",&HIt12_c10);
    //hltCh->SetBranchAddress("HLT_HIFullTrack18_L1Centrality010_v1",&HIt18_c10);
    hltCh->SetBranchAddress("HLT_HIFullTrack12_L1Centrality30100_v1",&HIt12_c30);
    hltCh->SetBranchAddress("HLT_HIFullTrack18_L1Centrality30100_v1",&HIt18_c30);
    hltCh->SetBranchAddress("HLT_HIFullTrack24_L1Centrality30100_v1",&HIt24_c30);
    hltCh->SetBranchAddress("HLT_HIFullTrack34_L1Centrality30100_v1",&HIt34_c30);
    hltCh->SetBranchAddress("HLT_HIFullTrack45_L1Centrality30100_v1",&HIt45_c30);
  }
  trkCh->AddFriend(hltCh);
//***********************************************************************************
//***********************************************************************
  std::cout << "starting event loop" << std::endl;
  std::cout << trkCh->GetEntries() << std::endl;
  for(int i = 0; i<trkCh->GetEntries(); i++)
  {
    //if(i%1000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<" "<<std::endl;
    if(i%10==0) std::cout << i<<"/"<<trkCh->GetEntries()<<" "<<std::endl;
    trkCh->GetEntry(i);
    //if(!NoiseFilter) continue;
    if(isPP && (!pVtx || !pBeamScrape)) continue;
    if(!isPP && (!pclusterCompatibilityFilter || !pprimaryVertexFilter || !phfCoincFilter3)) continue;
    bool MinBias = 0;
    for(int j = 0; j<20; j++) MinBias = MinBias || ((isPP)?(bool)MB[j]:(bool)HIMB[j]);
    if(isPP && !MinBias && !j40 && !j60 && !j80 && !t18 && !t24 && !t34 && !t45 && !t53) continue;
    //if(!isPP && !MinBias && !HIj40 && !HIj60 && !HIj80 && !HIj100 && !HIj40_c30 && !HIj60_c30 && !HIj80_c30 && !HIj100_c30&& !HIj40_c50 && !HIj60_c50 && !HIj80_c50 && !HIj100_c50 && !HIt12 && !HIt18 && !HIt24 && !HIt34 && !HIt45 && !HIt12_c10 && !HIt18_c10 && !HIt12_c30 && !HIt18_c30 && !HIt24_c30 && !HIt34_c30 && !HIt45_c30) continue;
    if(!isPP && !MinBias && !HIj40 && !HIj60 && !HIj80 && !HIj100 && !HIj40_c30 && !HIj60_c30 && !HIj80_c30 && !HIj100_c30&& !HIj40_c50 && !HIj60_c50 && !HIj80_c50 && !HIj100_c50 && !HIt12 && !HIt18 && !HIt24 && !HIt34 && !HIt45 && !HIt12_c30 && !HIt18_c30 && !HIt24_c30 && !HIt34_c30 && !HIt45_c30) continue;

    //**************************************************
    //for trigger combination with jet triggers
    float maxJtPt = 0;
    for(int j=0; j<nref; j++)
    {
      if(isPP && (chargedSum[j]/rawpt[j]<0.01 || TMath::Abs(jteta[j])>2)) continue;
      if(!isPP &&  (ecalSum[j]/(ecalSum[j]+hcalSum[j])<0.05 || hcalSum[j]/(ecalSum[j]+hcalSum[j])<0.1|| TMath::Abs(jteta[j])>2)) continue;
      if(jtpt[j]>maxJtPt) maxJtPt = jtpt[j];
    }

    float maxTrackPt = 0;
    for(int j=0; j<nTrk; j++)
    {
      if(TMath::Abs(trkEta[j])>1 || (int)trkAlgo[j]<4 || (int)trkAlgo[j]>8 || (int)trkNHit[j]<11 || trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>0.15 || !highPurity[j] || trkPtError[j]/trkPt[j]>0.1 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
      if(trkPt[j]>maxTrackPt) maxTrackPt = trkPt[j];
    }
    int PD = PDindx[trkCh->GetTreeNumber()];
    if(maxJtPt==0 && (PD==1 || PD==2)) continue;//remove jet events where no jets are in barrel  
    if(maxTrackPt==0 && ((isPP && PD==3) || (!isPP && PD==1))) continue;//remove jet events where no tracks are in barrel  
    if(MinBias && PD==0)
    {
      if(isPP){
        evtCount[0]->Fill(maxJtPt); 
        nVtxMB->Fill(nVtx);
        evtCount_trk[0]->Fill(maxTrackPt); 
        nVtxMB_trk->Fill(nVtx);
      }else{
        HIevtCount[0][hiBin/10]->Fill(maxJtPt);
        HInVtxMB->Fill(nVtx);
        HIevtCount_trk[0][hiBin/10]->Fill(maxTrackPt);
        HInVtxMB_trk->Fill(nVtx);
      }
    }
    if(j40 && PD==1) evtCount[1]->Fill(maxJtPt);  
    if(j60 && PD==1) evtCount[2]->Fill(maxJtPt);  
    if(j80 && PD==2) evtCount[3]->Fill(maxJtPt);  
    if(t18 && PD==3) evtCount_trk[1]->Fill(maxTrackPt);  
    if(t24 && PD==3) evtCount_trk[2]->Fill(maxTrackPt);  
    if(t34 && PD==3) evtCount_trk[3]->Fill(maxTrackPt);  
    if(t45 && PD==3) evtCount_trk[4]->Fill(maxTrackPt);  
    if(t53 && PD==3) evtCount_trk[5]->Fill(maxTrackPt);  
    if(HIj40 && PD==1)  HIevtCount[1][hiBin/10]->Fill(maxJtPt);  
    if(HIj60 && PD==1)  HIevtCount[2][hiBin/10]->Fill(maxJtPt);  
    if(HIj80 && PD==1)  HIevtCount[3][hiBin/10]->Fill(maxJtPt);  
    if(HIj100 && PD==1) HIevtCount[4][hiBin/10]->Fill(maxJtPt);  
    if(HIj40_c30  && !HIj40 && PD==2 && hiBin>=60)   HIevtCount[1][hiBin/10]->Fill(maxJtPt);  
    if(HIj60_c30  && !HIj60 && PD==2 && hiBin>=60)   HIevtCount[2][hiBin/10]->Fill(maxJtPt);  
    if(HIj80_c30  && !HIj80 && PD==2 && hiBin>=60)   HIevtCount[3][hiBin/10]->Fill(maxJtPt);  
    if(HIj100_c30 && !HIj100&& PD==2 && hiBin>=60)   HIevtCount[4][hiBin/10]->Fill(maxJtPt);  
    if(HIj40_c50  && !HIj40 && !HIj40_c30 && PD==2 && hiBin>=100)  HIevtCount[1][hiBin/10]->Fill(maxJtPt);  
    if(HIj60_c50  && !HIj60 && !HIj60_c30 && PD==2 && hiBin>=100)  HIevtCount[2][hiBin/10]->Fill(maxJtPt);  
    if(HIj80_c50  && !HIj80 && !HIj80_c30 && PD==2 && hiBin>=100)  HIevtCount[3][hiBin/10]->Fill(maxJtPt);  
    if(HIj100_c50 && !HIj100 && !HIj100_c30 && PD==2 && hiBin>=100)  HIevtCount[4][hiBin/10]->Fill(maxJtPt);  
    if(HIt12 && PD==1) HIevtCount_trk[1][hiBin/10]->Fill(maxTrackPt);  
    if(HIt18 && PD==1) HIevtCount_trk[2][hiBin/10]->Fill(maxTrackPt);  
    if(HIt24 && PD==1) HIevtCount_trk[3][hiBin/10]->Fill(maxTrackPt);  
    if(HIt34 && PD==1) HIevtCount_trk[4][hiBin/10]->Fill(maxTrackPt);  
    if(HIt45 && PD==1) HIevtCount_trk[5][hiBin/10]->Fill(maxTrackPt);  
    if(HIt12_c30 && !HIt12 && PD==1 && hiBin>=60) HIevtCount_trk[1][hiBin/10]->Fill(maxTrackPt);  
    if(HIt18_c30 && !HIt18 && PD==1 && hiBin>=60) HIevtCount_trk[2][hiBin/10]->Fill(maxTrackPt);  
    if(HIt24_c30 && !HIt34 && PD==1 && hiBin>=60) HIevtCount_trk[3][hiBin/10]->Fill(maxTrackPt);  
    if(HIt34_c30 && !HIt34 && PD==1 && hiBin>=60) HIevtCount_trk[4][hiBin/10]->Fill(maxTrackPt);  
    if(HIt45_c30 && !HIt45 && PD==1 && hiBin>=60) HIevtCount_trk[5][hiBin/10]->Fill(maxTrackPt);  
    //if(HIt12_c10 && !HIt12 && PD==1 && hiBin<20) HIevtCount_trk[1][hiBin/10]->Fill(maxTrackPt);  
    //if(HIt18_c10 && !HIt18 && PD==1 && hiBin<20) HIevtCount_trk[2][hiBin/10]->Fill(maxTrackPt);  
    
    if(PD!=3)
    {
      for(int j = 0; j<nTrk; j++)
      { 
        if(TMath::Abs(trkEta[j])>1) continue;
        if(trkPt[j]<0.5 || trkPt[j]>=400) continue;
        if(highPurity[j]!=1) continue;
        if((trkMVA[j]<0.5 && trkMVA[j]!=-99) || (int)trkNHit[j]<8 || trkPtError[j]/trkPt[j]>0.3 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
  
        //if((trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>15 && (trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>pfHcal[j]+pfEcal[j]) continue;} //Calo Matching 

        float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j]);
        //dividing by pt at bin center instead of track by track pt (just a convention)
        float binCenter;
        if(isPP) binCenter = spec[0]->GetYaxis()->GetBinCenter(spec[0]->GetYaxis()->FindBin(trkPt[j]));
        else     binCenter = HIspec[0][0]->GetYaxis()->GetBinCenter(HIspec[0][0]->GetYaxis()->FindBin(trkPt[j]));
        if(isPP){
          if(MinBias && PD==0) spec[0]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(j40 && PD==1)     spec[1]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(j60 && PD==1)     spec[2]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(j80 && PD==2)     spec[3]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
        }else{
          if(MinBias && PD==0) HIspec[0][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(HIj40 && PD==1)     HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(HIj60 && PD==1)     HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(HIj80 && PD==1)     HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(HIj100 && PD==1)    HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          if(HIj40_c30  && !HIj40 && PD==2 && hiBin>=60)   HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
          if(HIj60_c30  && !HIj60 && PD==2 && hiBin>=60)   HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
          if(HIj80_c30  && !HIj80 && PD==2 && hiBin>=60)   HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
          if(HIj100_c30 && !HIj100&& PD==2 && hiBin>=60)   HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
          if(HIj40_c50  && !HIj40 && !HIj40_c30 && PD==2 && hiBin>=100)   HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
          if(HIj60_c50  && !HIj60 && !HIj60_c30 && PD==2 && hiBin>=100)   HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
          if(HIj80_c50  && !HIj80 && !HIj80_c30 && PD==2 && hiBin>=100)   HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
          if(HIj100_c50 && !HIj100 && !HIj100_c30 && PD==2 && hiBin>=100) HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);  
        }
      } //end trk loop
    }//end if statement  
 
    if((!isPP && PD!=2) || (isPP && PD!=1 && PD!=2))
    { 
      for(int j = 0; j<nTrk; j++)
      {
        if(trkPt[j]<0.5 || trkPt[j]>=400) continue;
        if(TMath::Abs(trkEta[j])>1 || (int)trkAlgo[j]<4 || (int)trkAlgo[j]>8 || (int)trkNHit[j]<11 || trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>0.15 || !highPurity[j] || trkPtError[j]/trkPt[j]>0.1 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
        //if((trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>15 && (trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>pfHcal[j]+pfEcal[j]) continue;} //Calo Matching 

        float binCenter;
        if(isPP) binCenter = spec_trk[0]->GetYaxis()->GetBinCenter(spec[0]->GetYaxis()->FindBin(trkPt[j]));
        else     binCenter = HIspec_trk[0][0]->GetYaxis()->GetBinCenter(HIspec[0][0]->GetYaxis()->FindBin(trkPt[j]));
        float correction = trkCorr_trkTriggered->getTrkCorr(trkPt[j],trkEta[j]);
        //dividing by pt at bin center instead of track by track pt (just a convention)
        if(isPP){
          if(MinBias && PD==0) spec_trk[0]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(t18 && PD==3)     spec_trk[1]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(t24 && PD==3)     spec_trk[2]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(t34 && PD==3)     spec_trk[3]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(t45 && PD==3)     spec_trk[4]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(t53 && PD==3)     spec_trk[5]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
        }else{
          if(MinBias && PD==0) HIspec_trk[0][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(HIt12 && PD==1)     HIspec_trk[1][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(HIt18 && PD==1)     HIspec_trk[2][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(HIt24 && PD==1)     HIspec_trk[3][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(HIt34 && PD==1)     HIspec_trk[4][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          if(HIt45 && PD==1)     HIspec_trk[5][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          if(HIt12_c30 && !HIt12 && PD==1 && hiBin>=60) HIspec_trk[1][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          if(HIt18_c30 && !HIt18 && PD==1 && hiBin>=60) HIspec_trk[2][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          if(HIt24_c30 && !HIt34 && PD==1 && hiBin>=60) HIspec_trk[3][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          if(HIt34_c30 && !HIt34 && PD==1 && hiBin>=60) HIspec_trk[4][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          if(HIt45_c30 && !HIt45 && PD==1 && hiBin>=60) HIspec_trk[5][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          //if(HIt12_c10 && !HIt12 && PD==1 && hiBin<20) HIspec_trk[1][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
          //if(HIt18_c10 && !HIt18 && PD==1 && hiBin<20) HIspec_trk[2][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
        } 
      }//end trk loop
    }//end if statement
  }//end event loop

  //for pp
  TFile * outF = TFile::Open(Form("RAA_output_%d.root",jobNum),"recreate");
  outF->cd();
  if(isPP){
    for(int i = 0; i<s.nTriggers; i++)
    {
      spec[i]->Write();
      evtCount[i]->Write();
    }
    for(int i = 0; i<s.nTriggers_trk; i++)
    {
      spec_trk[i]->Write();
      evtCount_trk[i]->Write();
    }
    nVtxMB->Write();
    nVtxMB_trk->Write();
  }else{
    for(int i = 0; i<s.HInTriggers; i++)
    {
      for(int j = 0; j<20; j++){
        HIspec[i][j]->Write();
        HIevtCount[i][j]->Write();
      }
    }
    for(int i = 0; i<s.HInTriggers_trk; i++)
    {
      for(int j = 0; j<20; j++){
        HIspec_trk[i][j]->Write();
        HIevtCount_trk[i][j]->Write();
      }
    }
    HInVtxMB->Write();
    HInVtxMB_trk->Write();
    outF->Close(); 
  }
}



//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
int main(int argc, const char* argv[])
{
  if(argc != 5)
  {
    std::cout << "Usage: countTracks <fileList>  <job>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  int job = std::atoi(argv[2]);
  int totalJobs = std::atoi(argv[3]);
  int isPP = std::atoi(argv[4]);
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      if(line%totalJobs==job) listOfFiles.push_back(buffer);
      line++;
    }
  }
   
  countTracks(listOfFiles,job,isPP);
  return 0; 
}
