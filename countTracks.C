#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TBranch.h"
#include "TMath.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "getTrkCorr.h"
#include "TrkSettings.h"
#include "goldenJSON.h"
#include "EventSelectionCorrector.C"
#include "Settings.h"
#include "chi2Corrector/Chi2Corrector_PbPb.C"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void countTracks(std::vector<std::string> inputFiles, int jobNum, int isPP, bool isTest = false)
{
  //TODO
  //better corrections
  //better trig combo (replace triggers)
  //more MB
  //fix 19 tracks?
  //figure out pp EPOS difference
  //pp alignement
  //cancelations?

  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  bool doOnly1Vertex = false;
  bool doevtSelCorrection = true;
  bool dotrkcorr = true;
  bool useTrkCorrEverywhere =false;
  float caloMatchValue = 0.5;
  float caloMatchStart = 20;
  float jetEtaSelection = 2;
  float jetTrackCutThreshhold = 50;
  float trkBufferSize = 0.0;
  bool removePbPbPU = false;
  bool doChi2Shift = false;
  //TF1 * chiF = new TF1("chiF","(x<=20)*(0.976-0.159*TMath::Log10(x))+(x>20)*0.768",0.1,500); 
  //std::cout << chiF->Eval(1) << " " << chiF->Eval(18) << " " << chiF->Eval(21) << " " << chiF->Eval(40)  << std::endl;

  TH1D * hiHFDist = new TH1D("hiHFDist",";hiHF;Counts",200,0,10000);
 
  Settings s; 
  if(isPP){
    for(int i = 0; i<s.nTriggers; i++)
    {
      s.spec[i] = new TH2D(Form("spectrum_trigger%d",i),"",s.njetBins,0,s.maxJetBin,s.ntrkBins,s.xtrkbins);
      s.evtCount[i] = new TH1D(Form("evtCount%d",i),";max jet p_{T};N",s.njetBins,0,s.maxJetBin);
      s.evtCount[i]->SetMarkerColor(i);
      s.evtCount_JetVars[i] = new TH2D(Form("evtCount_JetVars%d",i),"max jet #eta;max jet p_{T};N",10,-2,2,16,40,120);
    }
    s.nVtxMB = new TH1D("nVtxMB","nVtx;N Events",12,0,12);
    for(int i = 0; i<s.nTriggers_trk; i++)
    {
      s.spec_trk[i] = new TH2D(Form("spectrum_trigger%d_trk",i),"",s.nTrktriggerBins,0,s.maxTrktriggerBin,s.ntrkBins,s.xtrkbins);
      s.evtCount_trk[i] = new TH1D(Form("evtCount%d_trk",i),";max jet p_{T};N",s.nTrktriggerBins,0,s.maxTrktriggerBin);
      s.evtCount_trk[i]->SetMarkerColor(i);
    }
    s.nVtxMB_trk = new TH1D("nVtxMB_trk","nVtx;N Events",12,0,12);
  }else{  //end of pp loop, start of PbPb loop
    for(int i = 0; i<s.HInTriggers; i++)
    {
      for(int j = 0; j<20; j++){
        s.HIspec[i][j] = new TH2D(Form("HI_spectrum_trigger%d_cent%d",i,j),"",s.njetBins,0,s.maxJetBin,s.ntrkBins,s.xtrkbins);
        s.HIevtCount[i][j] = new TH1D(Form("HI_evtCount%d_cent%d",i,j),";max jet p_{T};N",s.njetBins,0,s.maxJetBin);
        s.HIevtCount[i][j]->SetMarkerColor(i);
        s.HIevtCount_JetVars[i][j] = new TH2D(Form("HI_evtCount_JetVars%d_cent%d",i,j),"max jet #eta;max jet p_{T};N",10,-2,2,16,40,120);
      }
    }
    for(int j = 0; j<20; j++) s.HInVtxMB[j] = new TH1D(Form("HI_nVtxMB_%d",j),"nVtx;N Events",12,0,12);
    for(int i = 0; i<s.HInTriggers_trk; i++)
    {
      for(int j = 0; j<20; j++){
        s.HIspec_trk[i][j] = new TH2D(Form("HI_spectrum_trigger%d_cent%d_trk",i,j),"",s.nTrktriggerBins,0,s.maxTrktriggerBin,s.ntrkBins,s.xtrkbins);
        s.HIevtCount_trk[i][j] = new TH1D(Form("HI_evtCount%d_cent%d_trk",i,j),";max jet p_{T};N",s.nTrktriggerBins,0,s.maxTrktriggerBin);
        s.HIevtCount_trk[i][j]->SetMarkerColor(i);
      }
    }
    for(int j = 0; j<20; j++) s.HInVtxMB_trk[j] = new TH1D(Form("HI_nVtxMB_%d_trk",j),"nVtx;N Events",12,0,12);
  }   //end of PbPb loop
//******************************************************************************************************************************
//******************************************************************************************************************************
//******************************************************************************************************************************
  int nTrk;
  int nVtx;
  int nTrkTimesnVtx;
  bool highPurity[50000];
  float trkPt[50000];
  float trkPtError[50000];
  float trkEta[50000];
  float trkPhi[50000];
  float trkMVA[50000];
  float trkDxy1[50000];
  float trkDxyError1[50000];
  float trkDz1[50000];
  float trkDzError1[50000];
  float trkDzOverDzError[500000];
  float trkDxyOverDxyError[500000];
  float pfEcal[50000];
  float pfHcal[50000];
  float trkChi2[50000];
  float zVtx[20];
  unsigned char trkNHit[50000];
  unsigned char trkNlayer[50000];
  unsigned char trkNdof[50000];
  unsigned char trkAlgo[50000];
  unsigned char trkOriginalAlgo[50000];
  int trkCharge[50000];

  unsigned int run=0;
  unsigned int lumi=0;
  unsigned long long evt=0;

  int pVtx;
  int pBeamScrape;
  //int NoiseFilter; 
  int pclusterCompatibilityFilter; 
  int pprimaryVertexFilter;  
  int phfCoincFilter3;
  int hiBin = 1;
  float hiHF = 0;

  int nref;
  float jtpt[200];
  float jteta[200];
  float jtphi[200];
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
  int HIMB[21]={0};
  int HIj40_v1=0, HIj40_v2=0;
  int HIj40=0, HIj40_c30=0, HIj40_c50=0;
  int HIj60=0, HIj60_c30=0, HIj60_c50=0;
  int HIj80=0, HIj80_c30=0, HIj80_c50=0;
  int HIj100=0,HIj100_c30=0,HIj100_c50=0;
  int HIt12=0,              HIt12_c30=0;
  int HIt18=0,              HIt18_c30=0;
  int HIt24=0,              HIt24_c30=0;
  int HIt34=0,              HIt34_c30=0;
  int HI_Muon_L2Mu20=0;
 
  TrkCorr* trkCorr;
  TrkCorr* trkCorr_trk;
  TrkCorr* trkCorr_loosepp;
  if(isPP){  
    //trkCorr = new TrkCorr("TrkCorr_Mar15_Iterative_pp/");
    //trkCorr = new TrkCorr("TrkCorr_Feb16_Iterative_pp/");
    //trkCorr_trk = new TrkCorr("TrkCorr_Mar4_Iterative_pp_TrkTrig/");
    trkCorr = new TrkCorr("TrkCorr_May6_Iterative_pp/");
    //trkCorr = new TrkCorr("TrkCorr_May12_Iterative_RpPbCuts/");
    trkCorr_trk = new TrkCorr("TrkCorr_Mar4_Iterative_pp_TrkTrig/");
    //trkCorr_trk = new TrkCorr("TrkCorr_May6_Iterative_pp/");
    //trkCorr_loosepp = new TrkCorr("TrkCorr_Feb16_Iterative_pp/");
  }else{
    trkCorr = new TrkCorr("TrkCorr_Jun7_Iterative_PbPb/");
    trkCorr_trk = new TrkCorr("TrkCorr_Jun7_Iterative_PbPb/");
    //trkCorr = new TrkCorr("TrkCorr_May6_Iterative_PbPb/");
    //trkCorr_trk = new TrkCorr("TrkCorr_May6_Iterative_PbPb/");
  }
  EventSelectionCorrector corrEvSel;
  Chi2Corrector_PbPb * chi2corr = new Chi2Corrector_PbPb();

  TFile * inputFile;
  TTree * trkCh;
  TTree * jetCh;
  TTree * evtCh;
  TTree * hltCh;
  TTree * hiCh;

  //Ntuple for looking at specific tracks
  std::string trkSkimVars;
  trkSkimVars=   "run:lumi:evt:isPP:trkPt:trkEta:trkPhi:hiBin:hiHF:rmin:correction:maxjtpt:maxjteta:maxjtphi:inConeJetPt:inConeJetEta:inConeJetPhi:passesJetID:ecalSum:hcalSum:maxTrackPt:PD:trkNHit:trkChi2:trkMVA:highPurity:trkPtError:trkDxy1:trkDxyError1:trkDz1:trkDzError1:pfEcal:pfHcal:trkNlayer:trkNdof:trkAlgo:trkOriginalAlgo:isMB:isj40:isj60:isj80:isj100:ist12:ist18:ist24:ist34:isMu20";
  TNtuple * trkSkim  = new TNtuple("trkSkim","",trkSkimVars.data()); 

  //for documenting which PD a file comes out of to avoid overlaps between PDs
  //0 is MB, 1 is jet40/60, 2 is jet80
  int PDindx[5000];
  int MBPDindx[5000] = {0};
  for(unsigned int i = 0; i<inputFiles.size(); i++)
  {
    if(isPP){
      if((inputFiles.at(i).find("MinimumBias") != std::string::npos) || (inputFiles.at(i).find("MinBias") != std::string::npos)) PDindx[i]=0;
      else if(inputFiles.at(i).find("HighPtLowerJets") != std::string::npos) PDindx[i]=1;
      else if(inputFiles.at(i).find("HighPtJet80") != std::string::npos) PDindx[i]=2;
      else if(inputFiles.at(i).find("FullTrack") != std::string::npos) PDindx[i]=3;
      else PDindx[i]=-1;
    }else{
      if((inputFiles.at(i).find("MinimumBias") != std::string::npos) || (inputFiles.at(i).find("MinBias") != std::string::npos)){
        PDindx[i]=0;
        if((inputFiles.at(i).find("MinimumBias2") != std::string::npos) || (inputFiles.at(i).find("MinBias2") != std::string::npos))  MBPDindx[i]=2;
        if((inputFiles.at(i).find("MinimumBias3") != std::string::npos) || (inputFiles.at(i).find("MinBias3") != std::string::npos))  MBPDindx[i]=3;
        if((inputFiles.at(i).find("MinimumBias4") != std::string::npos) || (inputFiles.at(i).find("MinBias4") != std::string::npos))  MBPDindx[i]=4;
      }
      else if(inputFiles.at(i).find("HIHardProbes-") != std::string::npos) PDindx[i]=1;
      else if(inputFiles.at(i).find("HIHardProbesPeripheral") != std::string::npos) PDindx[i]=2;
      else PDindx[i]=-1;
    }
  }

  for(int nFile = 0; nFile<inputFiles.size(); nFile++){
    inputFile = TFile::Open(inputFiles.at(nFile).c_str(),"read");
    if(isPP) trkCh = (TTree*)inputFile->Get("ppTrack/trackTree");
    else     trkCh = (TTree*)inputFile->Get("anaTrack/trackTree");
    
    trkCh->SetBranchAddress("nTrk",&nTrk);
    trkCh->SetBranchAddress("nVtx",&nVtx);
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
    trkCh->SetBranchAddress("trkChi2",&trkChi2);
    trkCh->SetBranchAddress("trkNlayer",&trkNlayer);
    trkCh->SetBranchAddress("trkNdof",&trkNdof);
    trkCh->SetBranchAddress("trkAlgo",&trkAlgo);
    trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo);
    trkCh->SetBranchAddress("zVtx",&zVtx);
    trkCh->SetBranchAddress("trkCharge",&trkCharge);
    if(isPP){
      trkCh->SetBranchAddress("nTrkTimesnVtx",&nTrkTimesnVtx);
      trkCh->SetBranchAddress("trkDzOverDzError",&trkDzOverDzError);
      trkCh->SetBranchAddress("trkDxyOverDxyError",&trkDxyOverDxyError); 
    }
  
    if(isPP) jetCh = (TTree*)inputFile->Get("ak4CaloJetAnalyzer/t");
    else     jetCh = (TTree*)inputFile->Get("akPu4CaloJetAnalyzer/t");
    jetCh->SetBranchAddress("nref",&nref);
    jetCh->SetBranchAddress("jtpt",&jtpt);
    jetCh->SetBranchAddress("jteta",&jteta);  
    jetCh->SetBranchAddress("jtphi",&jtphi);  
    jetCh->SetBranchAddress("rawpt",&rawpt);
    jetCh->SetBranchAddress("chargedSum",&chargedSum);  
    jetCh->SetBranchAddress("ecalSum",&ecalSum);
    jetCh->SetBranchAddress("hcalSum",&hcalSum);  
    trkCh->AddFriend(jetCh);
  
    evtCh = (TTree*)inputFile->Get("skimanalysis/HltTree");
    if(isPP){
      evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
      evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
      //evtCh->SetBranchAddress("pHBHENoiseFilterResultProducer",&NoiseFilter);
    }else{
      evtCh->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);  
      evtCh->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);  
      evtCh->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3); 
    } 
   
    if(!isPP){
      hiCh = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
      hiCh->SetBranchAddress("hiBin",&hiBin);
      hiCh->SetBranchAddress("hiHF",&hiHF);
      hiCh->SetBranchAddress("run",&run);
      hiCh->SetBranchAddress("lumi",&lumi);
      hiCh->SetBranchAddress("evt",&evt);
      evtCh->AddFriend(hiCh);
    }
   
    hltCh = (TTree*)inputFile->Get("hltanalysis/HltTree");
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
      hltCh->SetBranchAddress("HLT_HIL1MinimumBiasHF2AND_v1",&(HIMB[20]));
      hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_v1",&HIj40_v1);
      hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_v2",&HIj40_v2);
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
      //hltCh->SetBranchAddress("HLT_HIFullTrack12_L1Centrality30100_v1",&HIt12_c30);
      //hltCh->SetBranchAddress("HLT_HIFullTrack18_L1Centrality30100_v1",&HIt18_c30);
      hltCh->SetBranchAddress("HLT_HIFullTrack24_L1Centrality30100_v1",&HIt24_c30);
      hltCh->SetBranchAddress("HLT_HIFullTrack34_L1Centrality30100_v1",&HIt34_c30);
      hltCh->SetBranchAddress("HLT_HIL2Mu20_v1",&HI_Muon_L2Mu20);
    }
    trkCh->AddFriend(hltCh);
  //***********************************************************************************
  //***********************************************************************
    std::cout << "starting event loop" << std::endl;
    std::cout << trkCh->GetEntries() << std::endl;
    for(int i = 0; i<trkCh->GetEntries(); i++)
    {
      HIj40_v1=0; HIj40_v2=0; //zeroing these triggers as precaution against a rare bug seen in 7 TeV pp versioned triggers, will be reset on next GetEntry
      for(int trig = 0; trig<21; trig++) HIMB[trig]=0;//same as above

      //if(i%1000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<" "<<std::endl;
      evtCh->GetEntry(i);
      //if(!NoiseFilter) continue;
      if(isPP && (!pVtx || !pBeamScrape)) continue;
      if(!isPP && (!pclusterCompatibilityFilter || !pprimaryVertexFilter || !phfCoincFilter3 || !isInGoldenJSON(run,lumi))) continue;

      trkCh->GetEntry(i);
      if(!isPP && removePbPbPU==true && hiHF>5400) continue;

      if(doOnly1Vertex && nVtx!=1) continue;
      
      bool hasGoodVtx = 0; 
      for(int vtx = 0; vtx<nVtx; vtx++){
        if(TMath::Abs(zVtx[vtx])<15) hasGoodVtx = true;
      }
      if(hasGoodVtx==false) continue;

      bool MinBias = 0;
      for(int j = 0; j<21; j++) MinBias = MinBias || ((isPP)?(MB[j]==1):(HIMB[j]==1 && isGoodMB(MBPDindx[nFile],run)));
      HIj40 = (HIj40_v1==1) || (HIj40_v2==1);
      if(isPP && !MinBias && !j40 && !j60 && !j80 && !t18 && !t24 && !t34 && !t45 && !t53) continue;
      if(!isPP && !MinBias && !HIj40 && !HIj60 && !HIj80 && !HIj100 && !HIj40_c30 && !HIj60_c30 && !HIj80_c30 && !HIj100_c30&& !HIj40_c50 && !HIj60_c50 && !HIj80_c50 && !HIj100_c50 && !HIt12 && !HIt18 && !HIt24 && !HIt34 && !HIt12_c30 && !HIt18_c30 && !HIt24_c30 && !HIt34_c30) continue;
  
      //**************************************************
      //for trigger combination with jet triggers
      float maxJtPt = 0;
      float maxJtEta = -99;
      float maxJtPhi = -99;
      for(int j=0; j<nref; j++)
      {
        if((chargedSum[j]/rawpt[j]<0.01 || TMath::Abs(jteta[j])>jetEtaSelection)) continue;
        //if(!isPP &&  ((ecalSum[j]/(ecalSum[j]+hcalSum[j])<0.05) || (hcalSum[j]/(ecalSum[j]+hcalSum[j])<0.1)|| TMath::Abs(jteta[j])>jetEtaSelection)) continue;
        if(jtpt[j]>maxJtPt){
          maxJtPt = jtpt[j];
          maxJtEta = jteta[j];
          maxJtPhi = jtphi[j];
        }
      }//end maxJt
  
      float maxTrackPt = 0;
      float eventMultiplicity = 0;
      for(int j=0; j<nTrk; j++)
      {
        if(TMath::Abs(trkEta[j])>2.4) continue;
        if(!highPurity[j]) continue;
        if(trkPt[j]<0.5 || trkPt[j]>400) continue;
      
        if(trkPtError[j]/trkPt[j]>0.1) continue; 

        if(trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>0.15) continue;      
        if(trkNHit[j]<11 && trkPt[j]>0.7) continue; 
        if(isPP){
          bool isCompatibleWithVertex = false;
          for(int v = 0; v<nVtx; v++){
            if(TMath::Abs(zVtx[v])>15) continue;
            if(TMath::Abs(trkDxyOverDxyError[j*nVtx+v])<3 && TMath::Abs(trkDzOverDzError[j*nVtx+v])<3){
              isCompatibleWithVertex = true;
              break;
            }
          } 
          if(!isCompatibleWithVertex) continue;
        }else if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
  
        float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
        if(!(trkPt[j]<caloMatchStart || (Et>caloMatchValue*trkPt[j]))) continue; //Calo Matching
        //if((maxJtPt>jetTrackCutThreshhold && trkPt[j]>maxJtPt+trkBufferSize) || (maxJtPt<=jetTrackCutThreshhold && trkPt[j]>jetTrackCutThreshhold+trkBufferSize)) continue;//upper boundary on track pt
        eventMultiplicity++;         

        //applying some tighter cuts before doing the maxTrackPt calculaiton
        if(TMath::Abs(trkEta[j])>1) continue;
        if(trkNHit[j]<11  || (int)trkOriginalAlgo[j]<4 || (int)trkOriginalAlgo[j]>7) continue; //track trigger cuts
        if(trkPt[j]>maxTrackPt) maxTrackPt = trkPt[j];
      }//end maxTrck
 
      float evtSelCorrection = 1;
      if(doevtSelCorrection){
        if(isPP){
          evtSelCorrection = corrEvSel.getEventWeightFromData(eventMultiplicity,nVtx);
        }
        else evtSelCorrection = 1;
      }
 
      int PD = PDindx[nFile];
      if(MinBias==1 && PD==0)
      {
        if(isPP){
          s.evtCount[0]->Fill(maxJtPt); 
          s.evtCount_JetVars[0]->Fill(maxJtEta,maxJtPt);
          s.nVtxMB->Fill(nVtx);
          s.evtCount_trk[0]->Fill(maxTrackPt); 
          s.nVtxMB_trk->Fill(nVtx);
        }else{
          s.HIevtCount[0][hiBin/10]->Fill(maxJtPt);
          s.HIevtCount_JetVars[0][hiBin/10]->Fill(maxJtEta,maxJtPt);
          s.HInVtxMB[hiBin/10]->Fill(nVtx);
          s.HIevtCount_trk[0][hiBin/10]->Fill(maxTrackPt);
          s.HInVtxMB_trk[hiBin/10]->Fill(nVtx);
          hiHFDist->Fill(hiHF);
        }
      }
      if(j40 && PD==1){s.evtCount_JetVars[1]->Fill(maxJtEta,maxJtPt); s.evtCount[1]->Fill(maxJtPt);}  
      if(j60 && PD==1){s.evtCount_JetVars[2]->Fill(maxJtEta,maxJtPt); s.evtCount[2]->Fill(maxJtPt);}  
      if(j80 && PD==2){s.evtCount_JetVars[3]->Fill(maxJtEta,maxJtPt); s.evtCount[3]->Fill(maxJtPt);} 
      if(t18 && PD==3) s.evtCount_trk[1]->Fill(maxTrackPt);  
      if(t24 && PD==3) s.evtCount_trk[2]->Fill(maxTrackPt);  
      if(t34 && PD==3) s.evtCount_trk[3]->Fill(maxTrackPt);  
      if(t45 && PD==3) s.evtCount_trk[4]->Fill(maxTrackPt);  
      if(t53 && PD==3) s.evtCount_trk[5]->Fill(maxTrackPt);  
      if(HIj40==1 && PD==1){ s.HIevtCount_JetVars[1][hiBin/10]->Fill(maxJtEta,maxJtPt);  s.HIevtCount[1][hiBin/10]->Fill(maxJtPt);} 
      if(HIj60 && PD==1){ s.HIevtCount_JetVars[2][hiBin/10]->Fill(maxJtEta,maxJtPt);  s.HIevtCount[2][hiBin/10]->Fill(maxJtPt);} 
      if(HIj80 && PD==1){ s.HIevtCount_JetVars[3][hiBin/10]->Fill(maxJtEta,maxJtPt);  s.HIevtCount[3][hiBin/10]->Fill(maxJtPt);} 
      if(HIj100 && PD==1){ s.HIevtCount_JetVars[4][hiBin/10]->Fill(maxJtEta,maxJtPt); s.HIevtCount[4][hiBin/10]->Fill(maxJtPt);} 
      if(HIj40_c30  && HIj40!=1 && PD==2 && hiBin>=60)   s.HIevtCount[1][hiBin/10]->Fill(maxJtPt);  
      if(HIj60_c30  && !HIj60 && PD==2 && hiBin>=60)   s.HIevtCount[2][hiBin/10]->Fill(maxJtPt);  
      if(HIj80_c30  && !HIj80 && PD==2 && hiBin>=60)   s.HIevtCount[3][hiBin/10]->Fill(maxJtPt);  
      if(HIj100_c30 && !HIj100&& PD==2 && hiBin>=60)   s.HIevtCount[4][hiBin/10]->Fill(maxJtPt);  
      if(HIj40_c50  && HIj40!=1 && !HIj40_c30 && PD==2 && hiBin>=100)  s.HIevtCount[1][hiBin/10]->Fill(maxJtPt);  
      if(HIj60_c50  && !HIj60 && !HIj60_c30 && PD==2 && hiBin>=100)  s.HIevtCount[2][hiBin/10]->Fill(maxJtPt);  
      if(HIj80_c50  && !HIj80 && !HIj80_c30 && PD==2 && hiBin>=100)  s.HIevtCount[3][hiBin/10]->Fill(maxJtPt);  
      if(HIj100_c50 && !HIj100 && !HIj100_c30 && PD==2 && hiBin>=100)  s.HIevtCount[4][hiBin/10]->Fill(maxJtPt);  
      if(HIt12 && PD==1) s.HIevtCount_trk[1][hiBin/10]->Fill(maxTrackPt);  
      if(HIt18 && PD==1) s.HIevtCount_trk[2][hiBin/10]->Fill(maxTrackPt);  
      if(HIt24 && PD==1) s.HIevtCount_trk[3][hiBin/10]->Fill(maxTrackPt);  
      if(HIt34 && PD==1) s.HIevtCount_trk[4][hiBin/10]->Fill(maxTrackPt);  
      //if(HIt12_c30 && !HIt12 && PD==1 && hiBin>=60) s.HIevtCount_trk[1][hiBin/10]->Fill(maxTrackPt);  
      //if(HIt18_c30 && !HIt18 && PD==1 && hiBin>=60) s.HIevtCount_trk[2][hiBin/10]->Fill(maxTrackPt);  
      if(HIt24_c30 && !HIt24 && PD==1 && hiBin>=60) s.HIevtCount_trk[3][hiBin/10]->Fill(maxTrackPt);  
      if(HIt34_c30 && !HIt34 && PD==1 && hiBin>=60) s.HIevtCount_trk[4][hiBin/10]->Fill(maxTrackPt);  
      
      if(PD!=3)
      {
        for(int j = 0; j<nTrk; j++)
        { 
          if(trkPt[j]<0.5 || trkPt[j]>=400) continue;
          if(TMath::Abs(trkEta[j])>1) continue;
          if(highPurity[j]!=1) continue;
          /*//FIXME
          float correction = 1;
          float rmin = 999;
          float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
          if(!(trkPt[j]<caloMatchStart || (Et>caloMatchValue*trkPt[j]))){
            if(!isPP){
               float inConeJetPt = 0;
               float inConeJetEta = 0;
               float inConeJetPhi = 0;
               float passesJetID = 0;
               float ecSum = 0;
               float hcSum = 0;
               for(int jt = 0; jt<nref; jt++)
               {
                 if((chargedSum[jt]/rawpt[jt]<0.01 || TMath::Abs(jteta[jt])>2)) continue;
                 //if((ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05 || hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1 || TMath::Abs(jteta[jt])>jetEtaSelection)) continue;
                 if(jtpt[jt]>inConeJetPt && TMath::Power(TMath::Power(trkEta[j]-jteta[jt],2)+TMath::Power(TMath::ACos(TMath::Cos(trkPhi[j]-jtphi[jt])),2),0.5)<0.4){
                   inConeJetPt =  jtpt[jt];
                   inConeJetEta = jteta[jt];
                   inConeJetPhi = jtphi[jt];
                   ecSum = ecalSum[jt];
                   hcSum = hcalSum[jt];
                   passesJetID = !((ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05) || (hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1));
                 }
               }
             float skimEntry[] = {(float)isPP,trkPt[j],trkEta[j],trkPhi[j],(float)hiBin,hiHF,rmin,correction,maxJtPt,maxJtEta,maxJtPhi,inConeJetPt,inConeJetEta,inConeJetPhi,passesJetID,ecSum,hcSum,maxTrackPt,(float)PD,(float)trkNHit[j],trkChi2[j],trkMVA[j],(float)highPurity[j],trkPtError[j],trkDxy1[j],trkDxyError1[j],trkDz1[j],trkDzError1[j],pfEcal[j],pfHcal[j],(float)trkNlayer[j],trkNdof[j],(float)trkAlgo[j],(float)trkOriginalAlgo[j],(float)MinBias,(float)HIj40,(float)HIj60,(float)HIj80,(float)HIj100,(float)HIt12,(float)HIt18,(float)HIt24,(float)HIt34,(float)HI_Muon_L2Mu20};
              trkSkim->Fill(skimEntry);
            }
            continue; //Calo Matching
          }
          //End FIXME*/

          if( trkPtError[j]/trkPt[j]>0.3) continue;       
         
          if(isPP){
            bool isCompatibleWithVertex = false;
            for(int v = 0; v<nVtx; v++){
              if(TMath::Abs(zVtx[v])>15) continue;
              if(TMath::Abs(trkDxyOverDxyError[j*nVtx+v])<3 && TMath::Abs(trkDzOverDzError[j*nVtx+v])<3){
                isCompatibleWithVertex = true;
                break;
              }
            }          
            if(!isCompatibleWithVertex) continue;
          }else{
            if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
            if(trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>((doChi2Shift)?chi2corr->getChi2Scale(hiBin,trkPt[j]):1)*0.15) continue; 
            if(trkPtError[j]/trkPt[j]>0.1) continue;       
            if(trkNHit[j]<11 && trkPt[j]>0.7) continue; 
          } 
   
          float rmin=999;
          for(int jt=0; jt<nref; jt++)
          {
            if((chargedSum[jt]/rawpt[jt]<0.01 || TMath::Abs(jteta[jt])>2)) continue;
            //if(!isPP &&  (ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05 || hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1|| TMath::Abs(jteta[jt])>2)) continue;
            if(jtpt[jt]<50) continue;
            float R = TMath::Power(jteta[jt]-trkEta[j],2) + TMath::Power(TMath::ACos(TMath::Cos(jtphi[jt]-trkPhi[j])),2);
            if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
          }
  
          float correction;
          if(!useTrkCorrEverywhere) correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin);
                                    //correction = trkCorr_loosepp->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin);
          else                      correction = trkCorr_trk->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin);
          correction = correction*evtSelCorrection;
 
          float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
          if(!(trkPt[j]<caloMatchStart || (Et>caloMatchValue*trkPt[j]))){
            /*if(!isPP){
               float inConeJetPt = 0;
               float inConeJetEta = 0;
               float inConeJetPhi = 0;
               float passesJetID = 0;
               float ecSum = 0;
               float hcSum = 0;
               for(int jt = 0; jt<nref; jt++)
               {
                 if((chargedSum[jt]/rawpt[jt]<0.01 || TMath::Abs(jteta[jt])>2)) continue;
                 //if((ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05 || hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1 || TMath::Abs(jteta[jt])>jetEtaSelection)) continue;
                 if(jtpt[jt]>inConeJetPt && TMath::Power(TMath::Power(trkEta[j]-jteta[jt],2)+TMath::Power(TMath::ACos(TMath::Cos(trkPhi[j]-jtphi[jt])),2),0.5)<0.4){
                   inConeJetPt =  jtpt[jt];
                   inConeJetEta = jteta[jt];
                   inConeJetPhi = jtphi[jt];
                   ecSum = ecalSum[jt];
                   hcSum = hcalSum[jt];
                   passesJetID = !((ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05) || (hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1));
                 }
               }
             float skimEntry[] = {(float)isPP,trkPt[j],trkEta[j],trkPhi[j],(float)hiBin,hiHF,rmin,correction,maxJtPt,maxJtEta,maxJtPhi,inConeJetPt,inConeJetEta,inConeJetPhi,passesJetID,ecSum,hcSum,maxTrackPt,(float)PD,(float)trkNHit[j],trkChi2[j],trkMVA[j],(float)highPurity[j],trkPtError[j],trkDxy1[j],trkDxyError1[j],trkDz1[j],trkDzError1[j],pfEcal[j],pfHcal[j],(float)trkNlayer[j],trkNdof[j],(float)trkAlgo[j],(float)trkOriginalAlgo[j],(float)MinBias,(float)HIj40,(float)HIj60,(float)HIj80,(float)HIj100,(float)HIt12,(float)HIt18,(float)HIt24,(float)HIt34,(float)HI_Muon_L2Mu20};
              trkSkim->Fill(skimEntry);
            }*/

            continue; //Calo Matching
          }

          if(useTrkCorrEverywhere && (trkNHit[j]<11 ||  (int)trkOriginalAlgo[j]<4 || (int)trkOriginalAlgo[j]>7)){
            /*if(!isPP){
              float skimEntry[] = {trkPt[j],trkEta[j],trkPhi[j],(float)hiBin,hiHF,rmin,correction,maxJtPt,maxTrackPt,(float)PD,(float)trkNHit[j],trkChi2[j],trkMVA[j],(float)highPurity[j],trkPtError[j],trkDxy1[j],trkDxyError1[j],trkDz1[j],trkDzError1[j],pfEcal[j],pfHcal[j],(float)trkNlayer[j],trkNdof[j],(float)trkAlgo[j],(float)trkOriginalAlgo[j],(float)MinBias,(float)HIj40,(float)HIj60,(float)HIj80,(float)HIj100,(float)HIt12,(float)HIt18,(float)HIt24,(float)HIt34};
              trkSkim->Fill(skimEntry);   
             }*///code for skimming tracks failing track cuts
            continue;
          }  

          if((maxJtPt>jetTrackCutThreshhold && trkPt[j]>maxJtPt+trkBufferSize) || (maxJtPt<=jetTrackCutThreshhold && trkPt[j]>jetTrackCutThreshhold+trkBufferSize)){
            if(!isPP){
               float inConeJetPt = 0;
               float inConeJetEta = 0;
               float inConeJetPhi = 0;
               float passesJetID = 0;
               float ecSum = 0;
               float hcSum = 0;
               for(int jt = 0; jt<nref; jt++)
               {
                 if((chargedSum[jt]/rawpt[jt]<0.01 || TMath::Abs(jteta[jt])>2)) continue;
                 //if((ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05 || hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1 || TMath::Abs(jteta[jt])>jetEtaSelection)) continue;
                 if(jtpt[jt]>inConeJetPt && TMath::Power(TMath::Power(trkEta[j]-jteta[jt],2)+TMath::Power(TMath::ACos(TMath::Cos(trkPhi[j]-jtphi[jt])),2),0.5)<0.4){
                   inConeJetPt =  jtpt[jt];
                   inConeJetEta = jteta[jt];
                   inConeJetPhi = jtphi[jt];
                   ecSum = ecalSum[jt];
                   hcSum = hcalSum[jt];
                   passesJetID = !((ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05) || (hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1));
                 }
               }
             float skimEntry[] = {(float)run,(float)lumi,(float)evt,(float)isPP,trkPt[j],trkEta[j],trkPhi[j],(float)hiBin,hiHF,rmin,correction,maxJtPt,maxJtEta,maxJtPhi,inConeJetPt,inConeJetEta,inConeJetPhi,passesJetID,ecSum,hcSum,maxTrackPt,(float)PD,(float)trkNHit[j],trkChi2[j],trkMVA[j],(float)highPurity[j],trkPtError[j],trkDxy1[j],trkDxyError1[j],trkDz1[j],trkDzError1[j],pfEcal[j],pfHcal[j],(float)trkNlayer[j],trkNdof[j],(float)trkAlgo[j],(float)trkOriginalAlgo[j],(float)MinBias,(float)HIj40,(float)HIj60,(float)HIj80,(float)HIj100,(float)HIt12,(float)HIt18,(float)HIt24,(float)HIt34,(float)HI_Muon_L2Mu20};
              trkSkim->Fill(skimEntry);
              //if(inConeJetPt==0) continue;//upper boundary on track pt
            }
          }

          //dividing by pt at bin center instead of track by track pt (just a convention)
          float binCenter;
          if(isPP) binCenter = s.spec[0]->GetYaxis()->GetBinCenter(s.spec[0]->GetYaxis()->FindBin(trkPt[j]));
          else     binCenter = s.HIspec[0][0]->GetYaxis()->GetBinCenter(s.HIspec[0][0]->GetYaxis()->FindBin(trkPt[j]));
          
          if(!dotrkcorr) correction=1;
          if(isPP){
            if(MinBias==1 && PD==0) s.spec[0]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(j40 && PD==1)     s.spec[1]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(j60 && PD==1)     s.spec[2]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(j80 && PD==2)     s.spec[3]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
          }else{
            if(MinBias==1 && PD==0) s.HIspec[0][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(HIj40==1 && PD==1)     s.HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(HIj60 && PD==1)     s.HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(HIj80 && PD==1)     s.HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(HIj100 && PD==1)    s.HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter); 
            if(HIj40_c30  && HIj40!=1 && PD==2 && hiBin>=60)   s.HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
            if(HIj60_c30  && !HIj60 && PD==2 && hiBin>=60)   s.HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
            if(HIj80_c30  && !HIj80 && PD==2 && hiBin>=60)   s.HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
            if(HIj100_c30 && !HIj100&& PD==2 && hiBin>=60)   s.HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
            if(HIj40_c50  && HIj40!=1 && !HIj40_c30 && PD==2 && hiBin>=100)   s.HIspec[1][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
            if(HIj60_c50  && !HIj60 && !HIj60_c30 && PD==2 && hiBin>=100)   s.HIspec[2][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
            if(HIj80_c50  && !HIj80 && !HIj80_c30 && PD==2 && hiBin>=100)   s.HIspec[3][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);
            if(HIj100_c50 && !HIj100 && !HIj100_c30 && PD==2 && hiBin>=100) s.HIspec[4][hiBin/10]->Fill(maxJtPt,trkPt[j],correction/binCenter);  
          }
        } //end trk loop
      }//end if statement  
   
      if((!isPP && PD!=2) || (isPP && PD!=1 && PD!=2))//trkTriggers
      { 
        for(int j = 0; j<nTrk; j++)
        {
          if(trkPt[j]<0.5 || trkPt[j]>=400) continue;
          if(TMath::Abs(trkEta[j])>1) continue;
          if(highPurity[j]!=1) continue;
          if(trkNHit[j]<11 || trkPtError[j]/trkPt[j]>0.1 || (int)trkOriginalAlgo[j]<4 || (int)trkOriginalAlgo[j]>7 || trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>((doChi2Shift)?chi2corr->getChi2Scale(hiBin,trkPt[j]):1)*0.15) continue; //track trigger cuts
          if(isPP){
            bool isCompatibleWithVertex = false;
            for(int v = 0; v<nVtx; v++){
              if(TMath::Abs(zVtx[v])>15) continue;
              if(TMath::Abs(trkDxyOverDxyError[j*nVtx+v])<3 && TMath::Abs(trkDzOverDzError[j*nVtx+v])<3){
                isCompatibleWithVertex = true;
                break;
              }
            }          
            if(!isCompatibleWithVertex) continue;
          }else if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
          
          float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
          if(!(trkPt[j]<caloMatchStart || (Et>caloMatchValue*trkPt[j]))) continue; //Calo Matching
          //if((maxJtPt>jetTrackCutThreshhold && trkPt[j]>maxJtPt+trkBufferSize) || (maxJtPt<=jetTrackCutThreshhold && trkPt[j]>jetTrackCutThreshhold+trkBufferSize)) continue;//upper boundary on track pt
  
          float rmin=999;
          for(int jt=0; jt<nref; jt++)
          {
            if((chargedSum[jt]/rawpt[jt]<0.01 || TMath::Abs(jteta[jt])>2)) continue;
            //if(!isPP &&  (ecalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.05 || hcalSum[jt]/(ecalSum[jt]+hcalSum[jt])<0.1|| TMath::Abs(jteta[jt])>2)) continue;
            if(jtpt[jt]<50) continue;
            float R = TMath::Power(jteta[jt]-trkEta[j],2) + TMath::Power(TMath::ACos(TMath::Cos(jtphi[jt]-trkPhi[j])),2);
            if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
          }
  
          float binCenter;
          if(isPP) binCenter = s.spec_trk[0]->GetYaxis()->GetBinCenter(s.spec[0]->GetYaxis()->FindBin(trkPt[j]));
          else     binCenter = s.HIspec_trk[0][0]->GetYaxis()->GetBinCenter(s.HIspec[0][0]->GetYaxis()->FindBin(trkPt[j]));
         

          float correction = trkCorr_trk->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j],hiBin,rmin);
          correction = correction*evtSelCorrection;
     
          //dividing by pt at bin center instead of track by track pt (just a convention)
          if(!dotrkcorr) correction=1;
          if(isPP){
            if(MinBias==1 && PD==0) s.spec_trk[0]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(t18 && PD==3)     s.spec_trk[1]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(t24 && PD==3)     s.spec_trk[2]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(t34 && PD==3)     s.spec_trk[3]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(t45 && PD==3)     s.spec_trk[4]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(t53 && PD==3)     s.spec_trk[5]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          }else{
            if(MinBias==1 && PD==0) s.HIspec_trk[0][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(HIt12 && PD==1)     s.HIspec_trk[1][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(HIt18 && PD==1)     s.HIspec_trk[2][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(HIt24 && PD==1)     s.HIspec_trk[3][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            if(HIt34 && PD==1)     s.HIspec_trk[4][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter); 
            //if(HIt12_c30 && !HIt12 && PD==1 && hiBin>=60) s.HIspec_trk[1][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
            //if(HIt18_c30 && !HIt18 && PD==1 && hiBin>=60) s.HIspec_trk[2][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
            if(HIt24_c30 && !HIt24 && PD==1 && hiBin>=60) s.HIspec_trk[3][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
            if(HIt34_c30 && !HIt34 && PD==1 && hiBin>=60) s.HIspec_trk[4][hiBin/10]->Fill(maxTrackPt,trkPt[j],correction/binCenter);
          } 
        }//end trk loop
      }//end if statement
    }//end event loop
    inputFile->Close();
  }//end file loop

  //for pp
  TFile * outF;
  if(!isPP) outF = TFile::Open(Form("PbPb_output_%d.root",jobNum),"recreate");
  else     outF = TFile::Open(Form("pp_output_%d.root",jobNum),"recreate");
  outF->cd();
  if(isPP){
    for(int i = 0; i<s.nTriggers; i++)
    {
      s.spec[i]->Write();
      s.evtCount[i]->Write();
      s.evtCount_JetVars[i]->Write();
    }
    for(int i = 0; i<s.nTriggers_trk; i++)
    {
      s.spec_trk[i]->Write();
      s.evtCount_trk[i]->Write();
    }
    s.nVtxMB->Write();
    s.nVtxMB_trk->Write();
    trkSkim->Write();
  }else{
    hiHFDist->Write();

    for(int i = 0; i<s.HInTriggers; i++)
    {
      for(int j = 0; j<20; j++){
        s.HIspec[i][j]->Write();
        s.HIevtCount[i][j]->Write();
        s.HIevtCount_JetVars[i][j]->Write();
      }
    }
    for(int i = 0; i<s.HInTriggers_trk; i++)
    {
      for(int j = 0; j<20; j++){
        s.HIspec_trk[i][j]->Write();
        s.HIevtCount_trk[i][j]->Write();
      }
    }
    for(int j = 0; j<20; j++){
      s.HInVtxMB[j]->Write();
      s.HInVtxMB_trk[j]->Write();
    }
    trkSkim->Write();
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
