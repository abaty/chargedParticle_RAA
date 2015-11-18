#include "TH1D.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include <iostream>

const int centBins = 6;
const int centBoundary[centBins+1] = {0,10,20,60,100,140,180}; //given in hibin

const int nTriggers = 5;
const float triggerRanges[nTriggers+1] = {-1000,50,70,90,110,10000};

int getCentBin(int hiBin)
{
  for(int i=0; i<centBins; i++)
  {
    if(hiBin >= centBoundary[i] && hiBin < centBoundary[i+1]) return i;
  }
  return -1;
}

void countTracks(bool isPP, ***inputfiles***)
{
  TH1D::SetDefaultSumw2();

  TH1D * s[centBins];
  TH1D * trigCounts[centBins][nTriggers];

  const int ptBins = 51;
  double ptAxis[ptBins];
  for(int x = 0; x<ptBins;x++) ptAxis[x] = TMath::Power(10,(x*(TMath::Log10(300)-TMath::Log10(0.5))/((float)(ptBins-1))) + TMath::Log10(0.5));

  for(int i = 0; i<centBins; i++)
  {
    s[i] = new TH1D(Form("spectrum_%d_%d",centBoundary[i],centBoundary[i+1]),Form("spectrum_%d_%d",centBoundary[i],centBoundary[i+1]),ptBins,ptAxis);
    for(int j=0; j<nTriggers; j++)
    {
       trigCounts[i][j] = new TH1D(Form("eventCounts_%d_%d_trigger_%d",centBoundary[i],centBoundary[i+1],j),Form("eventCounts_%d_%d_trigger_%d",centBoundary[i],centBoundary[i+1],j),nTriggers,triggerRanges);
    } 
  }

  int nTrk;
  float trkPt[75000];
  float trkEta[75000];
  float trkPhi[75000];
  bool highPurity[75000];
  float trkMVA[75000];

  float density = 0;
  int hiBin;
  float vz;

  int nref;
  float jtpt[100];
  float jtphi[100];
  float jteta[100];
  int pcoll;
  int NoiseFilter;

  int trig1=0, trig2=0, trig3=0, trig4=0, MB=0;

  TrkCorr* trkCorr = new TrkCorr();
  TChain * trkCh;
  TChain * centCh;
  TChain * evtCh;
  TChain * jet;
  TChain * hltCh;

  if(isPP) trkCh = new TChain("ppTrack/trackTree");
  else     trkCh = new TCHain("anaTrack/trackTree");
  for(int i = 0; i<s.nMC; i++)  trkCh->Add(s.MCFiles.at(i).c_str());
  trkCh->SetBranchAddress("nTrk",&nTrk);
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  trkCh->SetBranchAddress("trkMVA",&trkMVA);

  evtCh = new TChain("skimanalysis/HltTree");
  for(int i = 0; i<s.nMC; i++)  evtCh->Add(s.MCFiles.at(i).c_str());
  evtCh->SetBranchAddress("pcollisionEventSelection",&pcoll);
  evtCh->SetBranchAddress("pHBHENoiseFilterResultProducer",&NoiseFilter);
  trkCh->AddFriend(evtCh);
 
  jet = new TChain("ak4CaloJetAnalyzer/t");
  for(int i = 0; i<s.nMC; i++)  jet->Add(s.MCFiles.at(i).c_str());
  jet->SetBranchAddress("nref",&nref);
  jet->SetBranchAddress("jtpt",&jtpt);
  jet->SetBranchAddress("jteta",&jteta);
  jet->SetBranchAddress("jtphi",&jtphi);
  trkCh->AddFriend(jet);

  hltCh = new TChain("hltanalyzer/HltTree");
  for(int i = 0; i<s.nMC; i++)  hltCh->Add(s.MCFiles.at(i).c_str());
  hltCh->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part0_v1",&MB);
  hltCh->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&trig1);
  hltCh->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&trig2);
  hltCh->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&trig3);
  hltCh->SetBranchAddress("HLT_AK4CaloJet100_Eta5p1_v1",&trig4);
  trkCh->AddFriend(hltCh);
  
  if(!isPP)
  {
    centCh = new TChain("hiEvtAnalyzer/HiTree");
    for(int i = 0; i<s.nMC; i++)  centCh->Add(s.MCFiles.at(i).c_str());  
    centCh->SetBranchAddress("vz",&vz);
    centCh->SetBranchAddress("hiBin",&hiBin);
    trkCh->AddFriend(centCh); 
  }
  else hiBin=0;


//***********************************************************************
  std::cout << "starting event loop" << std::endl;
  for(int i = 0; i<trkCh->GetEntries(); i++)
  {
    if(i%50000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
    int firesTrigger[nTriggers] = {0};
    
    if(!isPP)  centCh->GetEntry(i);
    trkCh->GetEntry(i);
  
    if(hiBin>180 || !pcoll || !NoiseFilter ) continue;

    //insert trigger paths here
    if(!MB && !trig1 && !trig2 && !trig3 && !trig4) continue;
    if(MB==1) firesTrigger[0]==1;
    if(trig1==1) firesTrigger[1]==1;
    if(trig2==1) firesTrigger[2]==1;
    if(trig3==1) firesTrigger[3]==1;
    if(trig4==1) firesTrigger[4]==1;
 
    float maxJetPt = 1;
    for(int k = 0; k<nref; k++)
    {
      if(TMath::Abs(jteta[k])>2) continue;
      if(jtpt[k]>maxJetPt) maxJetPt=jtpt[k];
    }

    trkCorr->UpdateEventInfo(trkPt,trkEta,trkPhi,nTrk);

    //**************************************************
    maxTrkPt = 0;
    for(int j = 0; j<nTrk; j++)
    {
      if(trkPt[j]>1.2*maxJetPt && maxJetPt>25) continue;            //iterative good fix
      if(highPurity[j]!=1) continue;
      if(trkMVA[j]<0.5 && trkMVA[j]!=-99) continue;  //iterative good fix
      if(trkPt[j]>maxTrkPt) maxTrkPt = trkPt[j];

      if(TMath::Abs(trkEta[j])>2.4) continue;
      if(trkPt[j]<0.5 || trkPt[j]>=300) continue;

      
      float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j]);
      for(int k = 0; k<nTriggers; k++)
      {
        if(firesTrigger[k] && maxJetPt>=triggerRanges[k] && maxJetPt<triggerRanges[k+1]) s[getCentBin(hiBin)]->Fill(trkPt[j],correction);
      }
    }
    for(int k = 0; k<nTriggers; k++)
    {
      if(firesTrigger[k]==1) trigCounts[getCentBin(hiBin)][k]->Fill(maxJetPt);
    }
  } 
}
