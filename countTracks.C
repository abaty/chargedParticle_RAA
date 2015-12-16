#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "getTrkCorr_simple.h"
#include "Settings.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void countTracks(std::vector<std::string> inputFiles, int jobNum, bool isTest = false)
{
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
 
  Settings s; 
  TH2D * spec[s.nTriggers];
  TH1D * evtCount[s.nTriggers];
  TH1D * nVtxMB;

  for(int i = 0; i<s.nTriggers; i++)
  {
    spec[i] = new TH2D(Form("spectrum_trigger%d",i),"",s.njetBins,0,s.maxJetBin,s.ntrkBins,s.xtrkbins);
    evtCount[i] = new TH1D(Form("evtCount%d",i),";max jet p_{T};N",s.njetBins,0,s.maxJetBin);
    evtCount[i]->SetMarkerColor(i);
  }
  nVtxMB = new TH1D("nVtxMB","nVtx;N Events",12,0,12);

  int nTrk;
  int nVtx;
  float trkPt[75000];
  float trkPtError[75000];
  float trkEta[75000];
  //float trkPhi[75000];
  bool highPurity[75000];
  float trkMVA[75000];
  float pfHcal[75000];
  float trkDxy1[100000];
  float trkDxyError1[100000];
  float trkDz1[100000];
  float trkDzError1[100000];
  float pfEcal[75000];
  unsigned char trkNHit[75000];

  int pVtx;
  int pBeamScrape;
  int NoiseFilter;

  int nref;
  float jtpt[200];
  float jteta[200];

  int MB[20]={0};
  int j40=0;
  int j60=0;
  int j80=0;

  TrkCorr* trkCorr = new TrkCorr();
  TChain * trkCh;
  TChain * jetCh;
  TChain * evtCh;
  TChain * hltCh;

  //for documenting which PD a file comes out of to avoid overlaps between PDs
  //0 is MB, 1 is jet40/60, 2 is jet80
  int PDindx[1000];
  for(unsigned int i = 0; i<inputFiles.size(); i++)
  {
    if(inputFiles.at(i).find("MinimumBias") != std::string::npos) PDindx[i]=0;
    if(inputFiles.at(i).find("HighPtLowerJets") != std::string::npos) PDindx[i]=1;
    if(inputFiles.at(i).find("HighPtJet80") != std::string::npos) PDindx[i]=2;
  }

  trkCh = new TChain("ppTrack/trackTree");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  trkCh->Add(inputFiles.at(i).c_str());
  trkCh->SetBranchAddress("nTrk",&nTrk);
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  //trkCh->SetBranchAddress("trkPhi",&trkPhi);
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

  jetCh = new TChain("ak4CaloJetAnalyzer/t");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  jetCh->Add(inputFiles.at(i).c_str());
  jetCh->SetBranchAddress("nref",&nref);
  jetCh->SetBranchAddress("jtpt",&jtpt);
  jetCh->SetBranchAddress("jteta",&jteta);  
  trkCh->AddFriend(jetCh);

  evtCh = new TChain("skimanalysis/HltTree");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  evtCh->Add(inputFiles.at(i).c_str());
  evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
  evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
  evtCh->SetBranchAddress("pHBHENoiseFilterResultProducer",&NoiseFilter);
  trkCh->AddFriend(evtCh);
  
  hltCh = new TChain("hltanalysis/HltTree");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  hltCh->Add(inputFiles.at(i).c_str());
  for(int i = 0; i<20; i++) hltCh->SetBranchAddress(Form("HLT_L1MinimumBiasHF1OR_part%d_v1",i),&(MB[i]));
  hltCh->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&j40);
  hltCh->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&j60);
  hltCh->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&j80);
  trkCh->AddFriend(hltCh);

//***********************************************************************
  std::cout << "starting event loop" << std::endl;
  std::cout << trkCh->GetEntries() << std::endl;
  for(int i = 0; i<trkCh->GetEntries(); i++)
  {
    if(i%1000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<" "<<std::endl;
    trkCh->GetEntry(i);
    if(!NoiseFilter || !pVtx || !pBeamScrape) continue;

    bool MinBias = 0;
    for(int j = 0; j<20; j++) MinBias = MinBias || (bool)MB[j];
    if(!MinBias && !j40 && !j60 && !j80) continue;
    //**************************************************
    //for trigger combination with jet triggers
    float maxJtPt = 0;
    for(int j=0; j<nref; j++)
    {
      if(jtpt[j]>maxJtPt && TMath::Abs(jteta[j])<2) maxJtPt = jtpt[j];
    }

    int PD = PDindx[trkCh->GetTreeNumber()];
    if(maxJtPt==0 && PD!=0) continue;//remove jet events where no jets are in barrel  
    if(MinBias && PD==0)
    {
      evtCount[0]->Fill(maxJtPt); 
      nVtxMB->Fill(nVtx);
    }
    if(j40 && PD==1) evtCount[1]->Fill(maxJtPt);  
    if(j60 && PD==1) evtCount[2]->Fill(maxJtPt);  
    if(j80 && PD==2) evtCount[3]->Fill(maxJtPt);  

    for(int j = 0; j<nTrk; j++)
    {
      if(TMath::Abs(trkEta[j])>1) continue;
      if(trkPt[j]<0.5 || trkPt[j]>=300) continue;
      if(highPurity[j]!=1) continue;
      if((trkMVA[j]<0.5 && trkMVA[j]!=-99) || (int)trkNHit[j]<8 || trkPtError[j]/trkPt[j]>0.3 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
      //if((trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>15 && (trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>pfHcal[j]+pfEcal[j]) continue;} //Calo Matching 

      float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j]);
      if(MinBias && PD==0) spec[0]->Fill(maxJtPt,trkPt[j],correction/trkPt[j]); 
      if(j40 && PD==1)     spec[1]->Fill(maxJtPt,trkPt[j],correction/trkPt[j]); 
      if(j60 && PD==1)     spec[2]->Fill(maxJtPt,trkPt[j],correction/trkPt[j]); 
      if(j80 && PD==2)     spec[3]->Fill(maxJtPt,trkPt[j],correction/trkPt[j]); 
    }
  }

  //for pp
  TFile * outF = TFile::Open(Form("RAA_pp_output_%d.root",jobNum),"recreate");
  outF->cd();
  for(int i = 0; i<s.nTriggers; i++)
  {
    spec[i]->Write();
    evtCount[i]->Write();
  }
  nVtxMB->Write();
  outF->Close(); 
}



//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
int main(int argc, const char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Usage: countTracks <fileList>  <job>" << std::endl;
    return 1;
  }  

  std::string fList = argv[1];
  int job = std::atoi(argv[2]);
  int totalJobs = std::atoi(argv[3]);
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
   
  countTracks(listOfFiles,job);
  return 0; 
}
