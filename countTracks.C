#include "TH1D.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include "getTrkCorr_simple.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void countTracks(std::vector<std::string> inputFiles, int jobNum, bool isTest = false)
{
  TH1D::SetDefaultSumw2();
  const int nBins = 42;
  float xbins[nBins+1] = { 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 3.2 , 4.0 , 4.8 , 5.6 , 6.4 , 7.2 , 9.6 , 12.0, 14.4,19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8,73.6,86.4,103.6,120.8,140,165,190,220,250,280,310,350,400};
  const int nTriggers = 4;
  
  TH1D * spec[nTriggers];
  TH1D * evtCount[nTriggers];

  for(int i = 0; i<nTriggers; i++)
  {
    spec[i] = new TH1D(Form("spectrum_trigger%d",i),";1/N dN/dp_T;p_T",nBins,xbins);
    evtCount[i] = new TH1D(Form("evtCount%d",i),";N;max jet p_T",200,0,1000);
  }

  int nTrk;
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

  int MB=0;
  int j40=0;
  int j60=0;
  int j80=0;

  TrkCorr* trkCorr = new TrkCorr();
  TChain * trkCh;
  TChain * jetCh;
  TChain * evtCh;
  TChain * hltCh;

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

  jetCh = new TChain("ak4CaloJetAnalzer/t");
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
  hltCh->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part0_v1",&MB);
  hltCh->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&j40);
  hltCh->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&j60);
  hltCh->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&j80);
  trkCh->AddFriend(hltCh);

//***********************************************************************
  std::cout << "starting event loop" << std::endl;
  std::cout << trkCh->GetEntries() << std::endl;
  for(int i = 0; i<100000;i++)//trkCh->GetEntries(); i++)
  {
    if(i%50000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
    trkCh->GetEntry(i);
    if(!NoiseFilter || !pVtx || !pBeamScrape) continue;
    if(!MB && !j40 && !j60 && !j80) continue;
    //**************************************************
    //for trigger combination with jet triggers
    float maxJtPt = 0;
    for(int j=0; j<nref; j++)
    {
      if(jtpt[j]>maxJtPt && TMath::Abs(jteta[j])<5.1) maxJtPt = jtpt[j];
    }
    if(MB) evtCount[0]->Fill(maxJtPt); 
    if(j40) evtCount[1]->Fill(maxJtPt);  
    if(j60) evtCount[2]->Fill(maxJtPt);   
    if(j80) evtCount[3]->Fill(maxJtPt);  

    for(int j = 0; j<nTrk; j++)
    {
      if(TMath::Abs(trkEta[j])>1) continue;
      if(trkPt[j]<0.5 || trkPt[j]>=300) continue;
      if(highPurity[j]!=1) continue;
      if((trkMVA[j]<0.5 && trkMVA[j]!=-99) || (int)trkNHit[j]<8 || trkPtError[j]/trkPt[j]>0.3 || trkDz1[j]/trkDzError1[j]>3 || trkDxy1[j]/trkDxyError1[j]>3) continue;
      //if((trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>15 && (trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>pfHcal[j]+pfEcal[j]) continue; //Calo Matching 

      float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j]);
      if(MB) spec[0]->Fill(trkPt[j],correction/trkPt[j]);  //for minbias
      if(j40) spec[1]->Fill(trkPt[j],correction/trkPt[j]);  //for minbias
      if(j60) spec[2]->Fill(trkPt[j],correction/trkPt[j]);  //for minbias
      if(j80) spec[3]->Fill(trkPt[j],correction/trkPt[j]);  //for minbias
    }
  }

  for(int i = 0; i<nTriggers; i++)
  {
     for(int j=1; j<nBins+1; j++)
     {
       spec[i]->SetBinContent(j,spec[i]->GetBinContent(j)/(xbins[j]-xbins[j-1]));
       spec[i]->SetBinError(j,spec[i]->GetBinError(j)/(xbins[j]-xbins[j-1]));
     }
  }

  //for pp
  TFile * outF = TFile::Open(Form("output_%d.root",jobNum),"recreate");
  outF->cd();
  for(int i = 0; i<nTriggers; i++) spec[i]->Write();
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
