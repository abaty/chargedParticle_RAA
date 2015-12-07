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

void countTracks(std::vector<std::string> inputFiles, int jobNum)
{
///RpPb plots
  float xbins[34] = { 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 3.2 , 4.0 , 4.8 , 5.6 , 6.4 , 7.2 , 9.6 , 12.0, 14.4,19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8,73.6,86.4,103.6,120.8};

  //interpolation points
  float yvals[33] = {4.1202/0.5946*0.999, 2.8116/0.6211*1.003, 1.9778/0.6552*1.002, 1.4330/0.6984*0.994 , 1.0405/0.7219*1.002, 0.7719/0.7515*1.002, 0.5816/0.7809*1.002, 0.3893/0.825*1.001, 0.23381/0.866*1.000, 0.14395/0.901*1.004, 0.09103/0.925*1.005, 0.05937/0.965*0.997, 0.03906/0.984*0.998, 0.014787/1.023*1.062, 0.003806/1.052*1.042, 0.001181/1.056*1.033, 0.0004290/1.048*1.024, 0.0001787/1.054*1.018, 0.00008152/1.031*1.015, 0.00002216/1.023* 1.109, 0.000004653/1.036*1.061, 0.000001402/1.054*1.040, 0.0000003180/1.072*1.111, 0.00000006850/1.142*1.065, 0.00000001971/1.189*1.044, 0.000000006013/1.259*1.051, 0.000000001910/1.308*1.033, 0.0000000007181/1.342*1.024, 0.0000000002083/1.382*1.078,0.00000000005311/1.407*1.05,0.00000000001639/1.363*1.034386,0.000000000005354/1.381*1.047316,0.000000000001709/1.316*1.032760 };
  float yvals_pPb[33] = {4.1202*0.999, 2.8116*1.003, 1.9778*1.002, 1.4330*0.994 , 1.0405*1.002, 0.7719*1.002, 0.5816*1.002, 0.3893*1.001, 0.23381*1.000, 0.14395*1.004, 0.09103*1.005, 0.05937*0.997, 0.03906*0.998, 0.014787*1.062, 0.003806*1.042, 0.001181*1.033, 0.0004290*1.024, 0.0001787*1.018, 0.00008152*1.015, 0.00002216* 1.109, 0.000004653*1.061, 0.000001402*1.040, 0.0000003180*1.111, 0.00000006850*1.065, 0.00000001971*1.044, 0.000000006013*1.051, 0.000000001910*1.033, 0.0000000007181*1.024, 0.0000000002083*1.078,0.00000000005311*1.05,0.00000000001639*1.034386,0.000000000005354*1.047316,0.000000000001709*1.032760 };
  float yvals_RpPb[33] = {0.5946, 0.6211, 0.6552, 0.6984,0.7219,0.7515, 0.7809, 0.825, 0.866, 0.901, 0.925, 0.965, 0.984, 1.023, 1.052,1.056, 1.048, 1.054, 1.031, 1.023, 1.036, 1.054, 1.072, 1.142, 1.189, 1.259, 1.308, 1.342,1.382,1.407,1.363,1.381,1.316 };
  float yerrRpPb[33] = {0.0027/0.5946,0.0028/0.6211,0.0030/0.6552,0.0032/0.6984,0.0033/0.7219,0.0034/0.7515,0.0036/0.7809,0.003/0.825,0.003/0.866,0.003/0.901,0.003/0.925,0.003/0.965,0.003/0.984,0.002/1.023,0.002/1.052,0.002/1.056,0.002/1.048,0.003/1.054,0.003/1.031,0.003/1.023,0.005/1.036,0.008/1.054,0.001/1.072,0.002/1.142,0.001/1.189,0.002/1.259,0.002/1.308,0.004/1.342,0.004/1.382,0.008/1.407,0.013/1.363,0.019/1.381,0.030/1.316 };
  float yerr_pPb[33] = {0.0188/4.1202,0.0128/2.8116,0.0090/1.9778,0.0065/1.4330,0.0047/1.0405,0.0035/0.7719,0.0026/0.5816,0.0013/0.3893,0.00076/0.23381,0.00047/0.14395,0.00029/0.09103,0.00019/0.05937,0.00013/0.03906,0.000026/0.014787,0.000007/0.003806,0.000002/0.001181,0.0000009/0.0004290,0.0000004/0.0001787,0.00000025/0.00008152,0.00000025/0.00002216,0.00000006/0.000004653,0.000000023/0.000001402,0.000000011/0.0000003180,0.0000000004/0.00000006850,0.00000000011/0.00000001971,0.00000000002/0.000000006013,0.000000000009/0.000000001910,0.000000000004/0.0000000007181,0.0000000000020/0.0000000002083,0.0000000000007/0.00000000005311,0.00000000000031/0.00000000001639,0.00000000000016/0.000000000005354,0.000000000000073/0.000000000001709 };
  TH1D * RpPb = new TH1D("pp_interpolation",";pt;E#frac{d^{3}N}{dPtdyd#phi}",33,xbins);
  TH1D * RpPbSpectrum = new TH1D("RpPb_interp",";pt;RpPb",33,xbins);
  TH1D * pPbSpectrum = new TH1D("pPb_data",";pt;E#frac{d^{3}N}{dPtdyd#phi}",33,xbins);
  for(int i=1; i<34; i++) RpPb->SetBinContent(i,yvals[i-1]);
  for(int i=1; i<34; i++) RpPb->SetBinError(i,TMath::Power(TMath::Abs(TMath::Power(yerrRpPb[i-1],2)-TMath::Power(yerr_pPb[i-1],2)),0.5)*yvals[i-1]);
  for(int i=1; i<34; i++) RpPbSpectrum->SetBinContent(i,yvals_RpPb[i-1]);
  for(int i=1; i<34; i++) RpPbSpectrum->SetBinError(i,yvals_RpPb[i-1]*yerrRpPb[i-1]);
  for(int i=1; i<34; i++) pPbSpectrum->SetBinContent(i,yvals_pPb[i-1]);
  for(int i=1; i<34; i++) pPbSpectrum->SetBinError(i,yvals_pPb[i-1]*yerr_pPb[i-1]);

  TH1D::SetDefaultSumw2();

  TH1D * s[1];

    s[0] = (TH1D*)RpPb->Clone("spectrum");
    s[0]->Reset();

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

  //float density = 0;
  //float vz;

  int pVtx;
  int pBeamScrape;
  int NoiseFilter;

  int  MB=0;

  TrkCorr* trkCorr = new TrkCorr();
  TChain * trkCh;
  //TChain * centCh;
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
  std::cout << trkCh->GetEntries() << std::endl;


  evtCh = new TChain("skimanalysis/HltTree");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  evtCh->Add(inputFiles.at(i).c_str());
  evtCh->SetBranchAddress("pPAprimaryVertexFilter",&pVtx);
  evtCh->SetBranchAddress("pBeamScrapingFilter",&pBeamScrape);
  evtCh->SetBranchAddress("pHBHENoiseFilterResultProducer",&NoiseFilter);
  trkCh->AddFriend(evtCh);
  
  hltCh = new TChain("hltanalysis/HltTree");
  for(unsigned int i = 0; i<inputFiles.size(); i++)  hltCh->Add(inputFiles.at(i).c_str());
  hltCh->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part0_v1",&MB);
  trkCh->AddFriend(hltCh);

//***********************************************************************
  //TH1D * nMBTrig = new TH1D("nMBTrig","nMBTrig",2,-0.5,1.5);
  std::cout << "starting event loop" << std::endl;
  std::cout << trkCh->GetEntries() << std::endl;
  for(int i = 0; i<trkCh->GetEntries(); i++)
  {
    if(i%50000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
    trkCh->GetEntry(i);
    if(!MB) continue; 
    if(!NoiseFilter || !pVtx || !pBeamScrape) continue;

    //**************************************************
    for(int j = 0; j<nTrk; j++)
    {
      if(TMath::Abs(trkEta[j])>1) continue;
      if(trkPt[j]<0.5 || trkPt[j]>=300) continue;
      if(highPurity[j]!=1) continue;
      if((trkMVA[j]<0.5 && trkMVA[j]!=-99) || (int)trkNHit[j]<8 || trkPtError[j]/trkPt[j]>0.3 || trkDz1[j]/trkDzError1[j]>3 || trkDxy1[j]/trkDxyError1[j]>3) continue;
      if((trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>15 && (trkPt[j]-2*trkPtError[j])*TMath::CosH(trkEta[j])>pfHcal[j]+pfEcal[j]) continue; //Calo Matching 

      float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j]);
      s[0]->Fill(trkPt[j],correction/trkPt[j]);  //for minbias
    }
  }
  for(int i=1; i<34; i++) s[0]->SetBinContent(i,s[0]->GetBinContent(i)/(xbins[i]-xbins[i-1]));
  for(int i=1; i<34; i++) s[0]->SetBinError(i,s[0]->GetBinError(i)/(xbins[i]-xbins[i-1]));

  //for pp
  TFile * outF = TFile::Open(Form("output_%d.root",jobNum),"recreate");
  outF->cd();
  if(jobNum==0)
  {
    RpPbSpectrum->Write();
    pPbSpectrum->Write();
    RpPb->Write();
  }
  s[0]->Write();

  outF->Close(); 
}

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
