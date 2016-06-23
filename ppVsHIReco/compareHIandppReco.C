#include <map>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

inline long getHash(int evt,int run){
  return evt+(run-262620)*100000000000;
}

inline int getIndex(std::map<long,int> Map, int evt,int run){
  if(Map.find(getHash(evt,run))!=Map.end()) return Map.find(getHash(evt,run))->second;
  else return -1;
}

inline void clearIndex(std::map<long,int> Map, int evt,int run){
  Map.erase(getHash(evt,run));
  return;
}

void buildMap(std::map<long,int> Map){
  int nEv, nRun;

  TFile * inf = TFile::Open("/mnt/hadoop/cms/store/user/abaty/unmergedForests/2015PbPb_HIHardProbesPeripheralFiltered_ppReco_Forest/HIHardProbesPeripheral/crab_20160623_193101/160623_173111/0000/HiForestAOD_1.root","read");
  TTree * trkTree = (TTree*)inf->Get("ppTrack/trackTree"); 
  trkTree->SetBranchAddress("nEv",&nEv);
  trkTree->SetBranchAddress("nRun",&nRun);

  std::cout << "Building Centrality Map..." << std::endl;
  for(int i = 0; i<trkTree->GetEntries(); i++){
    trkTree->GetEntry(i);
    if(i%1000==0) std::cout << i <<"/" << trkTree->GetEntries()<< std::endl; 
    long hash = getHash(nEv,nRun);
    Map.insert(std::pair<long,int>(hash,i));
  }
  std::cout << "Done Building Map" << std::endl;  


  std::cout << getIndex(Map,nEv,nRun) << std::endl;

  inf->Close();
  return;
}

void compareHIandppReco(bool isTest = true){
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  std::map<long,int> Map;
  buildMap(Map);

  std::string fList;
  if(isTest) fList = "HIRecoFilesTest.txt";
  else       fList = "HIRecoFiles.txt";  
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return ;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
      line++;
    }
  } 
    
    TFile * ppFile = TFile::Open("/mnt/hadoop/cms/store/user/abaty/unmergedForests/2015PbPb_HIHardProbesPeripheralFiltered_ppReco_Forest/HIHardProbesPeripheral/crab_20160623_193101/160623_173111/0000/HiForestAOD_1.root","read");
    TTree * pptrkCh = (TTree*)ppFile->Get("ppTrack/trackTree");
    int nTrk2;
    int nVtx2;
    int nTrkTimesnVtx2;
    bool highPurity2[50000];
    float trkPt2[50000];
    float trkPtError2[50000];
    float trkEta2[50000];
    float trkPhi2[50000];
    float trkMVA2[50000];
    float trkDxy12[50000];
    float trkDxyError12[50000];
    float trkDz12[50000];
    float trkDzError12[50000];
    float trkDzOverDzError2[500000];
    float trkDxyOverDxyError2[500000];
    float pfEcal2[50000];
    float pfHcal2[50000];
    float trkChi22[50000];
    float zVtx2[20];
    unsigned char trkNHit2[50000];
    unsigned char trkNlayer2[50000];
    unsigned char trkNdof2[50000];
    unsigned char trkAlgo2[50000];
    unsigned char trkOriginalAlgo2[50000];
    pptrkCh = (TTree*)ppFile->Get("ppTrack/trackTree");
    pptrkCh->SetBranchAddress("nTrk",&nTrk2);
    pptrkCh->SetBranchAddress("nVtx",&nVtx2);
    pptrkCh->SetBranchAddress("trkPt",&trkPt2);
    pptrkCh->SetBranchAddress("trkEta",&trkEta2);
    pptrkCh->SetBranchAddress("trkPhi",&trkPhi2);
    pptrkCh->SetBranchAddress("highPurity",&highPurity2);
    pptrkCh->SetBranchAddress("trkMVA",&trkMVA2);
    pptrkCh->SetBranchAddress("trkNHit",&trkNHit2);
    pptrkCh->SetBranchAddress("trkPtError",&trkPtError2);
    pptrkCh->SetBranchAddress("pfHcal",&pfHcal2);
    pptrkCh->SetBranchAddress("pfEcal",&pfEcal2);
    pptrkCh->SetBranchAddress("trkDxy1",&trkDxy12);
    pptrkCh->SetBranchAddress("trkDxyError1",&trkDxyError12);
    pptrkCh->SetBranchAddress("trkDz1",&trkDz12);
    pptrkCh->SetBranchAddress("trkDzError1",&trkDzError12);
    pptrkCh->SetBranchAddress("trkChi2",&trkChi22);
    pptrkCh->SetBranchAddress("trkNlayer",&trkNlayer2);
    pptrkCh->SetBranchAddress("trkNdof",&trkNdof2);
    pptrkCh->SetBranchAddress("trkAlgo",&trkAlgo2);
    pptrkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo2);
    pptrkCh->SetBranchAddress("zVtx",&zVtx2);

  std::cout << "Number of HI Reco Files: " << listOfFiles.size() << std::endl;
  for(int f = 0; f<listOfFiles.size(); f++){
    std::cout << "File: " << f << "/" << listOfFiles.size() << std::endl;
    TFile * inputFile = TFile::Open(listOfFiles[f].data(),"read");
    TTree * hltCh = (TTree*)inputFile->Get("hltanalysis/HltTree");
    
    int HIj40_c30=0, HIj60_c30=0, HIj80_c30=0;
    int HIj40_c50=0, HIj60_c50=0, HIj80_c50=0;
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1",&HIj40_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1",&HIj60_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1",&HIj80_c30);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1",&HIj40_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1",&HIj60_c50);
    hltCh->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1",&HIj80_c50);

    int hiBin = -1;
    UInt_t run;
    ULong64_t evt;
    TTree * hiCh = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    hiCh->SetBranchAddress("hiBin",&hiBin);
    hiCh->SetBranchAddress("run",&run);
    hiCh->SetBranchAddress("evt",&evt);

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
    TTree * trkCh = (TTree*)inputFile->Get("anaTrack/trackTree");
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

    for(int i = 0; i<hltCh->GetEntries(); i++){
      hltCh->GetEntry(i);
      if(!(HIj80_c30 || HIj80_c50)) continue;
      hiCh->GetEntry(i); 
      int ppIndex = -1;
      if(getIndex(Map,evt,run)<0) continue; 
      else{
        ppIndex=getIndex(Map,evt,run);
        clearIndex(Map,evt,run);
      }
      trkCh->GetEntry(i);     
      pptrkCh->GetEntry(ppIndex);
      std::cout << i << " " << ppIndex << " " << hiBin << " " << nTrk << " " << nTrk2 << std::endl;
    }
 
    inputFile->Close();
  }
}
