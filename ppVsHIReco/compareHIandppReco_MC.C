#include <map>
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

bool passesTrigger(int HIj40_c30, int HIj40_c50, int HIj60_c30, int HIj60_c50, int HIj80_c30, int HIj80_c50){
  if(HIj40_c30 || HIj40_c50 || HIj60_c30 || HIj60_c50 || HIj80_c30 || HIj80_c50) return true;
  else return false;
}

long getHash(int evt,int run){
  return evt+(run-262620)*100000000000;
}

void clearIndex(std::map<long,int> Map, int evt,int run){
  Map.erase(getHash(evt,run));
  return;
}

void compareHIandppReco_MC(bool isTest = true, bool makeMap = false, bool isMC = true){
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  const int cent = 4;
  const float centBins[cent+1] = {25,30,50,70,90};
  const int pt = 12;
  const float ptBins[pt+1] = {0.7, 1.0 , 2.0 ,  4.0  , 7.2 , 12.0,19.2, 28.8, 41.6, 60.8,86.4,120.8,165};
  
  int HIj40_c30=0, HIj60_c30=0, HIj80_c30=0;
  int HIj40_c50=0, HIj60_c50=0, HIj80_c50=0;
 
  //cent map
  std::map<long,int> Map;

  if(makeMap){
    int nEv, nRun;
    TFile * inf;
    if(!isMC) inf = TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/ppReco_HIHardProbesPeripheral2015/0.root","read");
    else      inf = TFile::Open("/mnt/hadoop/cms/store/user/abaty/ppReco_diJet50PrivateSamples/ppRecoForest/Baty_ppReco_diJet50PrivateSamples_GENSIM/crab_20160627_154502/160627_134528/0000/HiForestAOD_10.root","read");
    TTree * trkTree = (TTree*)inf->Get("ppTrack/trackTree"); 
    trkTree->SetBranchAddress("nEv",&nEv);
    trkTree->SetBranchAddress("nRun",&nRun);
    TTree * hltTree = (TTree*)inf->Get("hltanalysis/HltTree");
    hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent30_100_v1",&HIj40_c30);
    hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent30_100_v1",&HIj60_c30);
    hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent30_100_v1",&HIj80_c30);
    hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet40_Eta5p1_Cent50_100_v1",&HIj40_c50);
    hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet60_Eta5p1_Cent50_100_v1",&HIj60_c50);
    hltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_Cent50_100_v1",&HIj80_c50);
    trkTree->AddFriend(hltTree);
  
    std::cout << "Building Centrality Map..." << std::endl;
    for(int i = 0; i<trkTree->GetEntries(); i++){
      if(i%1000==0) std::cout << i <<"/" << trkTree->GetEntries()<< std::endl; 
      trkTree->GetEntry(i);
      //if(!passesTrigger(HIj40_c30, HIj40_c50, HIj60_c30, HIj60_c50, HIj80_c30, HIj80_c50)) continue;
      long hash = getHash(nEv,nRun);
      Map.insert(std::pair<long,int>(hash,i));
    }
    std::cout << "Done Building Map" << std::endl;  
    inf->Close();
    
    //write map
    ofstream centMap;
    if(!isMC) centMap.open("centMap.txt");
    else      centMap.open("centMap_MC.txt");
    for (std::map<long,int>::iterator it=Map.begin(); it!=Map.end(); ++it){
      centMap << it->first << " " << it->second << '\n';
    }
    centMap.close();
  }else{
    std::cout << "Loading cent Map" << std::endl;
    std::fstream centMap;
    if(!isMC) centMap.open("centMap.txt");
    else      centMap.open("centMap_MC.txt");
    if (centMap.is_open())
    {
      int tmpIndex;
      long tmpHash;
      while (centMap >> tmpHash >> tmpIndex)
      {
        //std::cout << tmpIndex <<" "<< tmpHash << std::endl;
        Map.insert(std::pair<long,int>(tmpHash,tmpIndex));  
      }
      centMap.close();
    }
    
  /*int counter = 0;
    for (std::map<long,int>::iterator it=Map.begin(); it!=Map.end(); it++){
      counter++;
      if(counter>100) break;
      std::cout << it->first << " " << it->second << std::endl;
    } */
    std::cout << "Done Loading map, size: "<< Map.size() << std::endl;
  }
  //end map

  std::cout << "opening HI files" << std::endl;
  std::string fList;
  fList = "HIRecoFiles_MC.txt";
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
  std::cout << "opening pp files" << std::endl;
  std::string fListpp;
  fList = "ppRecoFiles_MC.txt";
  std::string bufferpp;
  std::vector<std::string> listOfFilespp;
  std::ifstream inFilepp(fList.data());

  if(!inFilepp.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return ;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFilepp >> bufferpp;
      if(inFilepp.eof()) break;
      listOfFilespp.push_back(bufferpp);
      line++;
    }
  } 
  
    TH2D * HIReco = new TH2D("HIReco","HIReco",cent,centBins,pt,ptBins);
    TH2D * ppReco = new TH2D("ppReco","ppReco",cent,centBins,pt,ptBins);
    HIReco->SetDirectory(0);
    ppReco->SetDirectory(0); 
   

  std::cout << "Number of HI Reco Files: " << listOfFiles.size() << std::endl;
  std::cout << "Number of pp Reco Files: " << listOfFilespp.size() << std::endl;
  for(int f = 0; f<listOfFiles.size(); f++){
    TFile * ppFile;
    ppFile = TFile::Open(listOfFilespp[f].data(),"read");
    TTree * pptrkCh = (TTree*)ppFile->Get("ppTrack/trackTree");
    int nTrk2;
    int nVtx2;
    int nTrkTimesnVtx2;
    bool highPurity2[15000];
    float trkPt2[15000];
    float trkPtError2[15000];
    float trkEta2[15000];
    //float trkPhi2[15000];
    //float trkMVA2[15000];
    float trkDxy12[15000];
    float trkDxyError12[15000];
    float trkDz12[15000];
    float trkDzError12[15000];
    float trkDzOverDzError2[150000];
    float trkDxyOverDxyError2[150000];
    float pfEcal2[15000];
    float pfHcal2[15000];
    //float trkChi22[15000];
    float zVtx2[20];
    //unsigned char trkNHit2[15000];
    //unsigned char trkNlayer2[15000];
    //unsigned char trkNdof2[15000];
    //unsigned char trkAlgo2[15000];
    //unsigned char trkOriginalAlgo2[15000];
    pptrkCh->SetBranchAddress("nTrk",&nTrk2);
    pptrkCh->SetBranchAddress("nVtx",&nVtx2);
    pptrkCh->SetBranchAddress("trkPt",&trkPt2);
    pptrkCh->SetBranchAddress("trkEta",&trkEta2);
    //pptrkCh->SetBranchAddress("trkPhi",&trkPhi2);
    pptrkCh->SetBranchAddress("highPurity",&highPurity2);
    //pptrkCh->SetBranchAddress("trkMVA",&trkMVA2);
    //pptrkCh->SetBranchAddress("trkNHit",&trkNHit2);
    pptrkCh->SetBranchAddress("trkPtError",&trkPtError2);
    pptrkCh->SetBranchAddress("pfHcal",&pfHcal2);
    pptrkCh->SetBranchAddress("pfEcal",&pfEcal2);
    pptrkCh->SetBranchAddress("trkDxy1",&trkDxy12);
    pptrkCh->SetBranchAddress("trkDxyError1",&trkDxyError12);
    pptrkCh->SetBranchAddress("trkDz1",&trkDz12);
    pptrkCh->SetBranchAddress("trkDzError1",&trkDzError12);
    //pptrkCh->SetBranchAddress("trkChi2",&trkChi22);
    //pptrkCh->SetBranchAddress("trkNlayer",&trkNlayer2);
    //pptrkCh->SetBranchAddress("trkNdof",&trkNdof2);
    //pptrkCh->SetBranchAddress("trkAlgo",&trkAlgo2);
    //pptrkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo2);
    pptrkCh->SetBranchAddress("zVtx",&zVtx2);
    std::cout << "here" << std::endl;
    

    std::cout << "File: " << f << "/" << listOfFiles.size() << std::endl;
    TFile * inputFile = TFile::Open(listOfFiles[f].data(),"read");
    TTree * hltCh = (TTree*)inputFile->Get("hltanalysis/HltTree");
    
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
    bool highPurity[15000];
    float trkPt[15000];
    float trkPtError[15000];
    float trkEta[15000];
    float trkPhi[15000];
    //float trkMVA[15000];
    float trkDxy1[15000];
    float trkDxyError1[15000];
    float trkDz1[15000];
    float trkDzError1[15000];
    float trkDzOverDzError[115000];
    float trkDxyOverDxyError[150000];
    float pfEcal[15000];
    float pfHcal[15000];
    float trkChi2[15000];
    float zVtx[20];
    unsigned char trkNHit[15000];
    unsigned char trkNlayer[15000];
    unsigned char trkNdof[15000];
    //unsigned char trkAlgo[15000];
    //unsigned char trkOriginalAlgo[15000];
    TTree * trkCh = (TTree*)inputFile->Get("anaTrack/trackTree");
    trkCh->SetBranchAddress("nTrk",&nTrk);
    trkCh->SetBranchAddress("nVtx",&nVtx);
    trkCh->SetBranchAddress("trkPt",&trkPt);
    trkCh->SetBranchAddress("trkEta",&trkEta);
    trkCh->SetBranchAddress("trkPhi",&trkPhi);
    trkCh->SetBranchAddress("highPurity",&highPurity);
    //trkCh->SetBranchAddress("trkMVA",&trkMVA);
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
    //trkCh->SetBranchAddress("trkAlgo",&trkAlgo);
    //trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo);
    trkCh->SetBranchAddress("zVtx",&zVtx);

    for(int i = 0; i<hltCh->GetEntries(); i++){
      hltCh->GetEntry(i);
      //if(!passesTrigger(HIj40_c30, HIj40_c50, HIj60_c30, HIj60_c50, HIj80_c30, HIj80_c50)) continue;

      //syncing for centrality
      hiCh->GetEntry(i);
      //if(Map.count((long)getHash(evt,run))<1) continue;
      //std::cout << i << " " << run << " " << evt << " " <<getHash((int)evt,(int)run) <<" " << Map.at((long)getHash(evt,run)) << std::endl;
      //int ppIndex = Map.at((long)getHash(evt,run));

      //track loops
      trkCh->GetEntry(i);     
      pptrkCh->GetEntry(i);
      //if(i%100==0) std::cout << i << " " << ppIndex  << " " << nTrk << " " << nTrk2<< " " << hiBin << std::endl;

      //HI
      for(int j = 0; j<nTrk; j++){
        if(!highPurity[j]) continue;
        if(TMath::Abs(trkEta[j])>1) continue;
        if(trkPtError[j]/trkPt[j]>0.1) continue;
        if(TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
        if(TMath::Abs(trkDz1[j]/trkDzError1[j])>3) continue;
        if(trkNHit[j]<11) continue;
        if(trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>0.15) continue;
        float Et = (pfEcal[j]+pfHcal[j])/TMath::CosH(trkEta[j]);
        if(!(trkPt[j]<20 || (Et>0.5*trkPt[j]))) continue;

        HIReco->Fill(hiBin/2.0,trkPt[j],1);
      }
       
      //pp
      for(int j = 0; j<nTrk2; j++){
        if(!highPurity2[j]) continue;
        if(TMath::Abs(trkEta2[j])>1) continue;
        if(trkPtError2[j]/trkPt2[j]>0.3) continue;
        if(TMath::Abs(trkDxy12[j]/trkDxyError12[j])>3) continue;
        if(TMath::Abs(trkDz12[j]/trkDzError12[j])>3) continue;
        float Et = (pfEcal2[j]+pfHcal2[j])/TMath::CosH(trkEta2[j]);
        if(!(trkPt2[j]<20 || (Et>0.5*trkPt2[j]))) continue;
   
        ppReco->Fill(hiBin/2.0,trkPt2[j],1);
      }
      clearIndex(Map,evt,run);
    }
    inputFile->Close();
    ppFile->Close();
  }

  TFile * out; 
  //if(!isMC) out = TFile::Open("outJetAll.root","recreate");
  out = TFile::Open("outJetAll_MC.root","recreate");
  //TH2D * test123 = new TH2D("test123","test123",10,0,10,10,0,10);
  //test123->Write();
  HIReco->Write();
  ppReco->Write();
  TH2D * Ratio = (TH2D*) HIReco->Clone("RecoRatio");
  Ratio->Divide(ppReco);
  Ratio->Write();
  std::cout << "Done" << std::endl;
  //out->Close();
}
