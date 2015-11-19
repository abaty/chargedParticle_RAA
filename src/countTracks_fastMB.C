#include "TH1D.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranch.h"
#include "TMath.h"
#include "getTrkCorr.h"
#include <iostream>
#include <vector>

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

void countTracks(bool isPP, int nFiles)
{
///RpPb plots
  float xbins[30] = { 0.5 , 0.6 , 0.7 , 0.8 , 0.9 , 1.0 , 1.1 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 3.2 , 4.0 , 4.8 , 5.6 , 6.4 , 7.2 , 9.6 , 12.0, 14.4,19.2, 24.0, 28.8, 35.2, 41.6, 48.0, 60.8};
  float yvals[29] = {4.1202/0.5946, 2.8116/0.6211, 1.9778/0.6552, 1.4330/0.6984 , 1.0405/0.7219, 0.7719/0.7515, 0.5816/0.7809, 0.3893/0.825, 0.23381/0.866, 0.14395/0.901, 0.09103/0.925, 0.05937/0.965, 0.03906/0.984, 0.014787/1.023, 0.003806/1.052, 0.001181/1.056, 0.0004290/1.048, 0.0001787/1.054, 0.00008152/1.031, 0.00002216/1.023, 0.000004653/1.036, 0.000001402/1.054, 0.0000003180/1.072, 0.00000006850/1.142, 0.00000001971/1.189, 0.000000006013/1.259, 0.000000001910/1.308, 0.0000000007181/1.342, 0.0000000002083/1.382 };
  TH1D * RpPb = new TH1D("RpPb",";1/N dN/dPt;pt",29,xbins);
  for(int i=1; i<30; i++) RpPb->SetBinContent(i,yvals[i-1]);


  TH1D::SetDefaultSumw2();

  TH1D * s[centBins];
  TH1D * trigCounts[centBins][nTriggers];

//  const int ptBins = 51;
//  double ptAxis[ptBins];
//  for(int x = 0; x<ptBins;x++) ptAxis[x] = TMath::Power(10,(x*(TMath::Log10(300)-TMath::Log10(0.5))/((float)(ptBins-1))) + TMath::Log10(0.5));

  for(int i = 0; i<centBins; i++)
  {
//    s[i] = new TH1D(Form("spectrum_%d_%d",centBoundary[i],centBoundary[i+1]),Form("spectrum_%d_%d:pt:N",centBoundary[i],centBoundary[i+1]),ptBins-1,ptAxis);
    s[i] = (TH1D*)RpPb->Clone("spectrum");// TH1D(Form("spectrum_%d_%d",centBoundary[i],centBoundary[i+1]),Form("spectrum_%d_%d:pt:N",centBoundary[i],centBoundary[i+1]),ptBins-1,ptAxis);
    s[i]->Reset();
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

  const char * inFile = "/net/hisrv0001/home/abaty/HiForest_262163_0_numEvent1000.root";
  //const char * inFile = "root://eoscms.cern.ch//store/group/phys_heavyions/velicanu/forest/Run2015E/ExpressPhysics/FEVT/000/262/163/mergedHiForest_v0.root";

  TrkCorr* trkCorr = new TrkCorr();
  TChain * trkCh;
  TChain * centCh;
  TChain * evtCh;
  TChain * jet;
  TChain * hltCh;

  if(isPP) trkCh = new TChain("ppTrack/trackTree");
  else     trkCh = new TChain("anaTrack/trackTree");
  //for(int i = 0; i<nFiles; i++)  trkCh->Add(s.MCFiles.at(i).c_str());
  trkCh->Add(inFile);
  trkCh->SetBranchAddress("nTrk",&nTrk);
  trkCh->SetBranchAddress("trkPt",&trkPt);
  trkCh->SetBranchAddress("trkEta",&trkEta);
  trkCh->SetBranchAddress("trkPhi",&trkPhi);
  trkCh->SetBranchAddress("highPurity",&highPurity);
  trkCh->SetBranchAddress("trkMVA",&trkMVA);

  evtCh = new TChain("skimanalysis/HltTree");
  //for(int i = 0; i<nFiles; i++)  evtCh->Add(s.MCFiles.at(i).c_str());
  evtCh->Add(inFile);
  evtCh->SetBranchAddress("PAcollisionEventSelection",&pcoll);
  evtCh->SetBranchAddress("pHBHENoiseFilterResultProducer",&NoiseFilter);
  trkCh->AddFriend(evtCh);
 
  jet = new TChain("ak4CaloJetAnalyzer/t");
  //for(int i = 0; i<nFiles; i++)  jet->Add(s.MCFiles.at(i).c_str());
  jet->Add(inFile);
  jet->SetBranchAddress("nref",&nref);
  jet->SetBranchAddress("jtpt",&jtpt);
  jet->SetBranchAddress("jteta",&jteta);
  jet->SetBranchAddress("jtphi",&jtphi);
  trkCh->AddFriend(jet);

  hltCh = new TChain("hltanalysis/HltTree");
 // for(int i = 0; i<nFiles; i++)  hltCh->Add(s.MCFiles.at(i).c_str());
  hltCh->Add(inFile);
  hltCh->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part1_v1",&MB);
  hltCh->SetBranchAddress("HLT_AK4CaloJet40_Eta5p1_v1",&trig1);
  hltCh->SetBranchAddress("HLT_AK4CaloJet60_Eta5p1_v1",&trig2);
  hltCh->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1_v1",&trig3);
  hltCh->SetBranchAddress("HLT_AK4CaloJet120_Eta5p1_v1",&trig4);
  trkCh->AddFriend(hltCh);
  
  if(!isPP)
  {
    centCh = new TChain("hiEvtAnalyzer/HiTree");
    //for(int i = 0; i<nFiles; i++)  centCh->Add(s.MCFiles.at(i).c_str());  
    centCh->Add(inFile);  
    centCh->SetBranchAddress("vz",&vz);
    centCh->SetBranchAddress("hiBin",&hiBin);
    trkCh->AddFriend(centCh); 
  }
  else hiBin=0;


//***********************************************************************
  TH1D * nMBTrig = new TH1D("nMBTrig","nMBTrig",2,-0.5,1.5);
  std::cout << "starting event loop" << std::endl;
  for(int i = 0; i<trkCh->GetEntries(); i++)
  {
    if(i%50000==0) std::cout << i<<"/"<<trkCh->GetEntries()<<std::endl;
 //   int firesTrigger[nTriggers] = {0};
    
    if(!isPP)  centCh->GetEntry(i);
    trkCh->GetEntry(i);
    if(!MB) continue;  


    if(hiBin>180 || !pcoll || !NoiseFilter ) continue;
    
    //insert trigger paths here
    /*if(!MB && !trig1 && !trig2 && !trig3 && !trig4) continue;
    if(MB==1) firesTrigger[0]=1;
    if(trig1==1) firesTrigger[1]=1;
    if(trig2==1) firesTrigger[2]=1;
    if(trig3==1) firesTrigger[3]=1;
    if(trig4==1) firesTrigger[4]=1;
/
    float maxJetPt = 1;
    for(int k = 0; k<nref; k++)
    {
      if(TMath::Abs(jteta[k])>2) continue;
      if(jtpt[k]>maxJetPt) maxJetPt=jtpt[k];
    }*/
    
    nMBTrig->Fill(1);
    trkCorr->UpdateEventInfo(trkPt,trkEta,trkPhi,nTrk);

    //**************************************************
    float maxTrkPt = 0;
    for(int j = 0; j<nTrk; j++)
    {
      //if(trkPt[j]>maxJetPt && maxJetPt>15) continue;            //iterative good fix
      if(highPurity[j]!=1) continue;
      if(trkMVA[j]<0.5 && trkMVA[j]!=-99) continue;  //iterative good fix
      if(trkPt[j]>maxTrkPt) maxTrkPt = trkPt[j];

      if(TMath::Abs(trkEta[j])>1) continue;
      if(trkPt[j]<0.5 || trkPt[j]>=300) continue;

      
      float correction = trkCorr->getTrkCorr(trkPt[j],trkEta[j],trkPhi[j]);
  //    for(int k = 0; k<nTriggers; k++)
  //    {
  //      if(firesTrigger[k] && maxJetPt>=triggerRanges[k] && maxJetPt<triggerRanges[k+1]) s[getCentBin(hiBin)]->Fill(trkPt[j],correction);
  //    }
        s[getCentBin(hiBin)]->Fill(trkPt[j],correction);  //for minbias
    }
  //  for(int k = 0; k<nTriggers; k++)
  //  {
  //    if(firesTrigger[k]==1) trigCounts[getCentBin(hiBin)][k]->Fill(maxJetPt);
  //  }
      
  }

  //for pp
  TFile * outF = TFile::Open("output.root","recreate");
  outF->cd();
  RpPb->Write();
  nMBTrig->Write();
  s[0]->Scale(1/((float)nMBTrig->GetEntries()));
  float scale = 1;
  scale = RpPb->GetBinContent(1)/s[0]->GetBinContent(1);
  std::cout << RpPb->GetBinContent(1) << " " << s[0]->GetBinContent(1) << " " << scale << std::endl;
  s[0]->Scale(scale);
  s[0]->Write();
  outF->Close(); 
}
