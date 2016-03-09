#include "TH1D.h"
#include "TAttLine.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TChain.h"

void hyperonStudy(){
  const int nBins = 19;
  const int maxPt = 20;
  int nEvts = 1000000; 
 
  //TFile * f =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/Pythia8_Dijet15_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_FOREST_758_PrivMC/HiForest_PYTHIA_QCD80_TuneCUETP8M1_cfi_5020GeV_tag_PPForestJECv6_merged/0.root","read");
  //TFile * f =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read");
  //TTree * trackTree = (TTree*)f->Get("anaTrack/trackTree");
  //TTree * trackTree = (TTree*)f->Get("ppTrack/trackTree");
  TChain * trackTree = new TChain("anaTrack/trackTree");
  for(int i = 0; i<10; i++)
  {
  trackTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9_merged_forest_%d.root",i));
  }
  
  TCanvas * c2 = new TCanvas("c2","c2",600,600);
  TH1D *num = new TH1D("num","num",nBins,0,maxPt);
  TH1D *dum = new TH1D("dum","dum",nBins,0,maxPt);
  TH1D *num1 = new TH1D("num1","num1",nBins,0,maxPt);
  TH1D *dum1 = new TH1D("dum1","dum1",nBins,0,maxPt);
  TH1D *num2 = new TH1D("num2","num2",nBins,0,maxPt);
  TH1D *dum2 = new TH1D("dum2","dum2",nBins,0,maxPt);
 
  trackTree->Draw("pPt>>num","(mhighPurity==1 && mtrkPtError/mtrkPt<0.3 && TMath::Abs(mtrkDxy1/mtrkDxyError1)<3 && TMath::Abs(mtrkDz1/mtrkDzError1)<3 && mtrkPt>0.5 && TMath::Abs(pEta)<2.4)","",nEvts);
  trackTree->Draw("pPt>>num1","((TMath::Abs(pPId)>3100 && TMath::Abs(pPId)!=3122 && TMath::Abs(pPId)<3350)?10:1)*(mhighPurity==1 && mtrkPtError/mtrkPt<0.3 && TMath::Abs(mtrkDxy1/mtrkDxyError1)<3 && TMath::Abs(mtrkDz1/mtrkDzError1)<3 && mtrkPt>0.5 && TMath::Abs(pEta)<2.4)","",nEvts);
  trackTree->Draw("pPt>>num2","((TMath::Abs(pPId)>3300 && TMath::Abs(pPId)<3350)?10:1)*(mhighPurity==1 && mtrkPtError/mtrkPt<0.3 && TMath::Abs(mtrkDxy1/mtrkDxyError1)<3 && TMath::Abs(mtrkDz1/mtrkDzError1)<3 && mtrkPt>0.5 && TMath::Abs(pEta)<2.4)","",nEvts);
  trackTree->Draw("pPt>>dum","(pPt>0.5 && TMath::Abs(pEta)<2.4)","",nEvts);
  trackTree->Draw("pPt>>dum1","((TMath::Abs(pPId)>3100 && TMath::Abs(pPId)!=3122 && TMath::Abs(pPId)<3350)?10:1)*(pPt>0.5 && TMath::Abs(pEta)<2.4)","",nEvts);
  trackTree->Draw("pPt>>dum2","((TMath::Abs(pPId)>3300 && TMath::Abs(pPId)<3350)?10:1)*(pPt>0.5 && TMath::Abs(pEta)<2.4)","",nEvts);
  
  num->Divide(dum);
  num1->Divide(dum1);
  num2->Divide(dum2);
  
  num1->SetLineColor(kBlue+1);
  num2->SetLineColor(kRed+1);
  num->GetYaxis()->SetRangeUser(0.3,0.8);
  num->GetYaxis()->SetTitle("Efficiency");
  num->GetXaxis()->SetTitle("p_{T}");

  num->Draw();
  num1->Draw("same");
  num2->Draw("same");

  TLegend * leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry((TObject*)0,"pthat 80 signal in P+H","");
  leg->AddEntry(num,"Nominal Efficiency","l");
  leg->AddEntry(num1,"#Sigma, #Xi, #Omega x10","l");
  leg->AddEntry(num2,"#Xi, #Omega x10","l");
  leg->Draw("same");
  c2->SaveAs("hyperonEff_hydjet.png");
  c2->SaveAs("hyperonEff_hydjet.pdf");
}




