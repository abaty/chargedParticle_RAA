#include "TCut.h"
#include "TH1D.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TChain.h"
#include "TColor.h"
#include <iostream>

double Quad(double a, double b)
{
  return TMath::Power(TMath::Power(a,2) + TMath::Power(b,2),0.5);
}

void hyperonFractions(int isPP = 0, int hLow = 0, int hHigh = 100){
  TH1::SetDefaultSumw2();
  /*int nEvts = 5000;
  int PHEvts = 10000;
  int EPOSEvts = 5000;*/
  int nEvts = 100000;
  int PHEvts = 200000;
  int EPOSEvts = 30000;
  if(isPP){
  int nEvts = 10000;
    PHEvts = 300000;
    EPOSEvts = 300000;
  }

  //hiBin actually hiHF
  float hiBinLow=0, hiBinHigh=0;
  if(hHigh>70)      hiBinLow=15;
  else if(hHigh>50) hiBinLow=123; 
  else if(hHigh>30) hiBinLow=477; 
  else if(hHigh>10) hiBinLow=1287; 
  else if(hHigh>5)  hiBinLow=3237; 
  else if(hHigh>0)  hiBinLow=4077; 
  else hiBinLow=10000; 
  
  hiBinHigh=15; 
  if(hLow<90)      hiBinHigh=123;
  if(hLow<70) hiBinHigh=477; 
  if(hLow<50) hiBinHigh=1287; 
  if(hLow<30) hiBinHigh=3237; 
  if(hLow<10)  hiBinHigh=4077; 
  if(hLow<5)  hiBinHigh=10000; 
  std::cout << hLow << " " << hHigh << " " << hiBinHigh << " " << hiBinLow << std::endl;

  TCut centCut, centCutEPOS;
  if(isPP){ centCut = "1"; centCutEPOS = "1";}
  if(!isPP){  centCutEPOS = Form("hiHF>=%f && hiHF<%f",hiBinLow,hiBinHigh); centCut = Form("hiBin>=%d*2 && hiBin<%d*2",hLow,hHigh);}
  const int nBins = 14;
  float ptBins[nBins+1] = {0.5,0.8,1.2,1.6,2,2.4,3.2,4,4.8,5.6,6.4,7.2,9.6,12,14.4};

  TFile *fP, *fEPOS, *f2, *f3, *fPH;
  TTree *PythiatrackTree, *EPOS, *PHtrackTree, *EPOShiBin, *PHtrackTreehiBin;
  TChain * trackTree, *trackTreehiBin;
  if(!isPP){ 
    fP =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read");
    PythiatrackTree = (TTree*)fP->Get("ppTrack/trackTree");

    fEPOS =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/PbPb_60kEPOS_v1/HiForest_Epos_merged_60k_v1.root","read");
    EPOS = (TTree*)fEPOS->Get("HiGenParticleAna/hi");
    EPOShiBin = (TTree*)fEPOS->Get("hiEvtAnalyzer/HiTree");
    EPOS->AddFriend(EPOShiBin);

    trackTree = new TChain("anaTrack/trackTree");
    trackTreehiBin = new TChain("hiEvtAnalyzer/HiTree");
  for(int i = 0; i<1; i++)
    {
    trackTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9_merged_forest_%d.root",i));
    trackTreehiBin->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9_merged_forest_%d.root",i));
    }
    trackTree->AddFriend(trackTreehiBin);
  }
  else{
    fP =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read");
    PythiatrackTree = (TTree*)fP->Get("ppTrack/trackTree");
    fEPOS =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/pp_300kEPOS_UseForGenStudyOnly/EPOS_pp300k.root","read");
    EPOS = (TTree*)fEPOS->Get("Events");
  } 

  TCanvas * c2 = new TCanvas("c2","c2",600,600);
  TH1D *ch = new TH1D("ch","ch",nBins,ptBins);
  TH1D *pi = new TH1D("pi","pi",nBins,ptBins);
  TH1D *k = new TH1D("k","k",nBins,ptBins);
  TH1D *sig = new TH1D("sig","sig",nBins,ptBins);
  TH1D *sigm = new TH1D("sigm","sigm",nBins,ptBins);
  TH1D *xi = new TH1D("xi","xi",nBins,ptBins);
  TH1D *omega = new TH1D("omega","omega",nBins,ptBins);
  TH1D *p = new TH1D("p","p",nBins,ptBins);
  TH1D *pbar = new TH1D("pbar","pbar",nBins,ptBins);
  TH1D *ppch = new TH1D("ppch","ppch",nBins,ptBins);
  TH1D *pppi = new TH1D("pppi","pppi",nBins,ptBins);
  TH1D *ppk = new TH1D("ppk","ppk",nBins,ptBins);
  TH1D *ppsig = new TH1D("ppsig","ppsig",nBins,ptBins);
  TH1D *ppsigm = new TH1D("ppsigm","ppsigm",nBins,ptBins);
  TH1D *ppxi = new TH1D("ppxi","ppxi",nBins,ptBins);
  TH1D *ppomega = new TH1D("ppomega","ppomega",nBins,ptBins);
  TH1D *ppp = new TH1D("ppp","ppp",nBins,ptBins);
  TH1D *pppbar = new TH1D("pppbar","pppbar",nBins,ptBins);
  TH1D *EPOSch = new TH1D("EPOSch","EPOSch",nBins,ptBins);
  TH1D *EPOSpi = new TH1D("EPOSpi","EPOSpi",nBins,ptBins);
  TH1D *EPOSk = new TH1D("EPOSk","EPOSk",nBins,ptBins);
  TH1D *EPOSsig = new TH1D("EPOSsig","EPOSsig",nBins,ptBins);
  TH1D *EPOSsigm = new TH1D("EPOSsigm","EPOSsigm",nBins,ptBins);
  TH1D *EPOSxi = new TH1D("EPOSxi","EPOSxi",nBins,ptBins);
  TH1D *EPOSomega = new TH1D("EPOSomega","EPOSomega",nBins,ptBins);
  TH1D *EPOSp = new TH1D("EPOSp","EPOSp",nBins,ptBins);
  TH1D *EPOSpbar = new TH1D("EPOSpbar","EPOSpbar",nBins,ptBins);
  TH1D *maxDiffpi = new TH1D("maxDiffpi","",nBins,ptBins);
  TH1D *maxDiffk = new TH1D("maxDiffk","",nBins,ptBins);
  TH1D *maxDiffsig = new TH1D("maxDiffsig","",nBins,ptBins);
  TH1D *maxDiffsigm = new TH1D("maxDiffsigm","",nBins,ptBins);
  TH1D *maxDiffxi = new TH1D("maxDiffxi","",nBins,ptBins);
  TH1D *maxDiffomega = new TH1D("maxDiffomega","",nBins,ptBins);
  TH1D *maxDiffp = new TH1D("maxDiffp","",nBins,ptBins);
  TH1D *maxDiffpbar = new TH1D("maxDiffpbar","",nBins,ptBins);


  if(!isPP){
    std::cout << "Hydjet" << std::endl; 
    trackTree->Draw("pPt>>ch",centCut && "pPt>0.5 && TMath::Abs(pEta)<1","",nEvts);
    trackTree->Draw("pPt>>pi",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==211 || TMath::Abs(pPId)==213)","",nEvts);
    trackTree->Draw("pPt>>k",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==321 || TMath::Abs(pPId)==323)","",nEvts);
    trackTree->Draw("pPt>>sig",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3222)","",nEvts);
    trackTree->Draw("pPt>>sigm",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3112)","",nEvts);
    trackTree->Draw("pPt>>xi",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(TMath::Abs(pPId)-3312)<=2","",nEvts);
    trackTree->Draw("pPt>>omega",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",nEvts);
    trackTree->Draw("pPt>>p",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && pPId==2212","",nEvts);
    trackTree->Draw("pPt>>pbar",centCut &&"pPt>0.5 && TMath::Abs(pEta)<1 && pPId==-2212","",nEvts);

    std::cout << "Pythia" << std::endl;
    PythiatrackTree->Draw("pPt>>ppch","pPt>0.5 && TMath::Abs(pEta)<1","",10*nEvts);
    PythiatrackTree->Draw("pPt>>pppi","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==211 || TMath::Abs(pPId)==213)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppk","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==321 || TMath::Abs(pPId)==323)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppsig","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3222)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppsigm","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3112)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppxi","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(TMath::Abs(pPId)-3312)<=2","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppomega","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppp","pPt>0.5 && TMath::Abs(pEta)<1 && (pPId)==2212","",10*nEvts);
    PythiatrackTree->Draw("pPt>>pppbar","pPt>0.5 && TMath::Abs(pEta)<1 && (pPId)==-2212","",10*nEvts);

    std::cout << "EPOS" << std::endl;
    EPOS->Draw("pt>>EPOSch",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(chg)>0","",EPOSEvts);
    EPOS->Draw("pt>>EPOSpi",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && (TMath::Abs(pdg)==211 || TMath::Abs(pdg)==213)","",EPOSEvts);
    EPOS->Draw("pt>>EPOSk",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && (TMath::Abs(pdg)==321 || TMath::Abs(pdg)==323)","",EPOSEvts);
    EPOS->Draw("pt>>EPOSsig",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && (TMath::Abs(pdg)==3222)","",EPOSEvts);
    EPOS->Draw("pt>>EPOSsigm",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && ( TMath::Abs(pdg)==3112)","",EPOSEvts);
    EPOS->Draw("pt>>EPOSxi",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(TMath::Abs(pdg)-3312)<=2","",EPOSEvts);
    EPOS->Draw("pt>>EPOSomega",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(pdg)==3334","",EPOSEvts);
    EPOS->Draw("pt>>EPOSp",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && (pdg)==2212","",EPOSEvts);
    EPOS->Draw("pt>>EPOSpbar",centCutEPOS &&"pt>0.5 && TMath::Abs(eta)<1 && (pdg)==-2212","",EPOSEvts);
  }else{
    std::cout << "Pythia" << std::endl;
    PythiatrackTree->Draw("pPt>>ppch","pPt>0.5 && TMath::Abs(pEta)<1","",10*nEvts);
    PythiatrackTree->Draw("pPt>>pppi","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==211 || TMath::Abs(pPId)==213)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppk","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==321 || TMath::Abs(pPId)==323)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppsig","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3222)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppsigm","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3112)","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppxi","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(TMath::Abs(pPId)-3312)<=2","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppomega","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",10*nEvts);
    PythiatrackTree->Draw("pPt>>ppp","pPt>0.5 && TMath::Abs(pEta)<1 && (pPId)==2212","",10*nEvts);
    PythiatrackTree->Draw("pPt>>pppbar","pPt>0.5 && TMath::Abs(pEta)<1 && (pPId)==-2212","",10*nEvts);

    std::cout << "EPOS" << std::endl;
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSch","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.charge())>0","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSpi","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && (TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId())==211 || TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId())==213)","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSk","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && (TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId())==321 || TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId())==323)","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSsig","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId())==3222","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSsigm","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId()==3112)","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSxi","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && TMath::Abs(TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId())-3312)<=2","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSomega","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId())==3334","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSp","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId()==2212","",EPOSEvts);
    EPOS->Draw("recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>>EPOSpbar","(Iteration$!=0 || Iteration$!=1) && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pt()>0.5 && TMath::Abs(recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.eta())<1 && recoGenParticles_genParticles__SIM.recoGenParticles_genParticles__SIM.obj.pdgId()==-2212","",EPOSEvts);
  }
  
  TFile *dataHyperonMeasurements;
  TH1D *dataHyperon[2];
  if(isPP){
    dataHyperonMeasurements = TFile::Open("msoverhRatio_ppData.root","read");
    dataHyperon[0] = (TH1D*)dataHyperonMeasurements->Get("XioverH_data");
    dataHyperon[1] = (TH1D*)dataHyperonMeasurements->Get("OmegaoverH_data");
    dataHyperon[0]->SetDirectory(0);
    dataHyperon[1]->SetDirectory(0);
    dataHyperonMeasurements->Close();
  }
  
  TLegend *leg;
  pppi->Divide(ppch);
  pppi->GetYaxis()->SetTitle("Ch. Particle Fraction");
  pppi->GetYaxis()->SetRangeUser(0,1);
  pppi->GetXaxis()->SetTitle("p_{T}");
  pppi->SetMarkerColor(kRed);
  pppi->SetLineColor(kRed);
  pppi->Draw();
  if(!isPP){
  pi->Divide(ch);
  pi->Draw("same");}
  EPOSpi->Divide(EPOSch);
  EPOSpi->SetMarkerColor(kBlue);
  EPOSpi->SetLineColor(kBlue);
  EPOSpi->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP) leg->AddEntry(pi,"MB Hydjet","p");
  leg->AddEntry(pppi,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSpi,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"#pi^{+/-},#rho^{+/-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_pi_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_pi_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_pi_%d_%d.C",isPP?"pp":"",hLow,hHigh));
  delete leg;

  ppk->Divide(ppch);
  ppk->SetMarkerColor(kRed);
  ppk->SetLineColor(kRed);
  ppk->GetYaxis()->SetTitle("Ch. Particle Fraction");
  ppk->GetYaxis()->SetRangeUser(0,1);
  ppk->GetXaxis()->SetTitle("p_{T}");
  ppk->Draw();
  if(!isPP){
  k->Divide(ch);
  k->Draw("same");}
  EPOSk->Divide(EPOSch);
  EPOSk->SetMarkerColor(kBlue);
  EPOSk->SetLineColor(kBlue);
  EPOSk->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP)leg->AddEntry(k,"MB Hydjet","p");
  leg->AddEntry(ppk,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSk,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"K^{+/-}, K*^{+/-} |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_k_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_k_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_k_%d_%d.C",isPP?"pp":"",hLow,hHigh));
  delete leg;

  ppsig->GetYaxis()->SetTitle("Ch. Particle Fraction");
  ppsig->GetYaxis()->SetRangeUser(0,0.3);
  ppsig->GetXaxis()->SetTitle("p_{T}");
  ppsig->Divide(ppch);
  ppsig->SetMarkerColor(kRed);
  ppsig->SetLineColor(kRed);
  ppsig->Draw();
  if(!isPP){
  sig->Divide(ch);
  sig->Draw("same");}
  EPOSsig->Divide(EPOSch);
  EPOSsig->SetMarkerColor(kBlue);
  EPOSsig->SetLineColor(kBlue);
  EPOSsig->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP)leg->AddEntry(sig,"MB Hydjet","p");
  leg->AddEntry(ppsig,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSsig,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"#Sigma^{+}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_sig_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_sig_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_sig_%d_%d.C",isPP?"pp":"",hLow,hHigh));
  delete leg;
  
  ppsigm->GetYaxis()->SetTitle("Ch. Particle Fraction");
  ppsigm->GetYaxis()->SetRangeUser(0,0.3);
  ppsigm->GetXaxis()->SetTitle("p_{T}");
  ppsigm->Divide(ppch);
  ppsigm->SetMarkerColor(kRed);
  ppsigm->SetLineColor(kRed);
  ppsigm->Draw();
  if(!isPP){
  sigm->Divide(ch);
  sigm->Draw("same");}
  EPOSsigm->Divide(EPOSch);
  EPOSsigm->SetMarkerColor(kBlue);
  EPOSsigm->SetLineColor(kBlue);
  EPOSsigm->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP)leg->AddEntry(sigm,"MB Hydjet","p");
  leg->AddEntry(ppsigm,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSsigm,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"#Sigma^{-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_sigm_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_sigm_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_sigm_%d_%d.C",isPP?"pp":"",hLow,hHigh));
  delete leg;

  ppxi->GetYaxis()->SetTitle("Ch. Particle Fraction");
  ppxi->GetYaxis()->SetRangeUser(0,0.1);
  ppxi->GetXaxis()->SetTitle("p_{T}");
  ppxi->SetMarkerColor(kRed);
  ppxi->SetLineColor(kRed);
  ppxi->Divide(ppch);
  ppxi->Draw();
  if(!isPP){
  xi->Divide(ch);
  xi->Draw("same");}
  EPOSxi->Divide(EPOSch);
  EPOSxi->SetMarkerColor(kBlue);
  EPOSxi->SetLineColor(kBlue);
  EPOSxi->Draw("same");
  if(isPP){
    dataHyperon[0]->SetMarkerColor(kBlack);
    dataHyperon[0]->SetLineColor(kBlack);
    dataHyperon[0]->SetLineWidth(2);
    dataHyperon[0]->SetMarkerSize(1);
    dataHyperon[0]->SetMarkerStyle(8);
    dataHyperon[0]->Draw("same");
  }
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP)leg->AddEntry(xi,"MB Hydjet","p");
  else     leg->AddEntry(dataHyperon[0],"pp data","p");
  leg->AddEntry(ppxi,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSxi,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"#Xi^{-}, #Xi*^{-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_xi_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_xi_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_xi_%d_%d.C",isPP?"pp":"",hLow,hHigh));
  delete leg;

  ppomega->GetYaxis()->SetTitle("Ch. Particle Fraction");
  ppomega->GetYaxis()->SetRangeUser(0,0.05);
  ppomega->GetXaxis()->SetTitle("p_{T}");
  ppomega->Divide(ppch);
  ppomega->SetMarkerColor(kRed);
  ppomega->SetLineColor(kRed);
  ppomega->Draw();
  if(!isPP){
  omega->Divide(ch);
  omega->Draw("same");}
  EPOSomega->Divide(EPOSch);
  EPOSomega->SetMarkerColor(kBlue);
  EPOSomega->SetLineColor(kBlue);
  EPOSomega->Draw("same");
  if(isPP){
    dataHyperon[1]->SetMarkerColor(kBlack);
    dataHyperon[1]->SetLineColor(kBlack);
    dataHyperon[1]->SetLineWidth(2);
    dataHyperon[1]->SetMarkerSize(1);
    dataHyperon[1]->SetMarkerStyle(8);
    dataHyperon[1]->Draw("same");
  }
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP)leg->AddEntry(omega,"MB Hydjet","p");
  else     leg->AddEntry(dataHyperon[1],"pp data","p");
  leg->AddEntry(ppomega,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSomega,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"#Omega^{-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_omega_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_omega_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_omega_%d_%d.C",isPP?"pp":"",hLow,hHigh));
  delete leg;
  
  ppp->GetYaxis()->SetTitle("Ch. Particle Fraction");
  ppp->GetYaxis()->SetRangeUser(0,0.3);
  ppp->GetXaxis()->SetTitle("p_{T}");
  ppp->Divide(ppch);
  ppp->SetMarkerColor(kRed);
  ppp->SetLineColor(kRed);
  ppp->Draw();
  if(!isPP){
  p->Divide(ch);
  p->Draw("same");}
  EPOSp->Divide(EPOSch);
  EPOSp->SetMarkerColor(kBlue);
  EPOSp->SetLineColor(kBlue);
  EPOSp->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP)leg->AddEntry(p,"MB Hydjet","p");
  leg->AddEntry(ppp,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSp,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"p, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_p_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_p_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_p_%d_%d.V",isPP?"pp":"",hLow,hHigh));
  delete leg;
  
  pppbar->GetYaxis()->SetTitle("Ch. Particle Fraction");
  pppbar->GetYaxis()->SetRangeUser(0,0.3);
  pppbar->GetXaxis()->SetTitle("p_{T}");
  pppbar->Divide(ppch);
  pppbar->SetMarkerColor(kRed);
  pppbar->SetLineColor(kRed);
  pppbar->Draw();
  if(!isPP){
  pbar->Divide(ch);
  pbar->Draw("same");}
  EPOSpbar->Divide(EPOSch);
  EPOSpbar->SetMarkerColor(kBlue);
  EPOSpbar->SetLineColor(kBlue);
  EPOSpbar->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  if(!isPP)leg->AddEntry(pbar,"MB Hydjet","p");
  leg->AddEntry(pppbar,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSpbar,Form("MB %s EPOS",isPP?"pp":"PbPb"),"p");
  leg->AddEntry((TObject*)0,"#bar{p}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs(Form("%shyperonFrac_pbar_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_pbar_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  c2->SaveAs(Form("%shyperonFrac_pbar_%d_%d.C",isPP?"pp":"",hLow,hHigh));
  delete leg;


  for(int i = 1; i<pi->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSpi->GetBinContent(i)-pi->GetBinContent(i));
    float pi_corrected_forXiOmega = pppi->GetBinContent(i);
    if(isPP){
      float binCent = ppxi->GetBinCenter(i);
      if(binCent > 10.8) binCent = 10.8;
      if(binCent < 1.3) binCent = 1.3;
      float diffXi = dataHyperon[0]->GetBinContent(dataHyperon[0]->FindBin(binCent))-ppxi->GetBinContent(i);
      float diffOmega = dataHyperon[1]->GetBinContent(dataHyperon[1]->FindBin(binCent))-ppomega->GetBinContent(i);
      
      pi_corrected_forXiOmega = pi_corrected_forXiOmega - diffXi - diffOmega; 
    }
    float diff2 = EPOSpi->GetBinContent(i)-pi_corrected_forXiOmega;
    //float diff3 = TMath::Abs(pppi->GetBinContent(i)-pi->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    maxDiffpi->SetBinContent(i,diff2);
  }
  for(int i = 1; i<k->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSk->GetBinContent(i)-k->GetBinContent(i));
    float diff2 = EPOSk->GetBinContent(i)-ppk->GetBinContent(i);
    //float diff3 = TMath::Abs(ppk->GetBinContent(i)-k->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    maxDiffk->SetBinContent(i,diff2);
  }
  for(int i = 1; i<sig->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSsig->GetBinContent(i)-sig->GetBinContent(i));
    float diff2 = EPOSsig->GetBinContent(i)-ppsig->GetBinContent(i);
    //float diff3 = TMath::Abs(ppsig->GetBinContent(i)-sig->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    maxDiffsig->SetBinContent(i,diff2);
  }
  for(int i = 1; i<sigm->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSsigm->GetBinContent(i)-sigm->GetBinContent(i));
    float diff2 = EPOSsigm->GetBinContent(i)-ppsigm->GetBinContent(i);
    //float diff3 = TMath::Abs(ppsigm->GetBinContent(i)-sigm->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    maxDiffsigm->SetBinContent(i,diff2);
  }
  for(int i = 1; i<xi->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSxi->GetBinContent(i)-xi->GetBinContent(i));
    float diff2 = EPOSxi->GetBinContent(i)-ppxi->GetBinContent(i);
    //float diff3 = TMath::Abs(ppxi->GetBinContent(i)-xi->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    if(isPP){
      float binCent = ppxi->GetBinCenter(i);
      if(binCent > 10.8) binCent = 10.8;
      if(binCent < 1.3) binCent = 1.3;
      diff2 = dataHyperon[0]->GetBinContent(dataHyperon[0]->FindBin(binCent))-ppxi->GetBinContent(i);
    }
    maxDiffxi->SetBinContent(i,diff2);
  }
  for(int i = 1; i<omega->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSomega->GetBinContent(i)-omega->GetBinContent(i));
    float diff2 = EPOSomega->GetBinContent(i)-ppomega->GetBinContent(i);
    //float diff3 = TMath::Abs(ppomega->GetBinContent(i)-omega->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    if(isPP){
      float binCent = ppomega->GetBinCenter(i);
      if(binCent > 10.8) binCent = 10.8;
      if(binCent < 1.3) binCent = 1.3;
      diff2 = dataHyperon[1]->GetBinContent(dataHyperon[1]->FindBin(binCent))-ppomega->GetBinContent(i);
    }
    maxDiffomega->SetBinContent(i,diff2);
  }
  for(int i = 1; i<p->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSp->GetBinContent(i)-p->GetBinContent(i));
    float diff2 = EPOSp->GetBinContent(i)-ppp->GetBinContent(i);
    //float diff3 = TMath::Abs(ppp->GetBinContent(i)-p->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    maxDiffp->SetBinContent(i,diff2);
  }
  for(int i = 1; i<pbar->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSp->GetBinContent(i)-p->GetBinContent(i));
    float diff2 = EPOSpbar->GetBinContent(i)-pppbar->GetBinContent(i);
    //float diff3 = TMath::Abs(ppp->GetBinContent(i)-p->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
    maxDiffpbar->SetBinContent(i,diff2);
  }
 
  if(!isPP) f2 = TFile::Open(Form("HyperonFractions_%d_%d.root",hLow,hHigh),"recreate");
  else f2 = TFile::Open("HyperonFractions_pp.root","recreate");
  pi->Write();
  k->Write();
  sig->Write();
  sigm->Write();
  xi->Write();
  omega->Write();
  p->Write();
  pbar->Write();
  pppi->Write();
  ppk->Write();
  ppsig->Write();
  ppsigm->Write();
  ppxi->Write();
  ppomega->Write();
  ppp->Write();
  pppbar->Write();
  EPOSpi->Write();
  EPOSk->Write();
  EPOSsig->Write();
  EPOSsigm->Write();
  EPOSxi->Write();
  EPOSomega->Write();
  EPOSp->Write();
  EPOSpbar->Write();
  maxDiffpi->Write();
  maxDiffk->Write();
  maxDiffsig->Write();
  maxDiffsigm->Write();
  maxDiffomega->Write();
  maxDiffxi->Write();
  maxDiffp->Write();
  maxDiffpbar->Write();
  f2->Close();
  std::cout << "finished writing step 1 of caluclation to file" << std::endl; 
 
  if(!isPP){
    fPH =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/Pythia8_Dijet15_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_FOREST_758_PrivMC/HiForest_PYTHIA_QCD80_TuneCUETP8M1_cfi_5020GeV_tag_PPForestJECv6_merged/0.root","read");
    PHtrackTree = (TTree*)fPH->Get("anaTrack/trackTree");
    PHtrackTreehiBin = (TTree*)fPH->Get("hiEvtAnalyzer/HiTree");
    PHtrackTree->AddFriend(PHtrackTreehiBin);
  }else{
    fPH =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet80_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read");
    PHtrackTree = (TTree*)fPH->Get("ppTrack/trackTree");
  }

 
  TCut ts = "!(mtrkPtError/mtrkPt>0.1 || TMath::Abs(mtrkDz1/mtrkDzError1)>3 || TMath::Abs(mtrkDxy1/mtrkDxyError1)>3 || mhighPurity==0 || (mtrkNHit<11 && mtrkPt>0.7) || mtrkChi2/mtrkNdof/mtrkNlayer>0.15)";
  //TCut ts = "!(mtrkPtError/mtrkPt>0.3 || TMath::Abs(mtrkDz1/mtrkDzError1)>3 || TMath::Abs(mtrkDxy1/mtrkDxyError1)>3 || mhighPurity==0)";
 
  TH1D *gen[9];
  TH1D *mgen[9];
  for(int i = 0; i<9; i++) gen[i] = new TH1D(Form("gen%d",i),Form("gen%d",i),nBins,ptBins);
  for(int i = 0; i<9; i++) mgen[i] = new TH1D(Form("mgen%d",i),Form("gen%d",i),nBins,ptBins);
  PHtrackTree->Draw("pPt>>gen0",centCut &&"TMath::Abs(pEta)<1","",PHEvts);
  PHtrackTree->Draw("pPt>>gen1",centCut &&"TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==211 || TMath::Abs(pPId)==213)","",PHEvts);
  PHtrackTree->Draw("pPt>>gen2",centCut &&"TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==321 || TMath::Abs(pPId)==323)","",PHEvts);
  PHtrackTree->Draw("pPt>>gen3",centCut &&"TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3222","",PHEvts);
  PHtrackTree->Draw("pPt>>gen4",centCut &&"TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3112","",PHEvts);
  PHtrackTree->Draw("pPt>>gen5",centCut &&"TMath::Abs(pEta)<1 && TMath::Abs(TMath::Abs(pPId)-3312)<=2","",PHEvts);
  PHtrackTree->Draw("pPt>>gen6",centCut &&"TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",PHEvts);
  PHtrackTree->Draw("pPt>>gen7",centCut &&"TMath::Abs(pEta)<1 && TMath::Abs(pPId)==2212","",PHEvts);
  PHtrackTree->Draw("pPt>>gen8",centCut &&"TMath::Abs(pEta)<1 && (pPId)==2212","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen0",centCut &&ts && "TMath::Abs(pEta)<1 ","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen1",centCut &&ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==211","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen2",centCut &&ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==321","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen3",centCut &&ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3222","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen4",centCut &&ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3112","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen5",centCut &&ts && "TMath::Abs(pEta)<1 && TMath::Abs(TMath::Abs(pPId)-3312)<=2","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen6",centCut &&ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen7",centCut &&ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==2212","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen8",centCut &&ts && "TMath::Abs(pEta)<1 && (pPId)==-2212","",PHEvts);

  TH1D * eff[9];
  for(int i = 0; i<9; i++){
    eff[i] = (TH1D*)mgen[i]->Clone(Form("eff%d",i));
    eff[i]->Divide(gen[i]);
  }


  TCanvas * eC = new TCanvas("eC","eC",800,600);
  TLegend * eCleg = new TLegend(0.2,0.65,0.4,0.85);
  eff[0]->GetYaxis()->SetTitle("Efficiency");
  eff[0]->GetYaxis()->SetRangeUser(0,1.4);
  eff[0]->GetXaxis()->SetTitle("p_{T}");
  eff[0]->Draw();
  int colorCount = 2;
  for(int i = 1; i<9; i++){
    if(i>2 && i<7) continue;
    if(colorCount==5) colorCount++;
    eff[i]->SetMarkerColor(colorCount);
    eff[i]->SetLineColor(colorCount);
    eff[i]->Draw("same");
    colorCount++;
  }
  eCleg->AddEntry(eff[0],"inclusive","p");
  eCleg->AddEntry(eff[1],"#pi^{+/-}","p");
  eCleg->AddEntry(eff[2],"K^{+/-}","p");
  eCleg->AddEntry(eff[7],"p","p");
  eCleg->AddEntry(eff[8],"#bar{p}","p");
  eCleg->Draw("same");
  eC->SaveAs(Form("%sefficiencyBySpecies1_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  eC->SaveAs(Form("%sefficiencyBySpecies1_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  eC->SaveAs(Form("%sefficiencyBySpecies1_%d_%d.C",isPP?"pp":"",hLow,hHigh));

  eCleg->Clear();
  eff[0]->Draw();
  colorCount = 2;
  for(int i = 3; i<7; i++){
    if(colorCount==5) colorCount++;
    eff[i]->SetMarkerColor(colorCount);
    eff[i]->SetLineColor(colorCount);
    eff[i]->Draw("same");
    colorCount++;
  }
  eCleg->AddEntry(eff[0],"inclusive","p");
  eCleg->AddEntry(eff[3],"#Sigma^{+}","p");
  eCleg->AddEntry(eff[4],"#Sigma^{-}","p");
  eCleg->AddEntry(eff[5],"#Xi^{-}","p");
  eCleg->AddEntry(eff[6],"#Omega^{-}","p");
  eCleg->Draw("same");
  eC->SaveAs(Form("%sefficiencyBySpecies2_%d_%d.png",isPP?"pp":"",hLow,hHigh));
  eC->SaveAs(Form("%sefficiencyBySpecies2_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));
  eC->SaveAs(Form("%sefficiencyBySpecies2_%d_%d.C",isPP?"pp":"",hLow,hHigh));


  TH1D * netSyst = (TH1D*)eff[1]->Clone("netSyst");
  for(int i = 1; i<gen[1]->GetSize()-1; i++)
  {
    float uncert[9] = {0};
    uncert[1] = eff[1]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffpi->GetBinContent(i)/2.0;
    uncert[2] = eff[2]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffk->GetBinContent(i)/2.0;
    uncert[3] = eff[3]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffsig->GetBinContent(i)/2.0;
    uncert[4] = eff[4]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffsigm->GetBinContent(i)/2.0;
    uncert[5] = eff[5]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffxi->GetBinContent(i)/2.0;
    uncert[6] = eff[6]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffomega->GetBinContent(i)/2.0;
    uncert[7] = eff[7]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffp->GetBinContent(i)/2.0;
    uncert[8] = eff[8]->GetBinContent(i)/eff[0]->GetBinContent(i)*maxDiffpbar->GetBinContent(i)/2.0;
    float uSum = 0;
    float quadrature = 0;
    for(int j = 1; j<9; j++){
      std::cout << j << " " << uncert[j] << std::endl; 
      uSum += uncert[j];
      quadrature = Quad(quadrature,uncert[j]);
    }
    netSyst->SetBinContent(i,1./(1+uSum)-1);
    std::cout << uSum << " " << quadrature << "\n" << std::endl;
  }

  netSyst->GetYaxis()->SetTitle("Species Fraction Correction");
  netSyst->GetXaxis()->SetTitle("p_{T}");
  c2->cd(0);
  netSyst->Draw();
  c2->SaveAs(Form("%sspeciesUncertaintyCorrection_%d_%d.png",isPP?"pp":"",hLow,hHigh));   
  c2->SaveAs(Form("%sspeciesUncertaintyCorrection_%d_%d.pdf",isPP?"pp":"",hLow,hHigh));   
  c2->SaveAs(Form("%sspeciesUncertaintyCorrection_%d_%d.C",isPP?"pp":"",hLow,hHigh));   

  if(!isPP) f3 = TFile::Open(Form("HyperonFractions_%d_%d.root",hLow,hHigh),"update");
  else      f3 = TFile::Open("HyperonFractions_pp.root","update");
  for(int i = 0; i<9; i++)
  {
    if(i==0) netSyst->Write();
    gen[i]->Write();
    mgen[i]->Write();
    eff[i]->Write();
  }
  f3->Close();
}




