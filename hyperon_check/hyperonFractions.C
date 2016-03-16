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
#include <iostream>

double Quad(double a, double b)
{
  return TMath::Power(TMath::Power(a,2) + TMath::Power(b,2),0.5);
}

void hyperonFractions(){
  TH1::SetDefaultSumw2();
  /*int nEvts = 500;
  int PHEvts = 1000;
  int EPOSEvts = 500; */
  int nEvts = 100000;
  int PHEvts = 200000;
  int EPOSEvts = 60000; 

  const int nBins = 14;
  float ptBins[nBins+1] = {0.5,0.8,1.2,1.6,2,2.4,3.2,4,4.8,5.6,6.4,7.2,9.6,12,14.4};
 
  TFile * fP =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet30_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read");
  TTree * PythiatrackTree = (TTree*)fP->Get("ppTrack/trackTree");

  TFile * fEPOS =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/PbPb_60kEPOS_v1/HiForest_Epos_merged_60k_v1.root","read");
  TTree * EPOS = (TTree*)fEPOS->Get("HiGenParticleAna/hi");

  TChain * trackTree = new TChain("anaTrack/trackTree");
  for(int i = 0; i<1; i++)
  {
  trackTree->Add(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9/HiForest_Centrality_Unpacker_Hydjet_Quenched_MinBias_5020GeV_750_RECODEBUG_v0_TAGHiForestPbPbJECv9_merged_forest_%d.root",i));
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


  std::cout << "Hydjet" << std::endl; 
  trackTree->Draw("pPt>>ch","pPt>0.5 && TMath::Abs(pEta)<1","",nEvts);
  trackTree->Draw("pPt>>pi","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==211","",nEvts);
  trackTree->Draw("pPt>>k","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==321","",nEvts);
  trackTree->Draw("pPt>>sig","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3222)","",nEvts);
  trackTree->Draw("pPt>>sigm","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3112)","",nEvts);
  trackTree->Draw("pPt>>xi","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3312","",nEvts);
  trackTree->Draw("pPt>>omega","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",nEvts);
  trackTree->Draw("pPt>>p","pPt>0.5 && TMath::Abs(pEta)<1 && pPId==2212","",nEvts);
  trackTree->Draw("pPt>>pbar","pPt>0.5 && TMath::Abs(pEta)<1 && pPId==-2212","",nEvts);

  std::cout << "Pythia" << std::endl;
  PythiatrackTree->Draw("pPt>>ppch","pPt>0.5 && TMath::Abs(pEta)<1","",10*nEvts);
  PythiatrackTree->Draw("pPt>>pppi","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==211","",10*nEvts);
  PythiatrackTree->Draw("pPt>>ppk","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==321","",10*nEvts);
  PythiatrackTree->Draw("pPt>>ppsig","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3222)","",10*nEvts);
  PythiatrackTree->Draw("pPt>>ppsigm","pPt>0.5 && TMath::Abs(pEta)<1 && (TMath::Abs(pPId)==3112)","",10*nEvts);
  PythiatrackTree->Draw("pPt>>ppxi","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3312","",10*nEvts);
  PythiatrackTree->Draw("pPt>>ppomega","pPt>0.5 && TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",10*nEvts);
  PythiatrackTree->Draw("pPt>>ppp","pPt>0.5 && TMath::Abs(pEta)<1 && (pPId)==2212","",10*nEvts);
  PythiatrackTree->Draw("pPt>>pppbar","pPt>0.5 && TMath::Abs(pEta)<1 && (pPId)==-2212","",10*nEvts);

  std::cout << "EPOS" << std::endl;
  EPOS->Draw("pt>>EPOSch","pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(chg)>0","",EPOSEvts);
  EPOS->Draw("pt>>EPOSpi","pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(pdg)==211","",EPOSEvts);
  EPOS->Draw("pt>>EPOSk","pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(pdg)==321","",EPOSEvts);
  EPOS->Draw("pt>>EPOSsig","pt>0.5 && TMath::Abs(eta)<1 && (TMath::Abs(pdg)==3222)","",EPOSEvts);
  EPOS->Draw("pt>>EPOSsigm","pt>0.5 && TMath::Abs(eta)<1 && ( TMath::Abs(pdg)==3112)","",EPOSEvts);
  EPOS->Draw("pt>>EPOSxi","pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(pdg)==3312","",EPOSEvts);
  EPOS->Draw("pt>>EPOSomega","pt>0.5 && TMath::Abs(eta)<1 && TMath::Abs(pdg)==3334","",EPOSEvts);
  EPOS->Draw("pt>>EPOSp","pt>0.5 && TMath::Abs(eta)<1 && (pdg)==2212","",EPOSEvts);
  EPOS->Draw("pt>>EPOSpbar","pt>0.5 && TMath::Abs(eta)<1 && (pdg)==-2212","",EPOSEvts);
  
  TLegend *leg;
  pi->Divide(ch);
  pi->GetYaxis()->SetTitle("Ch. Particle Fraction");
  pi->GetYaxis()->SetRangeUser(0,1);
  pi->GetXaxis()->SetTitle("p_{T}");
  pi->Draw();
  pppi->Divide(ppch);
  pppi->SetMarkerColor(kRed);
  pppi->SetLineColor(kRed);
  pppi->Draw("same");
  EPOSpi->Divide(EPOSch);
  EPOSpi->SetMarkerColor(kBlue);
  EPOSpi->SetLineColor(kBlue);
  EPOSpi->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(pi,"MB Hydjet","p");
  leg->AddEntry(pppi,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSpi,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"#pi^{+/-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_pi.png");
  c2->SaveAs("hyperonFrac_pi.pdf");
  delete leg;

  k->Divide(ch);
  k->GetYaxis()->SetTitle("Ch. Particle Fraction");
  k->GetYaxis()->SetRangeUser(0,1);
  k->GetXaxis()->SetTitle("p_{T}");
  k->Draw();
  ppk->Divide(ppch);
  ppk->SetMarkerColor(kRed);
  ppk->SetLineColor(kRed);
  ppk->Draw("same");
  EPOSk->Divide(EPOSch);
  EPOSk->SetMarkerColor(kBlue);
  EPOSk->SetLineColor(kBlue);
  EPOSk->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(k,"MB Hydjet","p");
  leg->AddEntry(ppk,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSk,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"K^{+/-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_k.png");
  c2->SaveAs("hyperonFrac_k.pdf");
  delete leg;

  sig->Divide(ch);
  sig->GetYaxis()->SetTitle("Ch. Particle Fraction");
  sig->GetYaxis()->SetRangeUser(0,0.3);
  sig->GetXaxis()->SetTitle("p_{T}");
  sig->Draw();
  ppsig->Divide(ppch);
  ppsig->SetMarkerColor(kRed);
  ppsig->SetLineColor(kRed);
  ppsig->Draw("same");
  EPOSsig->Divide(EPOSch);
  EPOSsig->SetMarkerColor(kBlue);
  EPOSsig->SetLineColor(kBlue);
  EPOSsig->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(sig,"MB Hydjet","p");
  leg->AddEntry(ppsig,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSsig,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"#Sigma^{+}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_sig.png");
  c2->SaveAs("hyperonFrac_sig.pdf");
  delete leg;
  
  sigm->Divide(ch);
  sigm->GetYaxis()->SetTitle("Ch. Particle Fraction");
  sigm->GetYaxis()->SetRangeUser(0,0.3);
  sigm->GetXaxis()->SetTitle("p_{T}");
  sigm->Draw();
  ppsigm->Divide(ppch);
  ppsigm->SetMarkerColor(kRed);
  ppsigm->SetLineColor(kRed);
  ppsigm->Draw("same");
  EPOSsigm->Divide(EPOSch);
  EPOSsigm->SetMarkerColor(kBlue);
  EPOSsigm->SetLineColor(kBlue);
  EPOSsigm->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(sigm,"MB Hydjet","p");
  leg->AddEntry(ppsigm,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSsigm,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"#Sigma^{-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_sigm.png");
  c2->SaveAs("hyperonFrac_sigm.pdf");
  delete leg;

  xi->Divide(ch);
  xi->GetYaxis()->SetTitle("Ch. Particle Fraction");
  xi->GetYaxis()->SetRangeUser(0,0.1);
  xi->GetXaxis()->SetTitle("p_{T}");
  xi->Draw();
  ppxi->Divide(ppch);
  ppxi->SetMarkerColor(kRed);
  ppxi->SetLineColor(kRed);
  ppxi->Draw("same");
  EPOSxi->Divide(EPOSch);
  EPOSxi->SetMarkerColor(kBlue);
  EPOSxi->SetLineColor(kBlue);
  EPOSxi->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(xi,"MB Hydjet","p");
  leg->AddEntry(ppxi,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSxi,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"#Xi^{-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_xi.png");
  c2->SaveAs("hyperonFrac_xi.pdf");
  delete leg;

  omega->Divide(ch);
  omega->GetYaxis()->SetTitle("Ch. Particle Fraction");
  omega->GetYaxis()->SetRangeUser(0,0.05);
  omega->GetXaxis()->SetTitle("p_{T}");
  omega->Draw();
  ppomega->Divide(ppch);
  ppomega->SetMarkerColor(kRed);
  ppomega->SetLineColor(kRed);
  ppomega->Draw("same");
  EPOSomega->Divide(EPOSch);
  EPOSomega->SetMarkerColor(kBlue);
  EPOSomega->SetLineColor(kBlue);
  EPOSomega->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(omega,"MB Hydjet","p");
  leg->AddEntry(ppomega,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSomega,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"#Omega^{-}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_omega.png");
  c2->SaveAs("hyperonFrac_omega.pdf");
  delete leg;
  
  p->Divide(ch);
  p->GetYaxis()->SetTitle("Ch. Particle Fraction");
  p->GetYaxis()->SetRangeUser(0,0.05);
  p->GetXaxis()->SetTitle("p_{T}");
  p->Draw();
  ppp->Divide(ppch);
  ppp->SetMarkerColor(kRed);
  ppp->SetLineColor(kRed);
  ppp->Draw("same");
  EPOSp->Divide(EPOSch);
  EPOSp->SetMarkerColor(kBlue);
  EPOSp->SetLineColor(kBlue);
  EPOSp->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(p,"MB Hydjet","p");
  leg->AddEntry(ppp,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSp,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"p, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_p.png");
  c2->SaveAs("hyperonFrac_p.pdf");
  delete leg;
  
  pbar->Divide(ch);
  pbar->GetYaxis()->SetTitle("Ch. Particle Fraction");
  pbar->GetYaxis()->SetRangeUser(0,0.05);
  pbar->GetXaxis()->SetTitle("p_{T}");
  pbar->Draw();
  pppbar->Divide(ppch);
  pppbar->SetMarkerColor(kRed);
  pppbar->SetLineColor(kRed);
  pppbar->Draw("same");
  EPOSpbar->Divide(EPOSch);
  EPOSpbar->SetMarkerColor(kBlue);
  EPOSpbar->SetLineColor(kBlue);
  EPOSpbar->Draw("same");
  leg = new TLegend(0.2,0.7,0.6,0.9);
  leg->AddEntry(pbar,"MB Hydjet","p");
  leg->AddEntry(pppbar,"PYTHIA Dijet 30","p");
  leg->AddEntry(EPOSpbar,"MB PbPb EPOS","p");
  leg->AddEntry((TObject*)0,"#bar{p}, |#eta|<1","");
  leg->Draw("same");
  c2->SaveAs("hyperonFrac_pbar.png");
  c2->SaveAs("hyperonFrac_pbar.pdf");
  delete leg;

  for(int i = 1; i<pi->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSpi->GetBinContent(i)-pi->GetBinContent(i));
    float diff2 = EPOSpi->GetBinContent(i)-pppi->GetBinContent(i);
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
    maxDiffxi->SetBinContent(i,diff2);
  }
  for(int i = 1; i<omega->GetSize()-1; i++){
    //float diff1 = TMath::Abs(EPOSomega->GetBinContent(i)-omega->GetBinContent(i));
    float diff2 = EPOSomega->GetBinContent(i)-ppomega->GetBinContent(i);
    //float diff3 = TMath::Abs(ppomega->GetBinContent(i)-omega->GetBinContent(i));
    //float diff = TMath::Max(diff3,TMath::Max(diff1,diff2));
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
  
  TFile *f2 = TFile::Open("HyperonFractions.root","recreate");
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
  
  TFile * fPH =  TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/Pythia8_Dijet15_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_FOREST_758_PrivMC/HiForest_PYTHIA_QCD80_TuneCUETP8M1_cfi_5020GeV_tag_PPForestJECv6_merged/0.root","read");
  TTree * PHtrackTree = (TTree*)fPH->Get("anaTrack/trackTree");

 
  TCut ts = "!(mtrkPtError/mtrkPt>0.1 || TMath::Abs(mtrkDz1/mtrkDzError1)>3 || TMath::Abs(mtrkDxy1/mtrkDxyError1)>3 || mhighPurity==0 || (mtrkNHit<11 && mtrkPt>0.7) || mtrkChi2/mtrkNdof/mtrkNlayer>0.15)";
 
  TH1D *gen[9];
  TH1D *mgen[9];
  for(int i = 0; i<9; i++) gen[i] = new TH1D(Form("gen%d",i),Form("gen%d",i),nBins,ptBins);
  for(int i = 0; i<9; i++) mgen[i] = new TH1D(Form("mgen%d",i),Form("gen%d",i),nBins,ptBins);
  PHtrackTree->Draw("pPt>>gen0","TMath::Abs(pEta)<1","",PHEvts);
  PHtrackTree->Draw("pPt>>gen1","TMath::Abs(pEta)<1 && TMath::Abs(pPId)==211","",PHEvts);
  PHtrackTree->Draw("pPt>>gen2","TMath::Abs(pEta)<1 && TMath::Abs(pPId)==321","",PHEvts);
  PHtrackTree->Draw("pPt>>gen3","TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3222","",PHEvts);
  PHtrackTree->Draw("pPt>>gen4","TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3112","",PHEvts);
  PHtrackTree->Draw("pPt>>gen5","TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3312","",PHEvts);
  PHtrackTree->Draw("pPt>>gen6","TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",PHEvts);
  PHtrackTree->Draw("pPt>>gen7","TMath::Abs(pEta)<1 && TMath::Abs(pPId)==2212","",PHEvts);
  PHtrackTree->Draw("pPt>>gen8","TMath::Abs(pEta)<1 && (pPId)==2212","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen0",ts && "TMath::Abs(pEta)<1 ","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen1",ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==211","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen2",ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==321","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen3",ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3222","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen4",ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3112","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen5",ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3312","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen6",ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==3334","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen7",ts && "TMath::Abs(pEta)<1 && TMath::Abs(pPId)==2212","",PHEvts);
  PHtrackTree->Draw("pPt>>mgen8",ts && "TMath::Abs(pEta)<1 && (pPId)==-2212","",PHEvts);

  TH1D * eff[9];
  for(int i = 0; i<9; i++){
    eff[i] = (TH1D*)mgen[i]->Clone(Form("eff%d",i));
    eff[i]->Divide(gen[i]);
  }


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
  netSyst->Draw();
  c2->SaveAs("speciesUncertaintyCorrection.png");   
  c2->SaveAs("speciesUncertaintyCorrection.pdf");   

  TFile *f3 = TFile::Open("HyperonFractions.root","update");
  for(int i = 0; i<9; i++)
  {
    if(i==0) netSyst->Write();
    gen[i]->Write();
    mgen[i]->Write();
    eff[i]->Write();
  }
  f3->Close();
}




