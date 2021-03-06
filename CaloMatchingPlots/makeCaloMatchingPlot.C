#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include "TLine.h"
#include "TLegend.h"
#include "TAttLine.h"

void makeCaloMatchingPlot(bool isPP = false){
  int nEvents = 200000;
  std::string legString = "PYTHIA dijet 120";
  std::string ppSaveString = "_pthat120";
  std::string SaveString2 = "LTT";

  TCut ptCut = "trkPt>50 && trkPt<80";
  //TCut selection = "highPurity";
  //TCut selection = "highPurity && trkChi2/trkNdof/trkNlayer<0.15 && trkNHit>10 && TMath::Abs(trkDxy1/trkDxyError1)<3 && TMath::Abs(trkDz1/trkDzError1)<3 && trkPtError/trkPt<0.1";
  TCut selection = "highPurity  && TMath::Abs(trkDxy1/trkDxyError1)<3 && TMath::Abs(trkDz1/trkDzError1)<3 && trkPtError/trkPt<0.1";

  TFile * f; 
  TTree * trk;
  TTree * hi;
  if(!isPP){
  f = TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/Pythia8_Dijet15_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_FOREST_758_PrivMC/HiForest_PYTHIA_QCD170_TuneCUETP8M1_cfi_5020GeV_tag_PPForestJECv6_merged/0.root","read");
  trk = (TTree*)f->Get("anaTrack/trackTree");
  hi = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  trk->AddFriend(hi);
  }else{
  f = TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet120_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read");
  trk = (TTree*)f->Get("ppTrack/trackTree");
  }

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  TH2D * correlation = new TH2D("corr",";p_{T};E_{T}",50,0,200,50,0,200);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta):trkPt>>corr","trkPt<400 &&((pfEcal+pfHcal)/TMath::CosH(trkEta)<400) && highPurity && TMath::Abs(trkEta)<1","colz",nEvents);
  c1->SetLogz();
  correlation->Draw("colz");
  TLine * cut1 = new TLine(20,0,20,10);
  cut1->SetLineColor(kRed);
  cut1->SetLineWidth(5);
  cut1->Draw("same");
  TLine * cut2 = new TLine(20,10,200,100);
  cut2->SetLineColor(kRed);
  cut2->SetLineWidth(5);
  cut2->Draw("same");
  /*if(isPP) c1->SaveAs("pdf/pp_pthat170_EtVsPt.pdf");
  else c1->SaveAs("pdf/pp_pthat170_EtVsPt.pdf");
  */
 
  if(!isPP){ 
  TH2D * correlation2 = new TH2D("corr2",";centrality;E_{T}/p_{T}",50,0,100,20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt:(hiBin/2)>>corr2","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2) && highPurity && trkPt>20 && TMath::Abs(trkEta)<1","colz",nEvents);
  c1->SetLogz();
  correlation2->Draw("colz");
  TLine * cut3 = new TLine(0,0.5,100,0.5);
  cut3->SetLineColor(kRed);
  cut3->SetLineWidth(5);
  cut3->Draw("same");
  c1->SaveAs("pdf/PbPb_pthat170_EtOverPtvsHibin.pdf");
  TH1D * centralProj = (TH1D*)correlation2->ProjectionY("centralProf",1,1);
  centralProj->Draw();
  c1->SaveAs("pdf/PbPb_pthat170_EtOverPtvsHibin_0_5.pdf");

  TH2D * correlation3 = new TH2D("corr3",";centrality;E_{T}/p_{T}",50,0,100,20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt:(hiBin/2)>>corr3","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2) && highPurity && trkPt>20 && trkStatus==-999 && TMath::Abs(trkEta)<1","colz",nEvents);
  c1->SetLogz();
  correlation3->Draw("colz");
  cut3->Draw("same");
  c1->SaveAs("pdf/PbPb_pthat170_EtOverPtvsHibin_Fake.pdf");
  centralProj = (TH1D*)correlation3->ProjectionY("centralProf_Fake",1,1);
  centralProj->Draw();
  c1->SaveAs("pdf/PbPb_pthat170_EtOverPtvsHibin_Fake_0_5.pdf");
  
  TH2D * correlation4 = new TH2D("corr4",";centrality;E_{T}/p_{T}",50,0,100,20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt:(hiBin/2)>>corr4","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2) && highPurity && trkPt>20 && trkStatus>0 && TMath::Abs(trkEta)<1","colz",nEvents);
  c1->SetLogz();
  correlation4->Draw("colz");
  cut3->Draw("same");
  c1->SaveAs("pdf/PbPb_pthat170_EtOverPtvsHibin_Real.pdf");
  centralProj = (TH1D*)correlation4->ProjectionY("centralProf_Real",1,1);
  centralProj->Draw();
  c1->SaveAs("pdf/PbPb_pthat170_EtOverPtvsHibin_Real_0_5.pdf");
  
  TH1D * numer = new TH1D("numer",";p_{T};Relative Efficiency",40,0,40);
  TH1D * den = new TH1D("den",";p_{T};Relative Efficiency",40,0,40);
  trk->Draw("mtrkPt>>numer","mhighPurity && TMath::Abs(pEta)<1 && mtrkPtError/mtrkPt<0.1 && mtrkNHit>10 && mtrkChi2/mtrkNlayer/mtrkNdof<0.15 && TMath::Abs(mtrkDxy1/mtrkDxyError1)<3 && TMath::Abs(mtrkDz1/mtrkDzError1)<3 && (mtrkPfEcal+mtrkPfHcal)/TMath::CosH(pEta)/mtrkPt>0.5","colz",nEvents);
  trk->Draw("mtrkPt>>den","mhighPurity && TMath::Abs(pEta)<1 && mtrkPtError/mtrkPt<0.1 && mtrkNHit>10 && mtrkChi2/mtrkNlayer/mtrkNdof<0.15 && TMath::Abs(mtrkDxy1/mtrkDxyError1)<3 && TMath::Abs(mtrkDz1/mtrkDzError1)<3","colz",nEvents);
  numer->Print("All");
  den->Print("All");
  numer->Divide(den);
  numer->Draw();
  TLegend * leg  = new TLegend(0.5,0.3,0.7,0.5);
  leg->AddEntry((TObject*)0,"PYTHIA+HYDJET","");
  leg->Draw("same");
  c1->SaveAs("pdf/PbPb_CaloMatchTurnon.pdf");
  delete leg;
  }
  else{
  TH1D * correlation2 = new TH1D("corr2",";E_{T}/p_{T}",20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt>>corr2","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2) &&  TMath::Abs(trkEta)<1" && ptCut && selection,"colz",nEvents);
  c1->SetLogz();
  c1->SetLogy();
  float scale = correlation2->GetEntries();
  correlation2->Scale(1.0/scale);
  correlation2->Draw();
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin.pdf",ppSaveString.c_str(),SaveString2.c_str()));
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin.png",ppSaveString.c_str(),SaveString2.c_str()));

  TH1D * correlation3 = new TH1D("corr3",";E_{T}/p_{T}",20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt>>corr3","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2)  && trkStatus==-999 && TMath::Abs(trkEta)<1" && ptCut && selection,"colz",nEvents);
  c1->SetLogz();
  correlation3->Scale(1.0/scale);
  correlation3->Draw();
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin_Fake.pdf",ppSaveString.c_str(),SaveString2.c_str()));
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin_Fake.png",ppSaveString.c_str(),SaveString2.c_str()));
  
  TH1D * correlation4 = new TH1D("corr4",";E_{T}/p_{T}",20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt>>corr4","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2) && trkStatus>0 && TMath::Abs(trkEta)<1" && ptCut && selection,"colz",nEvents);
  c1->SetLogz();
  correlation4->Scale(1.0/scale);
  correlation4->Draw();
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin_Real.pdf",ppSaveString.c_str(),SaveString2.c_str()));
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin_Real.png",ppSaveString.c_str(),SaveString2.c_str()));

  correlation2->GetYaxis()->SetTitle("Fraction of tracks");
  correlation2->GetYaxis()->SetRangeUser(0.00001,6);
  correlation2->Draw("");
  correlation3->SetMarkerColor(kBlue);
  correlation3->Draw("same p");
  correlation4->SetMarkerColor(kRed);
  correlation4->Draw("same p" );
  TLegend * leg2= new TLegend(0.6,0.7,0.8,0.9);
  leg2->AddEntry((TObject*)0,legString.c_str(),"");
  leg2->AddEntry(correlation2,"All Tracks","l");
  leg2->AddEntry(correlation3,"Fake","p");
  leg2->AddEntry(correlation4,"Matched","p");
  leg2->Draw("same");
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin_all.pdf",ppSaveString.c_str(),SaveString2.c_str()));
  c1->SaveAs(Form("pdf%s/%spp_EtOverPtvsHibin_all.png",ppSaveString.c_str(),SaveString2.c_str()));

  /*TH1D * numer = new TH1D("numer",";p_{T};Relative Efficiency",40,0,40);
  TH1D * den = new TH1D("den",";p_{T};Relative Efficiency",40,0,40);
  trk->Draw("mtrkPt>>numer","mhighPurity && TMath::Abs(pEta)<1 && mtrkPtError/mtrkPt<0.1 && mtrkNHit>10 && mtrkChi2/mtrkNlayer/mtrkNdof<0.15 && TMath::Abs(mtrkDxy1/mtrkDxyError1)<3 && TMath::Abs(mtrkDz1/mtrkDzError1)<3 && (mtrkPfEcal+mtrkPfHcal)/TMath::CosH(pEta)/mtrkPt>0.5","colz",nEvents);
  trk->Draw("mtrkPt>>den","mhighPurity && TMath::Abs(pEta)<1 && mtrkPtError/mtrkPt<0.1 && mtrkNHit>10 && mtrkChi2/mtrkNlayer/mtrkNdof<0.15 && TMath::Abs(mtrkDxy1/mtrkDxyError1)<3 && TMath::Abs(mtrkDz1/mtrkDzError1)<3","colz",nEvents);
  numer->Print("All");
  den->Print("All");
  numer->Divide(den);
  numer->Draw();
  TLegend * leg  = new TLegend(0.5,0.3,0.7,0.5);
  leg->AddEntry((TObject*)0,"PYTHIA","");
  leg->Draw("same");
  c1->SaveAs("pdf/pp_CaloMatchTurnon.pdf");
  delete leg;*/
  }
  return;
}
