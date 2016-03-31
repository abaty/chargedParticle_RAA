#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include "TLine.h"
#include "TLegend.h"
#include "TAttLine.h"

void makeCaloMatchingPlotData(bool isPP = false){
  int nEvents = 100000;

  TChain * trk; 
  TChain * h; 
  if(!isPP){
  trk = new TChain("anaTrack/trackTree");
  for(int i = 100; i<110; i++){
    trk->AddFile(Form("/mnt/hadoop/cms/store/user/velicanu/HIHardProbes/HIHardProbes-HIRun2015-PromptReco-v1-FOREST-2-v22/160126_203257/0000/HiForestAOD_%d.root",i));
  }
  h = new TChain("hiEvtAnalyzer/HiTree");
  for(int i = 100; i<110; i++){
    h->AddFile(Form("/mnt/hadoop/cms/store/user/velicanu/HIHardProbes/HIHardProbes-HIRun2015-PromptReco-v1-FOREST-2-v22/160126_203257/0000/HiForestAOD_%d.root",i));
  }
  trk->AddFriend(h);
  }else{
  trk = new TChain("ppTrack/trackTree");
  for(int i = 100; i<110; i++){
    trk->AddFile(Form("/mnt/hadoop/cms/store/user/abaty/transferTargetDirectories/2015pp_HighPtJet80/HiForestAOD_%d.root",i));
  }
  }

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  if(!isPP){ 
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
  c1->SaveAs("pdf/PbPb_data_EtVsPt.pdf");
 
  TH2D * correlation2 = new TH2D("corr2",";centrality;E_{T}/p_{T}",50,0,100,20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt:(hiBin/2)>>corr2","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2) && highPurity && trkPt>20 && TMath::Abs(trkEta)<1","colz",nEvents);
  c1->SetLogz();
  correlation2->Draw("colz");
  TLine * cut3 = new TLine(0,0.5,100,0.5);
  cut3->SetLineColor(kRed);
  cut3->SetLineWidth(5);
  cut3->Draw("same");
  c1->SaveAs("pdf/PbPb_data_EtOverPtvsHibin.pdf");
  centralProj = (TH1D*)correlation2->ProjectionY("centralProf_Data",1,1);
  centralProj->Draw();
  c1->SaveAs("pdf/PbPb_data_EtOverPtvsHibin_0_5.pdf");
  }
  else{
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
  c1->SaveAs("pdf/pp_data_EtVsPt.pdf");
  
  TH1D * correlation2 = new TH1D("corr2",";E_{T}/p_{T}",20,0,2);
  trk->Draw("(pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt>>corr2","((pfEcal+pfHcal)/TMath::CosH(trkEta)/trkPt<2) && highPurity && trkPt>20 && TMath::Abs(trkEta)<1","colz",nEvents);
  c1->SetLogz();
  correlation2->Draw();
  c1->SaveAs("pdf/pp_Data_EtOverPtvsHibin.pdf");
  }
}
