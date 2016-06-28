#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"

void compareRatios()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  
  const char * histName = "RecoRatio";
  TFile * fnew = TFile::Open("outJetAll.root","read");
  TFile * fold = TFile::Open("outJetAll_MC.root","read");
  TH2D * num = (TH2D*)fnew->Get(Form("%s",histName));
  TH2D * den = (TH2D*)fold->Get(Form("%s",histName));
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetZaxis()->SetTitle("#Delta_{HI}/#Delta_{pp}");
  num->GetZaxis()->SetRangeUser(0.9,1.1);
  num->GetYaxis()->SetRangeUser(0.7,60);
  num->Print("All");
  //num->Draw("colz");
  c1->SetLogy();
  c1->SaveAs("doubleRatio.png");
  c1->SaveAs("doubleRatio.pdf");
  //fnew->Close();
  //fold->Close();

  const int pt = 12;
  const float ptBins[pt+1] = {0.7, 1.0 , 2.0 ,  4.0  , 7.2 , 12.0,19.2, 28.8, 41.6, 60.8,86.4,120.8,165};
  TH1D * proj[3];
  for(int i = 0; i<3; i++){
    proj[i] = new TH1D(Form("proj%d",i),Form("proj%d",i),pt,ptBins);
    for(int j = 1; j<pt; j++){
      proj[i]->SetBinContent(j,num->GetBinContent(i+2,j));
      proj[i]->SetBinError(j,num->GetBinError(i+2,j));
    }
    proj[i]->Print("All");
  }

  /*const char * histName = "PbPbTrackSpectrum_0_5";
  TFile * fnew = TFile::Open("Spectra_Jun9_noChi2Cut.root","read");
  TFile * fold = TFile::Open("Spectra_Jun9_withChi2Cut.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetYaxis()->SetTitle("(No Chi2)/(With Chi2)");
  num->GetYaxis()->SetRangeUser(0.7,1.5);
  num->GetXaxis()->SetRangeUser(0.7,350);
  num->Print("All");
  num->Draw();
  c1->SetLogx();
  c1->SaveAs("plots/comparisonPlots/Chi2CutTest_PbPb_0_5.png");
  c1->SaveAs("plots/comparisonPlots/Chi2CutTest_PbPb_0_5.pdf");
  */

  //c1->SaveAs("plots/comparisonPlots/ppChargeFraction2.C");
}
