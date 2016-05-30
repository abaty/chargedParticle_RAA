#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  const char * histName = "PbPbMaxtrkSpectrum_0_10";
  TFile * fnew = TFile::Open("Spectra_PbPbChi2RW.root","read");
  TFile * fold = TFile::Open("Spectra_May7_withchi2NoTriCombo.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetYaxis()->SetTitle("(In Trigger Combo)/(Not in Trigger Combo)");
  num->GetYaxis()->SetRangeUser(0.7,1.3);
  num->GetXaxis()->SetRangeUser(1,200);
  num->Draw();
  c1->SetLogx(0);
  c1->SaveAs("plots/comparisonPlots/PbPb_MaxTrkSpectra_chi2_0_10.png");
  c1->SaveAs("plots/comparisonPlots/PbPb_MaxTrkSpectra_chi2_0_10.pdf");
  //c1->SaveAs("plots/comparisonPlots/ppChargeFraction2.C");
}
