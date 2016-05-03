#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  const char * histName = "PbPbTrackSpectrum_70_90";
  TFile * fnew = TFile::Open("Spectra_May3_noJetID_no50GeV.root","read");
  TFile * fold = TFile::Open("Spectra_April18_NominalResult.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetYaxis()->SetTitle("(70-90% No 50 GeV Cut)/(With Cut)");
  num->GetYaxis()->SetRangeUser(0.7,1.3);
  num->Draw();
  c1->SetLogx();
  c1->SaveAs("plots/comparisonPlots/50GeVCutComparison_noID_70_90.png");
  c1->SaveAs("plots/comparisonPlots/50GeVCutComparison_noID_70_90.pdf");
  //c1->SaveAs("plots/comparisonPlots/ppChargeFraction2.C");
}
