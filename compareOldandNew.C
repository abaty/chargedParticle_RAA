#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  const char * histName = "pp_NotperMBTrigger";
  TFile * fnew = TFile::Open("Spectra_April27_chi2Reweight.root","read");
  TFile * fold = TFile::Open("Spectra_April18_NominalResult.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  num->Divide(den);
  num->GetYaxis()->SetTitle("(#xi^{2} reweighting)/(No reweighting)");
  num->GetYaxis()->SetRangeUser(0.8,1.2);
  num->Draw();
  c1->SetLogx();
  c1->SaveAs("plots/comparisonPlots/ppchi2reweighting_Comparison.png");
  c1->SaveAs("plots/comparisonPlots/ppchi2reweighting_Comparison.pdf");
}
