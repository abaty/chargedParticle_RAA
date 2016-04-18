#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  const char * histName = "PbPbTrackSpectrum_trk_0_5";
  TFile * fnew = TFile::Open("Spectra_April18_NoPbPbPu.root","read");
  TFile * fold = TFile::Open("Spectra_April18_NominalResult.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  num->Divide(den);
  num->GetYaxis()->SetTitle("(No PU 0-5%)/(With PU 0-5%)");
  num->GetYaxis()->SetRangeUser(0.7,1.3);
  num->Draw();
  c1->SetLogx();
  c1->SaveAs("plots/comparisonPlots/PbPbPileupComparison_0-5_trk.png");
  c1->SaveAs("plots/comparisonPlots/PbPbPileupComparison_0-5_trk.pdf");
}
