#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  const char * histName = "PbPbTrackSpectrum_trk_0_5";
  TFile * fnew = TFile::Open("Spectra.root","read");
  TFile * fold = TFile::Open("Spectra_April5_TrkCorrEverywhere.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  num->Divide(den);
  num->GetYaxis()->SetTitle("(MB)/(MB+All trk Triggers)");
  num->GetYaxis()->SetRangeUser(0.7,1.3);
  num->Draw();
  c1->SetLogx();
}
