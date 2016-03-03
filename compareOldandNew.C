#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  const char * histName = "PbPbTrackSpectrum_10_30";
  TFile * fnew = TFile::Open("Spectra_March3_calo0p2_jteta1p1.root","read");
  TFile * fold = TFile::Open("Spectra_March3_calo0p2_trkcorr.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  num->Divide(den);
  num->GetYaxis()->SetTitle("new/old");
  num->GetYaxis()->SetRangeUser(0.8,1.2);
  num->Draw();
  c1->SetLogx();
}
