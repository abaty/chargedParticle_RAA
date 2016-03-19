#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  const char * histName = "ppTrackSpectrum";
  TFile * fnew = TFile::Open("Spectra_March17_evtselCorrData_NoPU.root","read");
  TFile * fold = TFile::Open("Spectra_March17_evtselCorrData.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  num->Divide(den);
  num->GetYaxis()->SetTitle("pp Spectra Ratio (1 Vertex)/(Inclusive)");
  num->GetYaxis()->SetRangeUser(0.7,1.3);
  num->Draw();
  c1->SetLogx();
}
