#include "Settings.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"

void makeSpectrum()
{
  Settings s;
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  
  TH2D * spec[s.nTriggers];
  TH1D * evtCount[s.nTriggers];
  TH1D * nVtxMB;

  //loading files
  TFile * inFile = TFile::Open("Dec11_tracks.root","read");
  for(int i = 0; i<s.nTriggers; i++)
  {
    spec[i] = (TH2D*) inFile->Clone(Form("spectrum_trigger%d",i));
    evtCount[i] = (TH1D*) inFile->Clone(Form("evtCount%d",i));
    spec[i]->SetDirectory(0);
    evtCount[i]->SetDirectory(0);
  }
  nVtxMB = (TH1D*) inFile->Clone("nVtxMB");
  nVtxMB->SetDirectory(0);
  inFile->Close();

  //calculation of overlaps
  float scale[s.nTriggers];
  //calculate total number of verticies from MB events
  int nVtx = 0;
  for(int i = 1; i<nVtx->GetSize+1;i++) nVtx = nVtx+i*nVtxMB->GetBinContents(nVtxMB->FindBin(i));
  
  for(int i = 0; i<s.nTriggers; i++)
  {
    float scale[i] = 68/((float)nVtx);//using 68mb as inelastic pp xsection
    for(int j = 0; j<i; j++){
      scale[i] = scale[i]*evtCount[j]->Integral(evtCount[j]->FindBin(s.triggerBins[j+1]),evtCount[j]->FindBin(s.triggerBins[j+2]))/(double)evtCount[j+1]->Integral(evtCount[j+1]->FindBin(s.triggerBins[j+1]),evtCount[j+1]->FindBin(s.triggerBins[j+2]))
    }
    spec[i]->Scale(scale[i]);
  }
}
