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
  TH1D * pp = new TH1D("ppTrackSpectrum","ppSpectrum",s.ntrkBins,s.xtrkbins);

  //loading files
  TFile * inFile = TFile::Open("countTracks.root","read");
  for(int i = 0; i<s.nTriggers; i++)
  {
    spec[i] = (TH2D*) inFile->Get(Form("spectrum_trigger%d",i));
    evtCount[i] = (TH1D*) inFile->Get(Form("evtCount%d",i));
    spec[i]->SetDirectory(0);
    evtCount[i]->SetDirectory(0);
  }
  nVtxMB = (TH1D*) inFile->Get("nVtxMB");
  nVtxMB->SetDirectory(0);
  inFile->Close();

  //calculation of overlaps
  float scale[s.nTriggers];
  //calculate total number of verticies from MB events
  int nVtx = 0;
  for(int i = 1; i<nVtxMB->GetSize()+1;i++) nVtx = nVtx+i*nVtxMB->GetBinContent(nVtxMB->FindBin(i));
  
  for(int i = 0; i<s.nTriggers; i++)
  {
    scale[i] = 68/((float)nVtx);//using 68mb as inelastic pp xsection
    for(int j = 0; j<i; j++){
      scale[i] = scale[i]*evtCount[j]->Integral(evtCount[j]->FindBin(s.triggerBins[j+1]),evtCount[j]->FindBin(s.triggerBins[j+2]))/(double)evtCount[j+1]->Integral(evtCount[j+1]->FindBin(s.triggerBins[j+1]),evtCount[j+1]->FindBin(s.triggerBins[j+2]));
    }
    spec[i]->Scale(scale[i]);

    for(int j = evtCount[i]->FindBin(s.triggerBins[i]); j<evtCount[i]->FindBin(s.triggerBins[i+1]); j++)
    {
      for(int k = 1; k<pp->GetSize()+1; k++)
      {
        pp->SetBinContent(k,pp->GetBinContent(k)+spec[i]->GetBinContent(j,k)); 
        pp->SetBinError(k,TMath::Power(TMath::Power(pp->GetBinError(k),2)+TMath::Power(spec[i]->GetBinError(j,k),2),0.5)); 
      }
    }
  }

  //Divide by bin width and jacobian stuff
  for(int i = 1; i<pp->GetSize()+1; i++)
  {
    pp->SetBinContent(i,pp->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
    pp->SetBinError(i,pp->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
  } 
  
  TFile * outF = TFile::Open("ppSpectrum.root","recreate");
  pp->Write();
  outF->Close();
}
