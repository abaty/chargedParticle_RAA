#ifndef HYPERON
#define HYPERON

#include "TH1D.h"
#include "TFile.h"
#include <iostream>
void applyHyperonCorrection(bool isPP, TH1D* h){

  if(isPP){

  }else{
    std::cout << "Opening PbPb Hyperon Correction file..." << std::endl;
    TFile * hyperFile = TFile::Open("HyperonFractions.root","read");
    TH1D * hyperonCorr = (TH1D*)hyperFile->Get("netSyst");

    //zeroing out the lower and upper overlfow bins
    hyperonCorr->SetBinContent(0,0);
    hyperonCorr->SetBinError(0,0);
    hyperonCorr->SetBinContent(hyperonCorr->GetSize()-1,0);
    hyperonCorr->SetBinError(hyperonCorr->GetSize()-1,0);

    //applying correction
    for(int i = 1; i<h->GetSize()-1; i++)
    {
      float scale = hyperonCorr->GetBinContent(hyperonCorr->FindBin(h->GetBinCenter(i)));
      scale = 1+scale;
      h->SetBinContent(i,h->GetBinContent(i)*scale); 
      h->SetBinError(i,h->GetBinError(i)*scale); 
    }
    hyperFile->Close();
  }
  std::cout << "Done applying hyperon correction!" << std::endl;
  return;
}
#endif
