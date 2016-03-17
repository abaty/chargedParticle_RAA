#ifndef HYPERON
#define HYPERON

#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <string>

//sign is -1 for removing correction for divinding by it
TH1D* returnHyperonCorrection(bool isPP, TH1D* h, std::string dirToCorrection=""){

  if(isPP){

  }else{
    std::cout << "Opening PbPb Hyperon Correction file..." << std::endl;
    TFile * hyperFile = TFile::Open(Form("%sHyperonFractions.root",dirToCorrection.c_str()),"read");
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
      h->SetBinContent(i,scale); 
      h->SetBinError(i,0); 
    }
    hyperFile->Close();
  }
  return h;
}
#endif
