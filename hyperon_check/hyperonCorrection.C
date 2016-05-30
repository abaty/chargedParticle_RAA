#ifndef HYPERON
#define HYPERON

#include "TH1D.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <string>

TH1D* returnHyperonCorrection(bool isPP, TH1D* h, int centindx, std::string dirToCorrection=""){

  TFile * hyperFile;
  TH1D * hyperonCorr;
  if(isPP){
    std::cout << "Opening pp Hyperon Correction file..." << std::endl;
    hyperFile = TFile::Open(Form("%s/HyperonFractions_pp.root",dirToCorrection.c_str()),"read");
    hyperonCorr = (TH1D*)hyperFile->Get("netSyst");
  }else{
    std::cout << "Opening PbPb Hyperon Correction file..." << std::endl;
    int hyperLow[8] = {0,5,10,30,50,70,0,0};
    int hyperHigh[8] = {5,10,30,50,90,90,10,100};
    if(centindx <0 || centindx>7){std::cout << "Bad hyperon indx" << std::endl;  centindx=7;}
    hyperFile = TFile::Open(Form("%s/HyperonFractions_%d_%d.root",dirToCorrection.c_str(),hyperLow[centindx],hyperHigh[centindx]),"read");
    hyperonCorr = (TH1D*)hyperFile->Get("netSyst");
  }
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
  return h;
}
#endif
