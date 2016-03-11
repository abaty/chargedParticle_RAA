#ifndef GETTRKCORR
#define GETTRKCORR

#include "TMath.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TrkSettings.h"
#include <iostream>
#include <vector>

class TrkCorr{
  public:
    double getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin=99, float jtpt=0, int correction=0);
    TrkCorr(std::string inputDirectory = "trkCorrections/");
    ~TrkCorr();    

  private:
    int nFiles;
    int nSteps;

    std::vector<std::vector<TH1D*> > eff;
    std::vector<std::vector<TH1D*> > fake;
    std::vector<std::vector<TH2D*> > eff2;
    std::vector<std::vector<TH2D*> > fake2;
    std::vector<TH2D*> secondary;
    std::vector<TH1D*> multiple;

    TrkSettings * s;
};

TrkCorr::TrkCorr(std::string inputDirectory)
{
  std::cout << "Initializing tracking correction files..." << std::endl;
  nFiles = 0;
  nSteps = 0;

  s = new TrkSettings(Form("%sTrkCorrInputFile.txt",inputDirectory.c_str()));
  for(int i = 0; i<s->nPtBinCoarse; i++)
  {
    nFiles += s->nCentPUBinCoarse[i];
  }
  nSteps = s->nStep;
  std::cout << "Correction tables reading " << nFiles << " Files, with " << nSteps << " iterations per file." << std::endl;

  TFile * f;
  for(int i = 0; i<nFiles; i++)
  {
    std::string isPP = "";
    if(s->nPb==0) isPP = "pp_";
    f = TFile::Open(Form("%s%scorrHists_job%d.root",inputDirectory.c_str(),isPP.c_str(),i),"read");
    std::vector<TH1D*> tempTH1EffVec;
    std::vector<TH2D*> tempTH2EffVec;
    std::vector<TH1D*> tempTH1FakeVec;
    std::vector<TH2D*> tempTH2FakeVec;
    for(int j = 0; j<nSteps; j++)
    {
      if(s->stepOrder.at(j)!=1 && s->stepOrder.at(j)!=7)
      {
	tempTH1EffVec.push_back((TH1D*)f->Get(Form("finalEff_type%d",j)));
        tempTH1EffVec.back()->SetDirectory(0);
        tempTH1FakeVec.push_back((TH1D*)f->Get(Form("finalFake_type%d",j)));
        tempTH1FakeVec.back()->SetDirectory(0);
      }
      else 
      {
	tempTH2EffVec.push_back((TH2D*)f->Get(Form("finalEff_type%d",j)));
        tempTH2EffVec.back()->SetDirectory(0);
        tempTH2FakeVec.push_back((TH2D*)f->Get(Form("finalFake_type%d",j)));
        tempTH2FakeVec.back()->SetDirectory(0);
      }
    }
    eff.push_back(tempTH1EffVec);
    fake.push_back(tempTH1FakeVec);
    eff2.push_back(tempTH2EffVec);
    fake2.push_back(tempTH2FakeVec);

    secondary.push_back((TH2D*)f->Get("SecondaryRate"));
    secondary.back()->SetDirectory(0);
    multiple.push_back((TH1D*)f->Get("MultipleRecoRate"));
    multiple.back()->SetDirectory(0);


    f->Close();
  }
  std::cout << "Initialization complete." << std::endl; 
}

//correction=0 is total, 1 is eff, 2 is fake, 3 is second, 4 is mult
double TrkCorr::getTrkCorr(float pt, float eta, float phi, float hiBin, float rmin, float jtpt, int correction)
{
  if(pt<0.5 || pt>=400){  std::cout << "\nPt of " << pt << " less than 500 MeV or > 300 GeV, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  if(eta<-2.4 || eta>2.4){  std::cout << "\nEta outside of |2.4|, please place a cut to prevent this. Returning a correction of 1" << std::endl; return 1;}
  if(hiBin<0 || hiBin>199){  std::cout << "\nhiBin not within 0 to 200, please place a cut to prevent this.  Returning a correction of 1" << std::endl; return 1;}
  
  //calculating what file to take corrections out of 
  int coarseBin = 0;
  float cent = hiBin;
  if(s->nPb==2) cent = cent/2.0;
  for(int i = 0; i<s->nPtBinCoarse; i++)
  {
    if(pt >= s->ptBinCoarse[i+1]) coarseBin+=s->nCentPUBinCoarse[i];
    else
    {
      for(int j = 0; j<s->nCentPUBinCoarse[i]; j++)
      {
        if(cent >= s->centPUBinCoarse[i][j+1]) coarseBin++;
      }
      break;
    }
  }
  //end bin calculation
 
  float netEff = 1;
  float netFake = 1;
  float netSec = 0;
  float netMult = 0; 

  int th1indx = 0;
  int th2indx = 0; 
  for(int j = 0; j<nSteps; j++)
  {
    if(s->stepOrder.at(j)==0)
    {
      netEff *= eff[coarseBin][th1indx]->GetBinContent(eff[coarseBin][th1indx]->FindBin(pt));
      netFake *= fake[coarseBin][th1indx]->GetBinContent(fake[coarseBi