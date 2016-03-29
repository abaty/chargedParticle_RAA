
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <string>
#include "TCut.h"
#include "TMath.h"
#include <iostream>

void pileupEfficiency(){
  TH1D::SetDefaultSumw2();
  TH1D * eff[3][2];
  TH1D * gen[3][2];

  for(int i = 0; i<6; i++){
    eff[i%3][i/3] = new TH1D(Form("eff_%d_PU%d",i/3,i%3+1),";p_{T};efficiency",100,0,50);
    gen[i%3][i/3] = new TH1D(Form("gen_%d_PU%d",i/3,i%3+1),";p_{T};efficiency",100,0,50);
  }

  TFile * file;
  TTree * trkCh;
  int nFiles[6] = {2,10,10,2,10,10};
  for(int i = 0; i<6; i++){
  for(int f = 1; f<nFiles[i]; f++){
    std::cout << i << " " << f << std::endl;
    if(i==0) file = TFile::Open("root://cmsxrootd.fnal.gov//store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet50_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read"); 
    if(i==1) file = TFile::Open(Form("root://eoscms.cern.ch//store/group/phys_heavyions/abaty/Pythia8_Dijet50_pp502_TuneCUETP8M1_Dijet50_pp_Pileup2_RECODEBUG/Pythia8_Dijet50_pp502_TuneCUETP8M1/crab_20160329_145557/160329_125617/0000/HiForestAOD_%d.root",f+1),"read");
    if(i==2) file = TFile::Open(Form("root://eoscms.cern.ch//store/group/phys_heavyions/abaty/Pythia8_Dijet50_pp502_TuneCUETP8M1_Dijet50_pp_Pileup3_RECODEBUG/Pythia8_Dijet50_pp502_TuneCUETP8M1/crab_20160329_151437/160329_131456/0000/HiForestAOD_%d.root",f+1),"read");
    if(i==3) file = TFile::Open("root://cmsxrootd.fnal.gov//store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet80_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root","read"); 
    if(i==4) file = TFile::Open(Form("root://eoscms.cern.ch//store/group/phys_heavyions/abaty/Pythia8_Dijet80_pp502_TuneCUETP8M1_Dijet80_pp_Pileup2_RECODEBUG/Pythia8_Dijet80_pp502_TuneCUETP8M1/crab_20160329_152851/160329_132910/0000/HiForestAOD_%d.root",f+1),"read");
    if(i==5) file = TFile::Open(Form("root://eoscms.cern.ch//store/group/phys_heavyions/abaty/Pythia8_Dijet80_pp502_TuneCUETP8M1_Dijet80_pp_Pileup3_RECODEBUG/Pythia8_Dijet80_pp502_TuneCUETP8M1/crab_20160329_152959/160329_133010/0000/HiForestAOD_%d.root",f+1),"read");
    trkCh = (TTree*)file->Get("ppTrack/trackTree");
  
    int nTrk;
    int nEv;
    float trkPt[75000];
    float trkPtError[50000];
    float trkEta[75000];
    float trkPhi[75000];
    float trkStatus[75000]; //for trkStatus, -999 = fake, -99 = secondary, 1 & 2 are matched tracks
    bool highPurity[75000];
    float trkMVA[75000];
    float pfHcal[75000];
    float pfEcal[75000];
    float trkDxy1[75000];
    float trkDxyError1[75000];
    float trkDz1[75000];
    float trkDzError1[75000];
    float trkChi2[60000];
    unsigned char trkNHit[60000];
    unsigned char trkNlayer[60000];
    unsigned char trkAlgo[60000];
    unsigned char trkOriginalAlgo[60000];
    unsigned char trkNdof[60000];
    int nVtx;
    float zVtx[100];
    int nParticle;
    int nParticleTimesnVtx;
    float genPt[75000];
    float mtrkPt[75000];
    float genEta[75000];
    float genPhi[75000];
    bool  mtrkQual[75000];//for 5.02 samples
    float mtrkPtError[50000];
    float mtrkMVA[75000];
    float mtrkDxy1[100000];
    float mtrkDxyError1[100000];
    float mtrkDz1[100000];
    float mtrkDzError1[100000];
    bool mhighPurity[75000];

    float mtrkDzOverDzError[100000]; 
    float mtrkDxyOverDxyError[100000];

    float mtrkPfHcal[100000];
    float mtrkPfEcal[100000];
    int   mtrkNHit[60000];
    int   mtrkAlgo[60000];
    int   mtrkOriginalAlgo[60000];
    int   mtrkNlayer[60000];
    float mtrkChi2[60000];
    int   mtrkNdof[60000];
    float pNRec[75000];

    /*trkCh->SetBranchAddress("nTrk",&nTrk);
    trkCh->SetBranchAddress("nEv",&nEv);
    trkCh->SetBranchAddress("trkPt",&trkPt);
    trkCh->SetBranchAddress("trkPtError",&trkPtError);
    trkCh->SetBranchAddress("trkEta",&trkEta);
    trkCh->SetBranchAddress("trkPhi",&trkPhi);
    trkCh->SetBranchAddress("highPurity",&highPurity);
    trkCh->SetBranchAddress("trkMVA",&trkMVA);
    trkCh->SetBranchAddress("trkStatus",&trkStatus);
    trkCh->SetBranchAddress("trkDxy1",&trkDxy1);
    trkCh->SetBranchAddress("trkDxyError1",&trkDxyError1);
    trkCh->SetBranchAddress("trkDz1",&trkDz1);
    trkCh->SetBranchAddress("trkDzError1",&trkDzError1);
    trkCh->SetBranchAddress("pfHcal",&pfHcal);
    trkCh->SetBranchAddress("pfEcal",&pfEcal);
    */

    trkCh->SetBranchAddress("nParticle",&nParticle);
    trkCh->SetBranchAddress("nParticleTimesnVtx",&nParticleTimesnVtx);
    trkCh->SetBranchAddress("pPt",&genPt);
    trkCh->SetBranchAddress("pEta",&genEta);
    trkCh->SetBranchAddress("pPhi",&genPhi);
    trkCh->SetBranchAddress("pNRec",&pNRec);
    trkCh->SetBranchAddress("mtrkPt",&mtrkPt);
    trkCh->SetBranchAddress("mtrkPtError",&mtrkPtError); 
    trkCh->SetBranchAddress("mhighPurity",&mtrkQual);  //for 5.02 samples
    trkCh->SetBranchAddress("mtrkMVA",&mtrkMVA);  //for 5.02 samples
    trkCh->SetBranchAddress("mtrkDxy1",&mtrkDxy1);
    trkCh->SetBranchAddress("mtrkDxyError1",&mtrkDxyError1);
    trkCh->SetBranchAddress("mtrkDz1",&mtrkDz1);
    trkCh->SetBranchAddress("mtrkDzError1",&mtrkDzError1);
    trkCh->SetBranchAddress("mtrkDzOverDzError",&mtrkDzOverDzError);
    trkCh->SetBranchAddress("mtrkDxyOverDxyError",&mtrkDxyOverDxyError);
    trkCh->SetBranchAddress("mtrkPfHcal",&mtrkPfHcal);
    trkCh->SetBranchAddress("mtrkPfEcal",&mtrkPfEcal);
    trkCh->SetBranchAddress("nVtx",&nVtx);
    trkCh->SetBranchAddress("zVtx",&zVtx);
    trkCh->SetBranchAddress("trkNHit",&trkNHit);
    trkCh->SetBranchAddress("trkChi2",&trkChi2);
    trkCh->SetBranchAddress("trkNlayer",&trkNlayer);
    trkCh->SetBranchAddress("trkAlgo",&trkAlgo);
    trkCh->SetBranchAddress("trkOriginalAlgo",&trkOriginalAlgo);
    trkCh->SetBranchAddress("trkNdof",&trkNdof);
    trkCh->SetBranchAddress("mtrkNHit",&mtrkNHit);
    trkCh->SetBranchAddress("mtrkChi2",&mtrkChi2);
    trkCh->SetBranchAddress("mtrkNlayer",&mtrkNlayer);
    trkCh->SetBranchAddress("mtrkAlgo",&mtrkAlgo);
    trkCh->SetBranchAddress("mtrkOriginalAlgo",&mtrkOriginalAlgo);
    trkCh->SetBranchAddress("mtrkNdof",&mtrkNdof);

    for(int evt = 0; evt<TMath::Min(40000,(int)trkCh->GetEntries()); evt++){
    trkCh->GetEntry(evt);
    for(int j = 0; j<nParticle; j++){
      if(TMath::Abs(genEta[j])>1) continue;
      if(genPt[j]<0.5)            continue;
      gen[i%3][i/3]->Fill(genPt[j]);
       
 
      if((float)mtrkNHit[j]<11)          continue;
      if(mtrkChi2[j]/(float)mtrkNlayer[j]/(float)mtrkNdof[j]>0.15) continue;
      if(mtrkPtError[j]/mtrkPt[j]>0.1) continue;
      if(!(mtrkQual[j])) continue;

      float Et = (mtrkPfHcal[j]+mtrkPfEcal[j])/TMath::CosH(genEta[j]);
      if(!(mtrkPt[j]<20 || (Et>0.5*mtrkPt[j]))) continue;

      if(i!=0 && i!=3){
        bool isCompatibleWithVertex = false;
        for(int v = 0; v<nVtx; v++){
          if(TMath::Abs(zVtx[v])>15) continue;
          if(TMath::Abs(mtrkDxyOverDxyError[j*nVtx+v])<3 && TMath::Abs(mtrkDzOverDzError[j*nVtx+v])<3){
            isCompatibleWithVertex = true;
            break;
          }
        }
        if(!isCompatibleWithVertex) continue;
      }
      else{
        if(TMath::Abs(mtrkDxy1[j]/mtrkDxyError1[j])>3 || TMath::Abs(mtrkDz1[j]/mtrkDzError1[j])>3 ) continue;
      }
      eff[i%3][i/3]->Fill(genPt[j]);
    }  
    }
    file->Close();
  }
  }

  TFile * out = TFile::Open("pileup_Efficiencies.root","recreate");
  for(int i = 0; i<6; i++){
    eff[i%3][i/3]->Divide(gen[i%3][i/3]);
    eff[i%3][i/3]->Print("All");
    gen[i%3][i/3]->Print("All");
    eff[i%3][i/3]->Write();
  }
  out->Close();

  return;
}
