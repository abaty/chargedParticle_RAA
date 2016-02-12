#include "Settings.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TAttFill.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TAttAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "THStack.h"
#include "ppPlotting.C"
#include "PbPbPlotting.C"

void makeSpectrum()
{
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  Settings s;

//output hist for jet triggered
  s.pp = new TH1D("ppTrackSpectrum",";p_{T} (GeV);E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xtrkbins);
  for(int i = 0; i<s.nTriggers; i++) s.ppByTrigger[i] = new TH1D(Form("ppTrackSpectrumByTrigger%d",i),"",s.ntrkBins,s.xtrkbins);
  for(int i = 0; i<s.nTriggers; i++) s.ppUsedByTrigger[i] = new TH1D(Form("ppUsedTrackSpectrumByTrigger%d",i),"",s.ntrkBins,s.xtrkbins);
  s.ppJets = new TH1D("ppJetSpectrum",";Leading Jet P_{T} (GeV);#sigma (mb)",s.njetBins,0,s.maxJetBin);
  for(int i = 0; i<s.nTriggers; i++) s.ppJetsByTrigger[i] = new TH1D(Form("ppJetSpectrumByTrigger%d",i),"",s.njetBins,0,s.maxJetBin);

//output hist for track triggers 
  s.pp_trk = new TH1D("ppTrackSpectrum_trk",";p_{T} (GeV);E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xtrkbins);
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppByTrigger_trk[i] = new TH1D(Form("ppTrackSpectrumByTrigger%d_trk",i),"",s.ntrkBins,s.xtrkbins);
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppUsedByTrigger_trk[i] = new TH1D(Form("ppUsedTrackSpectrumByTrigger%d_trk",i),"",s.ntrkBins,s.xtrkbins);
  s.ppMaxtrk = new TH1D("ppMaxtrkSpectrum",";Leading Track P_{T} (GeV);#sigma (mb)",s.nTrktriggerBins,0,s.maxTrktriggerBin);
  for(int i = 0; i<s.nTriggers_trk; i++) s.ppMaxtrkByTrigger[i] = new TH1D(Form("ppMaxtrkSpectrumByTrigger%d",i),"",s.nTrktriggerBins,0,s.maxTrktriggerBin);

//same for PbPb
  for(int c = 0 ; c<s.nCentBins; c++){
    s.HI[c] = new TH1D(Form("PbPbTrackSpectrum_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]),";p_{T} (GeV);E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xtrkbins);
    for(int i = 0; i<s.HInTriggers; i++) s.HIByTrigger[i][c] = new TH1D(Form("PbPbTrackSpectrumByTrigger%d_%d_%d",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.ntrkBins,s.xtrkbins);
    for(int i = 0; i<s.HInTriggers; i++) s.HIUsedByTrigger[i][c] = new TH1D(Form("PbPbUsedTrackSpectrumByTrigger%d_%d_%d",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.ntrkBins,s.xtrkbins);
    s.HIJets[c] = new TH1D(Form("PbPbJetSpectrum_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]),";Leading Jet P_{T} (GeV);#sigma (mb)",s.njetBins,0,s.maxJetBin);
    for(int i = 0; i<s.HInTriggers; i++) s.HIJetsByTrigger[i][c] = new TH1D(Form("PbPbJetSpectrumByTrigger%d_%d_%d",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.njetBins,0,s.maxJetBin);

    s.HI_trk[c] = new TH1D(Form("PbPbTrackSpectrum_trk_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]),";p_{T} (GeV);E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xtrkbins);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIByTrigger_trk[i][c] = new TH1D(Form("PbPbTrackSpectrumByTrigger%d_%d_%d_trk",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.ntrkBins,s.xtrkbins);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIUsedByTrigger_trk[i][c] = new TH1D(Form("PbPbUsedTrackSpectrumByTrigger%d_%d_%d_trk",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.ntrkBins,s.xtrkbins);
    s.HIMaxtrk[c] = new TH1D(Form("PbPbMaxtrkSpectrum_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]),";Leading Track P_{T} (GeV);#sigma (mb)",s.nTrktriggerBins,0,s.maxTrktriggerBin);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIMaxtrkByTrigger[i][c] = new TH1D(Form("PbPbMaxtrkSpectrumByTrigger%d_%d_%d",i,5*s.lowCentBin[c],5*s.highCentBin[c]),"",s.nTrktriggerBins,0,s.maxTrktriggerBin);
  }


  //loading files
  TFile * inFile = TFile::Open("countTracks.root","read");
  for(int i = 0; i<s.nTriggers; i++)
  {
    s.spec[i] = (TH2D*) inFile->Get(Form("spectrum_trigger%d",i));
    s.evtCount[i] = (TH1D*) inFile->Get(Form("evtCount%d",i));
    s.spec[i]->SetDirectory(0);
    s.evtCount[i]->SetDirectory(0);
  }
  for(int i = 0; i<s.nTriggers_trk; i++)
  {
    s.spec_trk[i] = (TH2D*) inFile->Get(Form("spectrum_trigger%d_trk",i));
    s.evtCount_trk[i] = (TH1D*) inFile->Get(Form("evtCount%d_trk",i));
    s.spec_trk[i]->SetDirectory(0);
    s.evtCount_trk[i]->SetDirectory(0);
  }
  s.nVtxMB = (TH1D*) inFile->Get("nVtxMB");
  s.nVtxMB->SetDirectory(0);
  s.nVtxMB_trk = (TH1D*) inFile->Get("nVtxMB_trk");
  s.nVtxMB_trk->SetDirectory(0);
  inFile->Close();
  
  inFile = TFile::Open("RAA_Feb10_testfile.root","read");
  for(int c = 0; c<20; c++)
  {
    for(int i = 0; i<s.HInTriggers; i++)
    {
      s.HIspec[i][c] = (TH2D*) inFile->Get(Form("HI_spectrum_trigger%d_cent%d",i,c));
      s.HIevtCount[i][c] = (TH1D*) inFile->Get(Form("HI_evtCount%d_cent%d",i,c));
      s.HIspec[i][c]->SetDirectory(0);
      s.HIevtCount[i][c]->SetDirectory(0);
    }
    for(int i = 0; i<s.HInTriggers_trk; i++)
    {
      s.HIspec_trk[i][c] = (TH2D*) inFile->Get(Form("HI_spectrum_trigger%d_cent%d_trk",i,c));
      s.HIevtCount_trk[i][c] = (TH1D*) inFile->Get(Form("HI_evtCount%d_cent%d_trk",i,c));
      s.HIspec_trk[i][c]->SetDirectory(0);
      s.HIevtCount_trk[i][c]->SetDirectory(0);
    }
  }
  s.HInVtxMB = (TH1D*) inFile->Get("HI_nVtxMB");
  s.HInVtxMB->SetDirectory(0);
  s.HInVtxMB_trk = (TH1D*) inFile->Get("HI_nVtxMB_trk");
  s.HInVtxMB_trk->SetDirectory(0);
  inFile->Close();

  //calculation of overlaps
  float scale[s.nTriggers];
  float scale_trk[s.nTriggers_trk];
  float HIscale[s.HInTriggers][3];
  float HIscale_trk[s.HInTriggers_trk][3];
  //calculate total number of verticies from MB events
  int nVtx = 0, nVtx_trk = 0, HInVtx = 0, HInVtx_trk = 0;
  for(int i = 1; i<s.nVtxMB->GetSize()+1;i++) nVtx = nVtx+i*s.nVtxMB->GetBinContent(s.nVtxMB->FindBin(i));
  for(int i = 1; i<s.nVtxMB_trk->GetSize()+1;i++) nVtx_trk = nVtx_trk+i*s.nVtxMB_trk->GetBinContent(s.nVtxMB_trk->FindBin(i));
  for(int i = 1; i<s.HInVtxMB->GetSize()+1;i++) HInVtx = nVtx+i*s.HInVtxMB->GetBinContent(s.HInVtxMB->FindBin(i));
  for(int i = 1; i<s.HInVtxMB_trk->GetSize()+1;i++) HInVtx_trk = nVtx_trk+i*s.HInVtxMB_trk->GetBinContent(s.HInVtxMB_trk->FindBin(i));
  
  //*****************************************************************************************************************************************************
  //pp jet triggers
  for(int i = 0; i<s.nTriggers; i++)
  {
    scale[i] = 68/((float)nVtx);//using 68mb as inelastic pp xsection
    for(int j = 0; j<i; j++){
      scale[i] = scale[i]*s.evtCount[j]->Integral(s.evtCount[j]->FindBin(s.triggerOverlapBins[j+1]),s.evtCount[j]->FindBin(s.maxJetBin))/(double)s.evtCount[j+1]->Integral(s.evtCount[j+1]->FindBin(s.triggerOverlapBins[j+1]),s.evtCount[j+1]->FindBin(s.maxJetBin));
    }
    std::cout << scale[i] << std::endl;
    s.spec[i]->Scale(scale[i]);

    //total spectrum
    for(int j = s.evtCount[i]->FindBin(s.triggerBins[i]); j<s.evtCount[i]->FindBin(s.triggerBins[i+1]); j++)
    {
      for(int k = 1; k<s.pp->GetSize()+1; k++)
      {
        s.pp->SetBinContent(k,s.pp->GetBinContent(k)+s.spec[i]->GetBinContent(j,k)); 
        s.pp->SetBinError(k,TMath::Power(TMath::Power(s.pp->GetBinError(k),2)+TMath::Power(s.spec[i]->GetBinError(j,k),2),0.5)); 
        s.ppUsedByTrigger[i]->SetBinContent(k,s.ppUsedByTrigger[i]->GetBinContent(k)+s.spec[i]->GetBinContent(j,k)); 
        s.ppUsedByTrigger[i]->SetBinError(k,TMath::Power(TMath::Power(s.ppUsedByTrigger[i]->GetBinError(k),2)+TMath::Power(s.spec[i]->GetBinError(j,k),2),0.5)); 
      }
      s.ppJets->SetBinContent(j,s.evtCount[i]->GetBinContent(j)*scale[i]);
      s.ppJets->SetBinError(j,s.evtCount[i]->GetBinError(j)*scale[i]);
    }
   
    //spectrum for each trigger
    for(int j = s.evtCount[i]->FindBin(0); j<s.evtCount[i]->FindBin(s.maxJetBin); j++)
    {
      for(int k = 1; k<s.ppByTrigger[i]->GetSize()+1; k++)
      {
        s.ppByTrigger[i]->SetBinContent(k,s.ppByTrigger[i]->GetBinContent(k)+s.spec[i]->GetBinContent(j,k)); 
        s.ppByTrigger[i]->SetBinError(k,TMath::Power(TMath::Power(s.ppByTrigger[i]->GetBinError(k),2)+TMath::Power(s.spec[i]->GetBinError(j,k),2),0.5)); 
      }
      s.ppJetsByTrigger[i]->SetBinContent(j,s.evtCount[i]->GetBinContent(j)*scale[i]);
      s.ppJetsByTrigger[i]->SetBinError(j,s.evtCount[i]->GetBinError(j)*scale[i]);
    }
  }
  
  //pp track triggers
  for(int i = 0; i<s.nTriggers_trk; i++)
  {
    scale_trk[i] = 68/((float)nVtx_trk);//using 68mb as inelastic pp xsection
    for(int j = 0; j<i; j++){
      scale_trk[i] = scale_trk[i]*s.evtCount_trk[j]->Integral(s.evtCount_trk[j]->FindBin(s.triggerOverlapBins_trk[j+1]),s.evtCount_trk[j]->FindBin(s.maxTrktriggerBin))/(double)s.evtCount_trk[j+1]->Integral(s.evtCount_trk[j+1]->FindBin(s.triggerOverlapBins_trk[j+1]),s.evtCount_trk[j+1]->FindBin(s.maxTrktriggerBin));
    }
    std::cout <<"pp scale: Trigger " << i  <<" "<< scale_trk[i] << std::endl;
    s.spec_trk[i]->Scale(scale_trk[i]);

    //total spectrum
    for(int j = s.evtCount_trk[i]->FindBin(s.triggerBins_trk[i]); j<s.evtCount_trk[i]->FindBin(s.triggerBins_trk[i+1]); j++)
    {
      for(int k = 1; k<s.pp_trk->GetSize()+1; k++)
      {
        s.pp_trk->SetBinContent(k,s.pp_trk->GetBinContent(k)+s.spec_trk[i]->GetBinContent(j,k)); 
        s.pp_trk->SetBinError(k,TMath::Power(TMath::Power(s.pp_trk->GetBinError(k),2)+TMath::Power(s.spec_trk[i]->GetBinError(j,k),2),0.5)); 
        s.ppUsedByTrigger_trk[i]->SetBinContent(k,s.ppUsedByTrigger_trk[i]->GetBinContent(k)+s.spec_trk[i]->GetBinContent(j,k)); 
        s.ppUsedByTrigger_trk[i]->SetBinError(k,TMath::Power(TMath::Power(s.ppUsedByTrigger_trk[i]->GetBinError(k),2)+TMath::Power(s.spec_trk[i]->GetBinError(j,k),2),0.5)); 
      }
      s.ppMaxtrk->SetBinContent(j,s.evtCount_trk[i]->GetBinContent(j)*scale_trk[i]);
      s.ppMaxtrk->SetBinError(j,s.evtCount_trk[i]->GetBinError(j)*scale_trk[i]);
    }
    
    //spec_trktrum for each trigger
    for(int j = s.evtCount_trk[i]->FindBin(0); j<s.evtCount_trk[i]->FindBin(s.maxTrktriggerBin); j++)
    {
      for(int k = 1; k<s.ppByTrigger_trk[i]->GetSize()+1; k++)
      {
        s.ppByTrigger_trk[i]->SetBinContent(k,s.ppByTrigger_trk[i]->GetBinContent(k)+s.spec_trk[i]->GetBinContent(j,k)); 
        s.ppByTrigger_trk[i]->SetBinError(k,TMath::Power(TMath::Power(s.ppByTrigger_trk[i]->GetBinError(k),2)+TMath::Power(s.spec_trk[i]->GetBinError(j,k),2),0.5)); 
      }
      s.ppMaxtrkByTrigger[i]->SetBinContent(j,s.evtCount_trk[i]->GetBinContent(j)*scale_trk[i]);
      s.ppMaxtrkByTrigger[i]->SetBinError(j,s.evtCount_trk[i]->GetBinError(j)*scale_trk[i]);
    }
  }

  TH1D * tempEvtCount[s.HInTriggers][3];
  int combinationCentUpperBoundary[3] = {5,9,19};
  int combinationCentLowerBoundary[3] = {0,6,10};
  for(int i = 0; i<s.HInTriggers; i++)
  {
    for(int m = 0; m<3; m++){
      HIscale[i][m] = 68/((float)HInVtx);//using 68mb as inelastic pp xsection
      tempEvtCount[i][m] = (TH1D*)s.HIevtCount[i][0]->Clone(Form("HItempEvtCount%d",m));
      tempEvtCount[i][m]->Reset();
      for(int c = combinationCentUpperBoundary[m]; c>=combinationCentLowerBoundary[m]; c--) tempEvtCount[i][m]->Add(s.HIevtCount[i][c]);
      for(int j = 0; j<i; j++){
        HIscale[i][m] = HIscale[i][m]*tempEvtCount[j][m]->Integral(tempEvtCount[j][m]->FindBin(s.HItriggerOverlapBins[j+1]),tempEvtCount[j][m]->FindBin(s.maxJetBin))/(double)tempEvtCount[j+1][m]->Integral(tempEvtCount[j+1][m]->FindBin(s.HItriggerOverlapBins[j+1]),tempEvtCount[j+1][m]->FindBin(s.maxJetBin));
      }
      std::cout <<"PbPb scale: Trigger and cent region "<< i<<" "<<m<<" "<<HIscale[i][m] << std::endl;
    }

    for(int c = 0; c<20; c++){
      if(c<6){
        s.HIspec[i][c]->Scale(HIscale[i][0]);
        s.HIevtCount[i][c]->Scale(HIscale[i][0]);
      }
      else if(c<10){
        s.HIspec[i][c]->Scale(HIscale[i][1]);
        s.HIevtCount[i][c]->Scale(HIscale[i][1]);
      }else{
        s.HIspec[i][c]->Scale(HIscale[i][2]);
        s.HIevtCount[i][c]->Scale(HIscale[i][2]);
      } 
    }
    for(int c = 20; c<s.nCentBins; c++){
      s.HIspec[i][c] = (TH2D*)s.HIspec[i][0]->Clone(Form("HI_spectrum_trigger%d_cent%d",i,c));
      s.HIspec[i][c]->Reset();
      s.HIevtCount[i][c] = (TH1D*)s.HIevtCount[i][0]->Clone(Form("HI_evtCount%d_cent%d",i,c));
      s.HIevtCount[i][c]->Reset();
      for(int cc = s.lowCentBin[c]; cc<s.highCentBin[c]; cc++){
        s.HIspec[i][c]->Add(s.HIspec[i][cc]);
        s.HIevtCount[i][c]->Add(s.HIevtCount[i][cc]);
      }
    }
     
    //total spectrum
    for(int c = 0; c<s.nCentBins; c++){
      for(int j = s.HIevtCount[i][c]->FindBin(s.HItriggerBins[i]); j<s.HIevtCount[i][c]->FindBin(s.HItriggerBins[i+1]); j++)
      {
        for(int k = 1; k<s.HI[c]->GetSize()+1; k++)
        {
          s.HI[c]->SetBinContent(k,s.HI[c]->GetBinContent(k)+s.HIspec[i][c]->GetBinContent(j,k)); 
          s.HI[c]->SetBinError(k,TMath::Power(TMath::Power(s.HI[c]->GetBinError(k),2)+TMath::Power(s.HIspec[i][c]->GetBinError(j,k),2),0.5)); 
          s.HIUsedByTrigger[i][c]->SetBinContent(k,s.HIUsedByTrigger[i][c]->GetBinContent(k)+s.HIspec[i][c]->GetBinContent(j,k)); 
          s.HIUsedByTrigger[i][c]->SetBinError(k,TMath::Power(TMath::Power(s.HIUsedByTrigger[i][c]->GetBinError(k),2)+TMath::Power(s.HIspec[i][c]->GetBinError(j,k),2),0.5)); 
        }
        s.HIJets[c]->SetBinContent(j,s.HIevtCount[i][c]->GetBinContent(j));
        s.HIJets[c]->SetBinError(j,s.HIevtCount[i][c]->GetBinError(j));
      }
     
      //spectrum for each trigger
      for(int j = s.HIevtCount[i][c]->FindBin(0); j<s.HIevtCount[i][c]->FindBin(s.maxJetBin); j++)
      {
        for(int k = 1; k<s.HIByTrigger[i][c]->GetSize()+1; k++)
        {
          s.HIByTrigger[i][c]->SetBinContent(k,s.HIByTrigger[i][c]->GetBinContent(k)+s.HIspec[i][c]->GetBinContent(j,k)); 
          s.HIByTrigger[i][c]->SetBinError(k,TMath::Power(TMath::Power(s.HIByTrigger[i][c]->GetBinError(k),2)+TMath::Power(s.HIspec[i][c]->GetBinError(j,k),2),0.5)); 
        }
        s.HIJetsByTrigger[i][c]->SetBinContent(j,s.HIevtCount[i][c]->GetBinContent(j));
        s.HIJetsByTrigger[i][c]->SetBinError(j,s.HIevtCount[i][c]->GetBinError(j));
      }
    }//cent loop closed
  }

  TH1D * tempEvtCount_trk[s.HInTriggers_trk][3];
  int combinationCentUpperBoundary_trk[3] = {5,9,19};
  int combinationCentLowerBoundary_trk[3] = {0,6,10};
  for(int i = 0; i<s.HInTriggers_trk; i++)
  {
    for(int m = 0; m<3; m++){
      HIscale_trk[i][m] = 68/((float)HInVtx_trk);//using 68mb as inelastic pp xsection
      tempEvtCount_trk[i][m] = (TH1D*)s.HIevtCount_trk[i][0]->Clone(Form("HItempEvtCount%d_trk",m));
      tempEvtCount_trk[i][m]->Reset();
      for(int c = combinationCentUpperBoundary_trk[m]; c>=combinationCentLowerBoundary_trk[m]; c--) tempEvtCount_trk[i][m]->Add(s.HIevtCount_trk[i][c]);
      for(int j = 0; j<i; j++){
        HIscale_trk[i][m] = HIscale_trk[i][m]*tempEvtCount_trk[j][m]->Integral(tempEvtCount_trk[j][m]->FindBin(s.HItriggerOverlapBins_trk[j+1]),tempEvtCount_trk[j][m]->FindBin(s.maxTrktriggerBin))/(double)tempEvtCount_trk[j+1][m]->Integral(tempEvtCount_trk[j+1][m]->FindBin(s.HItriggerOverlapBins_trk[j+1]),tempEvtCount_trk[j+1][m]->FindBin(s.maxTrktriggerBin));
      }
      std::cout <<"PbPb scale: Trigger and cent region "<< i<<" "<<m<<" "<<HIscale_trk[i][m] << std::endl;
    }

    for(int c = 0; c<20; c++){
      if(c<6){
        s.HIspec_trk[i][c]->Scale(HIscale_trk[i][0]);
        s.HIevtCount_trk[i][c]->Scale(HIscale_trk[i][0]);
      }
      else if(c<10){
        s.HIspec_trk[i][c]->Scale(HIscale_trk[i][1]);
        s.HIevtCount_trk[i][c]->Scale(HIscale_trk[i][1]);
      }else{
        s.HIspec_trk[i][c]->Scale(HIscale_trk[i][2]);
        s.HIevtCount_trk[i][c]->Scale(HIscale_trk[i][2]);
      } 
    }
    for(int c = 20; c<s.nCentBins; c++){
      s.HIspec_trk[i][c] = (TH2D*)s.HIspec_trk[i][0]->Clone(Form("HI_spectrum_trigger%d_cent%d_trk",i,c));
      s.HIspec_trk[i][c]->Reset();
      s.HIevtCount_trk[i][c] = (TH1D*)s.HIevtCount_trk[i][0]->Clone(Form("HI_evtCount%d_cent%d_trk",i,c));
      s.HIevtCount_trk[i][c]->Reset();
      for(int cc = s.lowCentBin[c]; cc<s.highCentBin[c]; cc++){
        s.HIspec_trk[i][c]->Add(s.HIspec_trk[i][cc]);
        s.HIevtCount_trk[i][c]->Add(s.HIevtCount_trk[i][cc]);
      }
    }
     
    //total spectrum
    for(int c = 0; c<s.nCentBins; c++){
      for(int j = s.HIevtCount_trk[i][c]->FindBin(s.HItriggerBins_trk[i]); j<s.HIevtCount_trk[i][c]->FindBin(s.HItriggerBins_trk[i+1]); j++)
      {
        for(int k = 1; k<s.HI_trk[c]->GetSize()+1; k++)
        {
          s.HI_trk[c]->SetBinContent(k,s.HI_trk[c]->GetBinContent(k)+s.HIspec_trk[i][c]->GetBinContent(j,k)); 
          s.HI_trk[c]->SetBinError(k,TMath::Power(TMath::Power(s.HI_trk[c]->GetBinError(k),2)+TMath::Power(s.HIspec_trk[i][c]->GetBinError(j,k),2),0.5)); 
          s.HIUsedByTrigger_trk[i][c]->SetBinContent(k,s.HIUsedByTrigger_trk[i][c]->GetBinContent(k)+s.HIspec_trk[i][c]->GetBinContent(j,k)); 
          s.HIUsedByTrigger_trk[i][c]->SetBinError(k,TMath::Power(TMath::Power(s.HIUsedByTrigger_trk[i][c]->GetBinError(k),2)+TMath::Power(s.HIspec_trk[i][c]->GetBinError(j,k),2),0.5)); 
        }
        s.HIMaxtrk[c]->SetBinContent(j,s.HIevtCount_trk[i][c]->GetBinContent(j));
        s.HIMaxtrk[c]->SetBinError(j,s.HIevtCount_trk[i][c]->GetBinError(j));
      }
     
      //spectrum for each trigger
      for(int j = s.HIevtCount_trk[i][c]->FindBin(0); j<s.HIevtCount_trk[i][c]->FindBin(s.maxTrktriggerBin); j++)
      {
        for(int k = 1; k<s.HIByTrigger_trk[i][c]->GetSize()+1; k++)
        {
          s.HIByTrigger_trk[i][c]->SetBinContent(k,s.HIByTrigger_trk[i][c]->GetBinContent(k)+s.HIspec_trk[i][c]->GetBinContent(j,k)); 
          s.HIByTrigger_trk[i][c]->SetBinError(k,TMath::Power(TMath::Power(s.HIByTrigger_trk[i][c]->GetBinError(k),2)+TMath::Power(s.HIspec_trk[i][c]->GetBinError(j,k),2),0.5)); 
        }
        s.HIMaxtrkByTrigger[i][c]->SetBinContent(j,s.HIevtCount_trk[i][c]->GetBinContent(j));
        s.HIMaxtrkByTrigger[i][c]->SetBinError(j,s.HIevtCount_trk[i][c]->GetBinError(j));
      }
    }//cent loop closed
  }

  //*********************************************************************************************************************************************88
  //Divide by bin width and jacobian stuff
  for(int i = 1; i<s.pp->GetSize()+1; i++)
  {
    s.pp->SetBinContent(i,s.pp->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
    s.pp->SetBinError(i,s.pp->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    s.pp_trk->SetBinContent(i,s.pp_trk->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
    s.pp_trk->SetBinError(i,s.pp_trk->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    for(int j = 0; j<s.nTriggers; j++)
    {
      s.ppByTrigger[j]->SetBinContent(i,s.ppByTrigger[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.ppByTrigger[j]->SetBinError(i,s.ppByTrigger[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      s.ppUsedByTrigger[j]->SetBinContent(i,s.ppUsedByTrigger[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.ppUsedByTrigger[j]->SetBinError(i,s.ppUsedByTrigger[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    } 
    for(int j = 0; j<s.nTriggers_trk; j++)
    {
      s.ppByTrigger_trk[j]->SetBinContent(i,s.ppByTrigger_trk[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.ppByTrigger_trk[j]->SetBinError(i,s.ppByTrigger_trk[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      s.ppUsedByTrigger_trk[j]->SetBinContent(i,s.ppUsedByTrigger_trk[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.ppUsedByTrigger_trk[j]->SetBinError(i,s.ppUsedByTrigger_trk[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    } 
  }

  for(int c = 0; c<s.nCentBins; c++)
  {
    for(int i = 1; i<s.HI[c]->GetSize()+1; i++)
    {
      s.HI[c]->SetBinContent(i,s.HI[c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.HI[c]->SetBinError(i,s.HI[c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      s.HI_trk[c]->SetBinContent(i,s.HI_trk[c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      s.HI_trk[c]->SetBinError(i,s.HI_trk[c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      for(int j = 0; j<s.HInTriggers; j++)
      {
        s.HIByTrigger[j][c]->SetBinContent(i,s.HIByTrigger[j][c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
        s.HIByTrigger[j][c]->SetBinError(i,s.HIByTrigger[j][c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
        s.HIUsedByTrigger[j][c]->SetBinContent(i,s.HIUsedByTrigger[j][c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
        s.HIUsedByTrigger[j][c]->SetBinError(i,s.HIUsedByTrigger[j][c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      } 
      for(int j = 0; j<s.HInTriggers_trk; j++)
      {
        s.HIByTrigger_trk[j][c]->SetBinContent(i,s.HIByTrigger_trk[j][c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
        s.HIByTrigger_trk[j][c]->SetBinError(i,s.HIByTrigger_trk[j][c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
        s.HIUsedByTrigger_trk[j][c]->SetBinContent(i,s.HIUsedByTrigger_trk[j][c]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
        s.HIUsedByTrigger_trk[j][c]->SetBinError(i,s.HIUsedByTrigger_trk[j][c]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      } 
    }
  }
 
  //*************************OUTPUT**********************************************************
  TFile * outF = TFile::Open("Spectra.root","recreate");
  s.pp->Write();
  s.pp_trk->Write();
  s.pp_perMBTrigger = (TH1D*)s.pp->Clone("pp_perMBTrigger");
  s.pp_perMBTrigger->Scale(1/scale[0]);
  s.pp_perMBTrigger->Write();
  s.pp_perMBTrigger_trk = (TH1D*)s.pp_trk->Clone("pp_perMBTrigger_trk");
  s.pp_perMBTrigger_trk->Scale(1/scale_trk[0]);
  s.pp_perMBTrigger_trk->Write();
  
  s.ppJets->SetMarkerSize(0);
  s.ppJets->Write();
  s.ppMaxtrk->SetMarkerSize(0);
  s.ppMaxtrk->Write();
  for(int j = 0; j<s.nTriggers; j++)
  {
    s.ppByTrigger[j]->Write();
    s.ppUsedByTrigger[j]->Write();
    s.ppJetsByTrigger[j]->Write();
  }
  for(int j = 0; j<s.nTriggers_trk; j++)
  {
    s.ppByTrigger_trk[j]->Write();
    s.ppUsedByTrigger_trk[j]->Write();
    s.ppMaxtrkByTrigger[j]->Write();
  }
  for(int c = 0; c<s.nCentBins; c++){
    s.HI[c]->Write();
    s.HI_trk[c]->Write();
    s.HI_perMBTrigger[c] = (TH1D*)s.HI[c]->Clone(Form("PbPb_perMBTrigger_%d_%d",5*s.lowCentBin[c],5*s.highCentBin[c]));
    s.HI_perMBTrigger[c]->Scale(1/HIscale[0][0]);
    s.HI_perMBTrigger[c]->Write();
    s.HI_perMBTrigger_trk[c] = (TH1D*)s.HI_trk[c]->Clone(Form("PbPb_perMBTrigger_%d_%d_trk",5*s.lowCentBin[c],5*s.highCentBin[c]));
    s.HI_perMBTrigger_trk[c]->Scale(1/HIscale_trk[0][0]);
    s.HI_perMBTrigger_trk[c]->Write();
    
    s.HIJets[c]->Write();
    s.HIMaxtrk[c]->Write();
    for(int j = 0; j<s.HInTriggers; j++){
      s.HIByTrigger[j][c]->Write();
      s.HIUsedByTrigger[j][c]->Write();
      s.HIJetsByTrigger[j][c]->Write();
    }
    for(int j = 0; j<s.HInTriggers_trk; j++){
      s.HIByTrigger_trk[j][c]->Write();
      s.HIUsedByTrigger_trk[j][c]->Write();
      s.HIMaxtrkByTrigger[j][c]->Write();
    }
  }
  outF->Close();
 
  makePlotsPbPb(s);
  makePlotsPP(s);
}
