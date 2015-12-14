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

void makeSpectrum()
{
  Settings s;
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
//input hist  
  TH2D * spec[s.nTriggers];
  TH1D * evtCount[s.nTriggers];
  TH1D * nVtxMB;

//output hist
  TH1D * pp = new TH1D("ppTrackSpectrum",";p_{T} (GeV);E#frac{d^{3}#sigma}{d^{3}p} (mb/GeV^{2})",s.ntrkBins,s.xtrkbins);
  TH1D * ppByTrigger[s.nTriggers];
  TH1D * ppUsedByTrigger[s.nTriggers];
  for(int i = 0; i<s.nTriggers; i++) ppByTrigger[i] = new TH1D(Form("ppTrackSpectrumByTrigger%d",i),"",s.ntrkBins,s.xtrkbins);
  for(int i = 0; i<s.nTriggers; i++) ppUsedByTrigger[i] = new TH1D(Form("ppUsedTrackSpectrumByTrigger%d",i),"",s.ntrkBins,s.xtrkbins);
  TH1D * ppJets = new TH1D("ppJetSpectrum",";Leading Jet P_{T} (GeV);#sigma (mb)",s.njetBins,0,s.maxJetBin);
  TH1D * ppJetsByTrigger[s.nTriggers];
  for(int i = 0; i<s.nTriggers; i++) ppJetsByTrigger[i] = new TH1D(Form("ppJetSpectrumByTrigger%d",i),"",s.njetBins,0,s.maxJetBin);
 

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
      scale[i] = scale[i]*evtCount[j]->Integral(evtCount[j]->FindBin(s.triggerOverlapBins[j+1]),evtCount[j]->FindBin(s.maxJetBin))/(double)evtCount[j+1]->Integral(evtCount[j+1]->FindBin(s.triggerOverlapBins[j+1]),evtCount[j+1]->FindBin(s.maxJetBin));
    }
    std::cout << scale[i] << std::endl;
    spec[i]->Scale(scale[i]);

    //total spectrum
    for(int j = evtCount[i]->FindBin(s.triggerBins[i]); j<evtCount[i]->FindBin(s.triggerBins[i+1]); j++)
    {
      for(int k = 1; k<pp->GetSize()+1; k++)
      {
        pp->SetBinContent(k,pp->GetBinContent(k)+spec[i]->GetBinContent(j,k)); 
        pp->SetBinError(k,TMath::Power(TMath::Power(pp->GetBinError(k),2)+TMath::Power(spec[i]->GetBinError(j,k),2),0.5)); 
        ppUsedByTrigger[i]->SetBinContent(k,ppUsedByTrigger[i]->GetBinContent(k)+spec[i]->GetBinContent(j,k)); 
        ppUsedByTrigger[i]->SetBinError(k,TMath::Power(TMath::Power(ppUsedByTrigger[i]->GetBinError(k),2)+TMath::Power(spec[i]->GetBinError(j,k),2),0.5)); 
      }
      ppJets->SetBinContent(j,evtCount[i]->GetBinContent(j)*scale[i]);
      ppJets->SetBinError(j,evtCount[i]->GetBinError(j)*scale[i]);
    }
    
    //spectrum for each trigger
    for(int j = evtCount[i]->FindBin(0); j<evtCount[i]->FindBin(s.maxJetBin); j++)
    {
      for(int k = 1; k<ppByTrigger[i]->GetSize()+1; k++)
      {
        ppByTrigger[i]->SetBinContent(k,ppByTrigger[i]->GetBinContent(k)+spec[i]->GetBinContent(j,k)); 
        ppByTrigger[i]->SetBinError(k,TMath::Power(TMath::Power(ppByTrigger[i]->GetBinError(k),2)+TMath::Power(spec[i]->GetBinError(j,k),2),0.5)); 
      }
      ppJetsByTrigger[i]->SetBinContent(j,evtCount[i]->GetBinContent(j)*scale[i]);
      ppJetsByTrigger[i]->SetBinError(j,evtCount[i]->GetBinError(j)*scale[i]);
    }
  }

  //Divide by bin width and jacobian stuff
  for(int i = 1; i<pp->GetSize()+1; i++)
  {
    pp->SetBinContent(i,pp->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
    pp->SetBinError(i,pp->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    for(int j = 0; j<s.nTriggers; j++)
    {
      ppByTrigger[j]->SetBinContent(i,ppByTrigger[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      ppByTrigger[j]->SetBinError(i,ppByTrigger[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
      ppUsedByTrigger[j]->SetBinContent(i,ppUsedByTrigger[j]->GetBinContent(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1]))); 
      ppUsedByTrigger[j]->SetBinError(i,ppUsedByTrigger[j]->GetBinError(i)/(4*TMath::Pi()*(s.xtrkbins[i]-s.xtrkbins[i-1])));
    } 
  } 
  
  TFile * outF = TFile::Open("ppSpectrum.root","recreate");
  pp->Write();
  TH1D * pp_perMBTrigger = (TH1D*)pp->Clone("pp_perMBTrigger");
  pp_perMBTrigger->Scale(1/scale[0]);
  pp_perMBTrigger->Write();
  ppJets->SetMarkerSize(0);
  ppJets->Write();
  for(int j = 0; j<s.nTriggers; j++)
  {
    ppByTrigger[j]->SetLineColor(j+1);
    ppByTrigger[j]->SetLineWidth(1);
    ppByTrigger[j]->SetMarkerColor(j+1);
    ppByTrigger[j]->SetFillColor(j+1);
    ppByTrigger[j]->Write();
    ppUsedByTrigger[j]->SetLineColor(kBlack);
    ppUsedByTrigger[j]->SetLineWidth(2);
    ppUsedByTrigger[j]->SetMarkerColor(j+1);
    ppUsedByTrigger[j]->SetMarkerSize(0);
    ppUsedByTrigger[j]->SetFillColor(j+1);
    ppUsedByTrigger[j]->Write();
    ppJetsByTrigger[j]->SetLineColor(j+1);
    ppJetsByTrigger[j]->SetLineWidth(1);
    ppJetsByTrigger[j]->SetMarkerColor(j+1);
    ppJetsByTrigger[j]->SetMarkerSize(0.8);
    ppJetsByTrigger[j]->SetFillColor(j+1);
    ppJetsByTrigger[j]->Write();
  }
  outF->Close();

//******************************************************************************************************************
//***************************************************PLOTS**********************************************************
//******************************************************************************************************************
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  c1->SetLogy();
  ppJets->Scale(100);
  float Ymin = 0.00000000001;
  float Ymax = 10;
  ppJets->GetYaxis()->SetRangeUser(Ymin,Ymax);
  ppJets->Draw("h");
  for(int i = 0; i<s.nTriggers; i++) ppJetsByTrigger[i]->Draw("same");
  TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(ppJets,"Jet Spectrum (x100)","l");
  leg->AddEntry(ppJetsByTrigger[0],"MB (I)","p");
  leg->AddEntry(ppJetsByTrigger[1],"jet 40 (II)","p");
  leg->AddEntry(ppJetsByTrigger[2],"jet 60 (III)","p");
  leg->AddEntry(ppJetsByTrigger[3],"jet 80 (IV)","p");
  leg->AddEntry((TObject*)0,"ak4Calo Jets, |#eta|<2","");
  leg->Draw("same");
  c1->SaveAs("plots/ppJets_FullSpectrum.png"); 
  c1->SaveAs("plots/ppJets_FullSpectrum.pdf"); 

  ppJets->GetXaxis()->SetRangeUser(20,200);
  Ymin = 0.0000001;
  Ymax = 100;
  ppJets->GetYaxis()->SetRangeUser(0.0000001,100);
  ppJets->Draw("h");
  for(int i = 0; i<s.nTriggers; i++) ppJetsByTrigger[i]->Draw("same");
  leg->Draw("same"); 
  TLine * l[3];
  for(int i = 0; i<3; i++) 
  {
    l[i] = new TLine(s.triggerBins[i+1],Ymin,s.triggerBins[i+1],Ymax); 
    l[i]->SetLineWidth(2);
    l[i]->SetLineStyle(2);
    l[i]->SetLineColor(1);
    l[i]->Draw("same");
  }
  TLatex * lat = new TLatex(1,1,"test");
  lat->DrawLatex(45,Ymin*3,"I");
  lat->DrawLatex(65,Ymin*3,"II");
  lat->DrawLatex(85,Ymin*3,"III");
  lat->DrawLatex(105,Ymin*3,"IV");
  c1->SaveAs("plots/ppJets_FullSpectrum_XZoom.png"); 
  c1->SaveAs("plots/ppJets_FullSpectrum_XZoom.pdf"); 

  c1->SetLogx();
  pp->SetMarkerSize(0.8);
  pp->Draw();
  ppUsedByTrigger[0]->SetFillColor(kGray);
  ppUsedByTrigger[3]->SetFillColor(kCyan+2);
  ppUsedByTrigger[1]->Add(ppUsedByTrigger[0]);
  ppUsedByTrigger[2]->Add(ppUsedByTrigger[1]);
  ppUsedByTrigger[3]->Add(ppUsedByTrigger[2]);
  for(int i = 0; i<s.nTriggers; i++) ppUsedByTrigger[s.nTriggers-1-i]->Draw("HIST same");
  pp->Draw("sameaxis");
  pp->Draw("same");
  leg->Clear();
  leg->AddEntry(pp,"pp track Spectrum","p");
  leg->AddEntry(ppUsedByTrigger[0],"MB trigger","f");
  leg->AddEntry(ppUsedByTrigger[1],"Jet40 trigger","f");
  leg->AddEntry(ppUsedByTrigger[2],"Jet60 trigger","f");
  leg->AddEntry(ppUsedByTrigger[3],"Jet80 trigger","f");
  leg->AddEntry((TObject*)0,"|#eta|<1","");
  leg->Draw("same");
   
  c1->SaveAs("plots/ppTrack_FullSpectrum.png");
  c1->SaveAs("plots/ppTrack_FullSpectrum.pdf");
}
