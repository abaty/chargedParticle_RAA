#ifndef PRETTYPLOT
#define PRETTYPLOT

#include "TAttMarker.h"
#include "TLine.h"
#include "TAttLine.h"
#include "TColor.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TAttAxis.h"
#include "TAttText.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TAttFill.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "Settings.h"
#include "hyperon_check/hyperonCorrection.C"
#include "TFrame.h"


double Quad(double a, double b)
{
  return TMath::Power(TMath::Power(a,2) + TMath::Power(b,2),0.5);
}

void prettyPlotting(Settings s){
  TFile * inputPlots = TFile::Open("Spectra.root","Update");
  TH1D * h[s.nCentBins];

  for(int c = 0; c<s.nCentBins; c++) h[c] = (TH1D*)inputPlots->Get(Form("RAA_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
  for(int c = 0; c<s.nCentBins; c++){
     h[c]->SetDirectory(0);
     h[c]->GetYaxis()->CenterTitle(true);
     h[c]->GetYaxis()->SetTitle("R_{AA}");
     h[c]->GetXaxis()->SetTitle("p_{T} (GeV)");
     h[c]->GetXaxis()->CenterTitle(true);
     s.RAA_totSyst[c] = (TH1D*)h[c]->Clone(Form("RAA_totSyst_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
     s.RAA_totSyst[c]->Reset();
     s.RAA_totSyst[c]->SetDirectory(inputPlots);
  }
  TH1D * hyperonPbPb = (TH1D*)h[0]->Clone("hyperonPbPb");
  TH1D * hyperonpp = (TH1D*)h[0]->Clone("hyperonpp");
  returnHyperonCorrection(0,hyperonPbPb,"hyperon_check/");
  returnHyperonCorrection(0,hyperonpp,"hyperon_check/");//need to change to pp
  hyperonPbPb->Print("All");
  hyperonpp->Print("All");
  for(int c = 0; c<s.nCentBins; c++){
    h[c]->Multiply(hyperonPbPb);
    h[c]->Divide(hyperonpp);
  }

  setTDRStyle();
  TLine * line1;
  TLatex * tex = new TLatex(0.1,0.1,"cent");
  TLatex * tex2 = new TLatex(0.1,0.1,"cent");
  TBox* bLumi = new TBox(0.1,0.1,0.2,0.2); 
  TBox* b[s.ntrkBins];
  for(int i = 0; i<s.ntrkBins; i++) b[i] = new TBox(0.1,0.1,0.2,0.2); 
  
 
  int W = 800;
  int H = 700;
  int H_ref = 700; 
  int W_ref = 800;
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.15*W_ref;
  float R = 0.04*W_ref;
    
  TCanvas* canv = new TCanvas("RAA","RAA",50,50,W,H);
  canv->SetLogx();
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
 
  gStyle->SetErrorX(0);

  for(int c = 0; c<s.nCentBins; c++){
    //adding up uncertainties
    for(int i = 1; i<s.RAA_totSyst[c]->GetSize()-1; i++){
      s.RAA_totSyst[c]->SetBinContent(i,0);
    
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));//4% difference in data/MC (PbPb)
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));//4% difference in data/MC (pp)
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.02));//2% for nonclosure PbPb
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.02));//2% for nonclosure pp 
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_HInormSyst[c]->GetBinContent(i)));//add in PbPb normalization uncert
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_normSyst->GetBinContent(i)));//add in pp normalization uncert
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),hyperonPbPb->GetBinContent(i)-1));//PbPb hyperon study
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),hyperonpp->GetBinContent(i)-1));//pp hyperon study
    }
    s.RAA_totSyst[c]->Write();

    //plotting
    canv->Clear();
    h[c]->GetXaxis()->SetRangeUser(0.7,350);
    h[c]->GetYaxis()->SetRangeUser(0,1.4);
    h[c]->Draw();
  
    float lumiUncert = 0.1;
    bLumi->SetFillColor(kGray);
    bLumi->SetLineWidth(0);
    bLumi->DrawBox(0.9,1-lumiUncert,1.5,1+lumiUncert);
  
    line1 = new TLine(h[c]->GetXaxis()->GetBinLowEdge(3),1,h[c]->GetXaxis()->GetBinUpEdge(h[c]->GetSize()-2),1);
    line1->SetLineWidth(2);
    line1->SetLineStyle(2);
    line1->Draw("same");
  
    tex2->DrawLatex(0.9,0.1,Form("%d-%d%s",5*s.lowCentBin[c],5*s.highCentBin[c],"%"));
    tex->SetTextFont(42);
    tex->SetTextSize(lumiTextSize*0.08);
    tex->DrawLatex(1.8,1.03,"T_{AA} and lumi. uncertainty");
    tex->DrawLatex(1.8,0.93,"|#eta|<1");
  
    for(int i = 1; i<h[c]->GetSize()-1; i++){
      if(i<3) continue;
      b[i-1]->SetFillColor(kOrange);
      b[i-1]->SetX1(h[c]->GetXaxis()->GetBinLowEdge(i));
      b[i-1]->SetX2(h[c]->GetXaxis()->GetBinUpEdge(i));
      b[i-1]->SetY1((h[c]->GetBinContent(i))*(1-s.RAA_totSyst[c]->GetBinContent(i)));
      b[i-1]->SetY2(h[c]->GetBinContent(i)*(1+s.RAA_totSyst[c]->GetBinContent(i)));
      b[i-1]->Draw("same");
    }
    h[c]->Draw("same");
  
    int iPeriod = 0;
    lumi_sqrtS = "25.8 pb^{-1} (5.02 TeV pp) + 404 #mub^{-1} (5.02 TeV PbPb)";
    writeExtraText = true;  
    extraText  = "Preliminary";
    CMS_lumi( canv, iPeriod, 11 );
 
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();    
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
    delete line1;
  }
  inputPlots->Close();
  
  return;
}
#endif
