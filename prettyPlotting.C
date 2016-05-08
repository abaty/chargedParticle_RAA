#ifndef PRETTYPLOT
#define PRETTYPLOT

#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TAttPad.h"
#include "TGraph.h"
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
#include "TAttCanvas.h"
#include "TBox.h"
#include "TAttFill.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "Settings.h"
#include "hyperon_check/hyperonCorrection.C"
#include "TFrame.h"
#include <vector>
#include <string>
#include <fstream>

void get276RAA(TCanvas * c276, Settings s, int centralityBin);
void gettheoryRAA(TCanvas * c_th, Settings s, int centralityBin, std::string saveString);

double Quad(double a, double b)
{
  return TMath::Power(TMath::Power(a,2) + TMath::Power(b,2),0.5);
}

void prettyPlotting(Settings s){
  TFile * inputPlots = TFile::Open("Spectra.root","Update");
  TH1D * h[s.nCentBins];
  TH1D * pbpbSpec[s.nCentBins];
  TH1D * ppSpec;

  for(int c = 0; c<s.nCentBins; c++) h[c] = (TH1D*)inputPlots->Get(Form("RAA_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
  for(int c = 0; c<s.nCentBins; c++) pbpbSpec[c] = (TH1D*)inputPlots->Get(Form("PbPbTrackSpectrum_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
  ppSpec = (TH1D*)inputPlots->Get(Form("pp_NotperMBTrigger"));
  for(int c = 0; c<s.nCentBins; c++){
     h[c]->SetDirectory(0);
     h[c]->GetYaxis()->CenterTitle(true);
     h[c]->GetYaxis()->SetTitle("R_{AA}");
     h[c]->GetXaxis()->SetTitle("p_{T} (GeV)");
     h[c]->GetXaxis()->CenterTitle(true);
     s.RAA_totSyst[c] = (TH1D*)h[c]->Clone(Form("RAA_totSyst_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
     s.RAA_totSyst[c]->Reset();
     s.RAA_totSyst[c]->SetDirectory(inputPlots);
     s.PbPb_totSyst[c] = (TH1D*)h[c]->Clone(Form("RAA_totSyst_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
     s.PbPb_totSyst[c]->Reset();
     s.PbPb_totSyst[c]->SetDirectory(inputPlots);
     if(c==0){
       s.pp_totSyst = (TH1D*)h[c]->Clone(Form("RAA_totSyst_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
       s.pp_totSyst->Reset();
       s.pp_totSyst->SetDirectory(inputPlots);
     }
  }
  TH1D * hyperonPbPb = (TH1D*)h[0]->Clone("hyperonPbPb");
  TH1D * hyperonpp = (TH1D*)h[0]->Clone("hyperonpp");
  returnHyperonCorrection(0,hyperonPbPb,"hyperon_check/");
  returnHyperonCorrection(1,hyperonpp,"hyperon_check/");//need to change to pp
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
  TBox* bLumi = new TBox(0.1,0.1,0.15,0.2); 
  TBox* bTAA = new TBox(0.15,0.1,0.2,0.2); 
  TBox* b[s.ntrkBins];
  for(int i = 0; i<s.ntrkBins; i++) b[i] = new TBox(0.1,0.1,0.2,0.2); 
  
 
  int W = 800;
  int H = 700;//700
  int H_ref = 700;//700
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
      s.PbPb_totSyst[c]->SetBinContent(i,0);
      if(c==0) s.pp_totSyst->SetBinContent(i,0);
    
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.05));//4% difference in data/MC (PbPb)
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.05));//4% difference in data/MC (PbPb)
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));//4% difference in data/MC (pp)
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.04));//4% difference in data/MC (pp)
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.05));//5% for nonclosure PbPb
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.05));//5% for nonclosure (PbPb)
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//1% for nonclosure pp 
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.01));//1% for nonclosure in pp
      
      //!this sytematic is largely bullshit since we don't know the data fake rate!
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.03));//3% for MC-based fake rate PbPb
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.03));//3% difference in data/MC (PbPb)
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//for MC-based fake rate pp 
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.01));//for MC-based faked rate pp
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//1% resolution for not unfolding
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.01));//1% resolution for not unfolding
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//1% resolution for not unfolding
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.01));//1% resolution for not unfolding
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_HInormSyst[c]->GetBinContent(i)));//add in PbPb normalization uncert
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),s.h_HInormSyst[c]->GetBinContent(i)));//add in PbPb normalization uncert
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_normSyst->GetBinContent(i)));//add in pp normalization uncert
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),s.h_normSyst->GetBinContent(i)));//add in pp normalization uncert
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),TMath::Max(hyperonPbPb->GetBinContent(i)-1,0.015)));//PbPb hyperon study
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),TMath::Max(hyperonPbPb->GetBinContent(i)-1,0.015)));//PbPb hyperon study
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),TMath::Max(hyperonpp->GetBinContent(i)-1,0.015)));//pp hyperon study
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),TMath::Max(hyperonpp->GetBinContent(i)-1,0.015)));//pp hyperon study
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.03));//pp uncertainty for pileup
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.03));//pp uncertainty for pileup
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.04));//tight selection data/MC
      s.PbPb_totSyst[c]->SetBinContent(i,Quad(s.PbPb_totSyst[c]->GetBinContent(i),0.04));//tight selection data/MC
      
      if(c==0)s.pp_totSyst->SetBinContent(i,Quad(s.pp_totSyst->GetBinContent(i),0.12));//pplumi uncertainty for spectrum
    }
    s.RAA_totSyst[c]->Write();
    s.PbPb_totSyst[c]->Write();
    if(c==0) s.pp_totSyst->Write();

    //plotting
    canv->Clear();
    h[c]->GetXaxis()->SetRangeUser(0.7,350);
    h[c]->GetXaxis()->SetLabelOffset(-0.005);
    h[c]->GetYaxis()->SetRangeUser(0,1.45);
    h[c]->Draw();
 
    float TAAUncert = s.TAAuncert[c]/100.0;
    float lumiUncert = 0.12;//12% for pp lumi
    bLumi->SetFillColor(kGray);
    bTAA->SetFillColor(kBlue-9);
    bLumi->SetLineWidth(0);
    bTAA->SetLineWidth(0);
    bTAA->DrawBox(0.9,1-TAAUncert,TMath::Power(10,TMath::Log10(0.9)+(TMath::Log10(1.5)-TMath::Log10(0.9))/2.0),1+TAAUncert);
    bLumi->DrawBox(TMath::Power(10,TMath::Log10(0.9)+(TMath::Log10(1.5)-TMath::Log10(0.9))/2.0),1-lumiUncert,1.5,1+lumiUncert);
  
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
    //extraText  = "Unpublished";
    CMS_lumi( canv, iPeriod, 11 );
 
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();    
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.png",5*s.lowCentBin[c],5*s.highCentBin[c]));
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.pdf",5*s.lowCentBin[c],5*s.highCentBin[c]));
    canv->SaveAs(Form("plots/prettyPlots/RAA_%d_%d.C",5*s.lowCentBin[c],5*s.highCentBin[c]));

    if(c==0 || c==1 || c==23 || c==24 || c==25 || c==30){
      TCanvas * canv_276 = (TCanvas*)canv->Clone("canv_276");
      get276RAA(canv_276,s,c);
      //get276RAA(canv,s,c);
    }
    if(c==0 || c==1  || c==24 || c==31){
      TCanvas * canv_th = (TCanvas*)canv->Clone("canv_th");
      gettheoryRAA(canv_th,s,c,"");
      //gettheoryRAA(canv,s,c,"");
    }
    delete line1;
  }

  TCanvas * canv2 = new TCanvas("canv2","canv2",700,800);
  canv2->SetBorderSize(0);
  TPad * pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0,0);
  TPad * pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3,0);
  canv2->SetLineWidth(0);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.15);
  pad1->SetTopMargin(0.08);
  pad1->SetBorderSize(0);
  pad1->Draw();
  pad2->SetTopMargin(0);
  pad2->SetLeftMargin(0.15);
  pad2->SetBottomMargin(0.3);
  pad2->SetBorderSize(0);
  pad2->Draw();
  pad1->cd();
  pad1->SetLogx();
  pad1->SetLogy();
  pbpbSpec[0]->GetXaxis()->SetRangeUser(0.7,390);
  pbpbSpec[0]->GetYaxis()->SetTitle("#frac{1}{N_{evt}}E#frac{d^{3}N}{dp^{3}} (GeV^{-2})");
  pbpbSpec[0]->GetYaxis()->SetTitleOffset(1.2);
  pbpbSpec[0]->GetYaxis()->SetTitleSize(0.05);
  pbpbSpec[0]->GetYaxis()->SetLabelSize(0.04);
  pbpbSpec[0]->SetMarkerStyle(24);
  pbpbSpec[0]->Scale(10);
  pbpbSpec[0]->GetYaxis()->SetRangeUser(1e-17,1e4);
  pbpbSpec[0]->Draw();
  pbpbSpec[1]->SetMarkerColor(kBlue);
  pbpbSpec[1]->SetLineColor(kBlue);
  pbpbSpec[1]->SetMarkerStyle(25);
  pbpbSpec[1]->Scale(3);
  pbpbSpec[1]->Draw("same");
  pbpbSpec[23]->SetMarkerColor(kRed);
  pbpbSpec[23]->SetLineColor(kRed);
  pbpbSpec[23]->SetMarkerStyle(28);
  pbpbSpec[23]->Draw("same");
  pbpbSpec[24]->SetMarkerStyle(20);
  pbpbSpec[25]->SetMarkerColor(kBlue);
  pbpbSpec[25]->SetLineColor(kBlue);
  pbpbSpec[24]->Draw("same");
  pbpbSpec[25]->SetMarkerStyle(21);
  pbpbSpec[25]->Draw("same");
  pbpbSpec[30]->SetMarkerColor(kRed);
  pbpbSpec[30]->SetLineColor(kRed);
  pbpbSpec[30]->SetMarkerStyle(34);
  pbpbSpec[30]->Draw("same");
  ppSpec->SetMarkerStyle(5);
  ppSpec->Scale(1/70.0);//scaled by inelastic xsection of 70 mb
  ppSpec->Print("All");
  ppSpec->Draw("same");
  TLegend * specLeg = new TLegend(0.25,0.1,0.45,0.5);
  specLeg->AddEntry((TObject*)0,"|#eta|<1",""); 
  specLeg->AddEntry(pbpbSpec[0],Form("0-5%s (x10)","%"),"p");  
  specLeg->AddEntry(pbpbSpec[1],Form("5-10%s (x3)","%"),"p");  
  specLeg->AddEntry(pbpbSpec[23],Form("10-30%s","%"),"p");  
  specLeg->AddEntry(pbpbSpec[24],Form("30-50%s","%"),"p");  
  specLeg->AddEntry(pbpbSpec[25],Form("50-70%s","%"),"p");  
  specLeg->AddEntry(pbpbSpec[30],Form("70-90%s","%"),"p");  
  specLeg->AddEntry(ppSpec,"pp","p"); 
  specLeg->Draw("same"); 
 
  pad2->cd();
  pad2->SetLogx();
  s.PbPb_totSyst[0]->GetYaxis()->SetTitleOffset(0.8);
  s.PbPb_totSyst[0]->GetYaxis()->SetTitleFont(42);
  s.PbPb_totSyst[0]->GetYaxis()->SetTitleSize(0.08);
  s.PbPb_totSyst[0]->GetYaxis()->SetLabelSize(0.08);
  s.PbPb_totSyst[0]->GetXaxis()->SetTitleFont(42);
  s.PbPb_totSyst[0]->GetYaxis()->SetTitle(Form("Syst. Uncert. (%s)","%"));
  s.PbPb_totSyst[0]->GetXaxis()->SetRangeUser(0.7,390);
  s.PbPb_totSyst[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
  s.PbPb_totSyst[0]->GetXaxis()->SetTitleSize(0.08);
  s.PbPb_totSyst[0]->GetXaxis()->SetLabelSize(0.08);
  s.PbPb_totSyst[0]->SetFillColor(kOrange);
  s.PbPb_totSyst[0]->Scale(100);
  s.PbPb_totSyst[0]->GetYaxis()->SetRangeUser(0.0,20);
  s.PbPb_totSyst[0]->Draw();
  s.PbPb_totSyst[30]->SetFillColor(kRed);
  s.PbPb_totSyst[30]->SetFillStyle(3004);
  s.PbPb_totSyst[30]->Scale(100);
  s.PbPb_totSyst[30]->Draw("same");
  s.pp_totSyst->SetFillColor(kBlack);
  s.pp_totSyst->SetFillStyle(3005);
  s.pp_totSyst->Scale(100);
  s.pp_totSyst->Draw("same");
  TLegend * systLeg = new TLegend(0.6,0.6,0.9,0.98);
  systLeg->AddEntry(s.PbPb_totSyst[0],Form("0-5%s","%"),"f");
  systLeg->AddEntry(s.PbPb_totSyst[30],Form("70-90%s","%"),"f");
  systLeg->AddEntry(s.pp_totSyst,"pp","f");
  systLeg->Draw("same");

  CMS_lumi( canv2, 0, 33 );
  canv2->Update();
  canv2->RedrawAxis();
  //canv2->GetFrame()->Draw();    
  canv2->SaveAs("plots/prettyPlots/Spectra_perEventYield.png");
  canv2->SaveAs("plots/prettyPlots/Spectra_perEventYield.pdf");
  canv2->SaveAs("plots/prettyPlots/Spectra_perEventYield.C");

  s.PbPb_totSyst[0]->Print("All");
  s.pp_totSyst->Print("All");
  s.RAA_totSyst[30]->Print("All");

  inputPlots->Close();
  
  return;
}


void get276RAA(TCanvas * c276, Settings s, int centralityBin){
  float raaaxis[28] = {1,1.1,1.2,1.4,1.6,1.8,2,2.2,2.4,3.2,4,4.8,5.6,6.4,7.2,9.6,12,14.4,19.2,24,28.8,35.2,41.6,48,60.8,73.6,86.4,103.6};
    int tempCentralityBin = centralityBin;
    if(centralityBin==23) tempCentralityBin=2;
    if(centralityBin==24) tempCentralityBin=3;
    if(centralityBin==25) tempCentralityBin=4;
    if(centralityBin==30) tempCentralityBin=5;
    float raaval[6][27] = {{3.591e-01, 3.701e-01, 3.887e-01, 4.059e-01, 4.181e-01, 4.204e-01, 4.202e-01, 4.104e-01, 3.622e-01, 2.640e-01, 1.882e-01, 1.487e-01, 1.3604e-01, 1.380e-01, 1.519e-01, 1.820e-01, 2.226e-01, 2.768e-01, 3.678e-01, 4.003e-01, 4.395e-01, 5.270e-01, 5.161e-01, 5.303e-01, 6.400e-01, 5.331e-01, 5.202e-01},
    {3.683e-01, 3.796e-01, 4.004e-01, 4.182e-01, 4.305e-01, 4.362e-01, 4.350e-01, 4.292e-01, 3.844e-01, 2.830e-01, 2.109e-01 , 1.717e-01, 1.616e-01, 1.628e-01 , 1.791e-01, 2.103e-01, 2.458e-01, 2.958e-01, 4.071e-01, 3.688e-01, 4.087e-01, 4.884e-01, 7.80e-01, 5.579e-01, 6.219e-01, 6.075e-01, 5.558e-01},
    { 3.999e-01 , 4.134e-01, 4.344e-01, 4.539e-01 , 4.673e-01, 4.731e-01,  4.748e-01, 4.693e-01, 4.290e-01 , 3.348e-01, 2.647e-01, 2.252e-01,  2.148e-01, 2.153e-01 , 2.358e-01, 2.682e-01, 3.176e-01, 3.950e-01, 4.970e-01, 4.893e-01, 4.702e-01, 6.024e-01, 8.08e-01, 7.215e-01, 7.242e-01, 5.498e-01 , 6.87e-01},
    {4.708e-01, 4.840e-01, 5.070e-01, 5.257e-01, 5.391e-01, 5.440e-01, 5.473e-01, 5.421e-01, 5.169e-01, 4.355e-01, 3.766e-01, 3.447e-01, 3.413e-01, 3.468e-01, 3.613e-01, 4.165e-01, 4.619e-01, 5.451e-01, 6.290e-01, 6.017e-01, 7.022e-01, 7.52e-01, 8.72e-01, 1.118e+00, 8.99e-01, 1.66e+00, 7.50e-01},
    {5.568e-01,  5.675e-01, 5.843e-01, 6.004e-01, 6.076e-01, 6.139e-01, 6.167e-01, 6.176e-01, 6.088e-01, 5.544e-01, 5.236e-01, 5.099e-01, 5.152e-01, 5.115e-01, 5.412e-01, 5.832e-01, 6.419e-01, 6.805e-01, 7.546e-01,  5.90e-01, 9.08e-01, 6.36e-01, 2.110e+00, 7.416e-01, 9.84e-01, 6.98e-01, 8.60e-01 },
    {6.053e-01, 6.102e-01, 6.154e-01, 6.250e-01, 6.297e-01, 6.341e-01, 6.406e-01, 6.456e-01, 6.492e-01, 6.390e-01 , 6.256e-01, 6.057e-01, 6.107e-01, 6.185e-01, 6.418e-01, 6.522e-01, 7.560e-01, 7.219e-01, 6.55e-01, 7.05e-01, 4.37e-01, 1.026e+00, 4.690e-01, 7.209e-01, 7.81e-01, 7.51e-01, 6.36e-01}};
    float raavalstat[6][27] = {{3e-04, 3e-04, 3e-04, 4e-04, 4e-04, 5e-04, 6e-04, 7e-04, 5e-04, 7e-04, 8e-04, 1.0e-03, 1.36e-03, 1.9e-03, 2.0e-03, 4.3e-03, 8.5e-03, 1.23e-02, 2.70e-02, 4.23e-02, 5.59e-02, 7.33e-02, 8.78e-02, 2.90e-02, 5.41e-02, 5.98e-02, 9.06e-02},
    {3e-04, 4e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 9e-04, 6e-04, 8e-04, 1.0e-03, 1.2e-03, 1.7e-03, 2.4e-03, 2.4e-03, 5.1e-03, 9.6e-03, 1.37e-02, 3.12e-02 ,  4.19e-02, 5.62e-02, 7.32e-02, 1.68e-01, 3.45e-02, 6.90e-02, 7.30e-02, 9.76e-02},
    {3e-04, 3e-04, 3e-04, 3e-04, 4e-04, 5e-04, 6e-04, 7e-04, 5e-04, 7e-04, 9e-04, 1.2e-03, 1.6e-03, 2.3e-03, 2.4e-03, 5.0e-03, 9.8e-03, 1.42e-02, 3.04e-02, 4.28e-02, 4.97e-02,  6.11e-02, 1.10e-01,  8.08e-02,  5.43e-02 , 5.70e-02, 1.09e-01},
    {4e-04, 5e-04, 4e-04, 5e-04,  7e-04, 8e-04, 1.0e-03, 1.1e-03, 8e-04, 1.1e-03, 1.6e-03, 2.2e-03 , 3.1e-03, 4.3e-03, 4.1e-03, 8.8e-03, 1.62e-02, 2.25e-02, 4.54e-02, 6.47e-02, 9.08e-02, 1.10e-01, 1.75e-01, 2.16e-01, 7.2e-02, 1.01e+00, 1.30e-01},
    {8e-04, 9e-04, 8e-04 , 1.0e-03, 1.2e-03, 1.5e-03, 1.7e-03 , 2.0e-03, 1.4e-03, 2.0e-03, 2.9e-03, 4.3e-03, 6.1e-03, 8.3e-03, 8.1e-03, 1.65e-02, 3.02e-02, 3.96e-02, 7.98e-02, 1.03e-01, 1.84e-01, 1.79e-01, 7.38e-01, 4.25e-02, 9.7e-02, 1.12e-01, 1.93e-01},
    {1.7e-03, 1.9e-03, 1.5e-03, 1.8e-03, 2.1e-03, 2.5e-03, 2.9e-03, 3.3e-03, 2.3e-03,  3.8e-03,  6.0e-03, 8.9e-03, 1.30e-02, 1.83e-02 , 1.81e-02, 3.63e-02,  6.86e-02, 8.55e-02, 1.56e-01, 2.66e-01, 2.50e-01, 7.14e-01, 5.06e-02, 7.99e-02, 1.49e-01, 2.08e-01, 2.77e-01}};
    float raavalsyst[6][27] = {{2.62e-02, 2.70e-02, 2.84e-02, 2.96e-02, 3.05e-02, 3.07e-02, 3.07e-02, 3.00e-02, 2.65e-02, 1.93e-02, 1.38e-02, 1.09e-02, 1.000e-02, 1.02e-02, 1.12e-02, 1.37e-02, 1.70e-02, 2.22e-02, 3.12e-02, 3.53e-02, 4.06e-02, 5.35e-02, 5.72e-02, 6.10e-02, 7.73e-02, 6.47e-02, 6.36e-02},
    {2.69e-02, 2.77e-02, 2.92e-02 , 3.05e-02,  3.14e-02, 3.19e-02, 3.18e-02,  3.14e-02  , 2.81e-02   , 2.07e-02, 1.55e-02, 1.26e-02, 1.19e-02, 1.20e-02, 1.33e-02, 1.58e-02, 1.88e-02, 2.37e-02,  3.45e-02,  3.26e-02,  3.77e-02, 4.96e-02, 8.6e-02, 6.42e-02, 7.51e-02 ,7.38e-02,  6.79e-02 },
    {2.84e-02,  2.93e-02,  3.08e-02, 3.22e-02, 3.32e-02, 3.36e-02, 3.37e-02, 3.33e-02, 3.05e-02,  2.38e-02 , 1.88e-02, 1.60e-02, 1.53e-02, 1.54e-02, 1.69e-02, 1.93e-02, 2.30e-02, 2.93e-02, 3.80e-02, 3.87e-02, 3.86e-02, 5.54e-02, 8.2e-02, 7.68e-02, 8.15e-02,  6.23e-02, 7.8e-02},
    {3.34e-02,  3.44e-02, 3.60e-02, 3.73e-02, 3.83e-02, 3.86e-02, 3.89e-02, 3.85e-02, 3.67e-02, 3.10e-02, 2.68e-02, 2.45e-02, 2.43e-02, 2.47e-02, 2.58e-02, 3.00e-02, 3.34e-02, 4.04e-02, 4.81e-02 , 4.76e-02, 5.76e-02 , 6.9e-02, 8.9e-02, 1.19e-01, 1.01e-01 , 1.9e-01, 8.6e-02},
    {3.95e-02, 4.03e-02, 4.15e-02, 4.26e-02, 4.31e-02 , 4.36e-02 , 4.38e-02, 4.39e-02,  4.33e-02, 3.94e-02, 3.72e-02, 3.63e-02, 3.67e-02, 3.65e-02, 3.87e-02 , 4.20e-02, 4.65e-02, 5.04e-02, 5.77e-02, 4.7e-02,  7.5e-02, 5.8e-02, 2.15e-01, 7.89e-02, 1.11e-01, 7.9e-02, 9.8e-02},
    {4.30e-02, 4.33e-02, 4.37e-02, 4.44e-02, 4.47e-02, 4.50e-02, 4.55e-02, 4.59e-02, 4.61e-02, 4.54e-02, 4.45e-02, 4.31e-02, 4.35e-02, 4.41e-02, 4.59e-02, 4.69e-02, 5.47e-02 , 5.35e-02, 5.0e-02, 5.6e-02, 3.6e-02, 9.4e-02, 4.78e-02, 7.67e-02, 8.8e-02, 8.5e-02,  7.3e-02 }};
//ATLAS Data below
double p8800_d40x1y1_xval[] = { 0.5365, 0.615, 0.7050000000000001, 0.808, 0.9259999999999999, 1.0594999999999999, 1.21, 1.385, 1.5899999999999999, 
    1.825, 2.095, 2.4050000000000002, 2.755, 3.1550000000000002, 3.62, 4.15, 4.755, 5.455, 6.255, 
    7.17, 8.219999999999999, 9.39, 10.75, 12.35, 14.149999999999999, 16.2, 18.6, 21.35, 24.450000000000003, 
    28.05, 33.85, 42.6, 53.65, 67.55, 85.05, 106.9, 134.5 };
  double p8800_d40x1y1_xerrminus[] = { 0.03649999999999998, 0.04200000000000004, 0.04800000000000004, 0.05500000000000005, 0.06299999999999994, 0.0704999999999999, 0.08000000000000007, 0.09499999999999997, 0.10999999999999988, 
    0.125, 0.14500000000000024, 0.16500000000000004, 0.18500000000000005, 0.2150000000000003, 0.25, 0.28000000000000025, 0.3250000000000002, 0.375, 0.4249999999999998, 
    0.4900000000000002, 0.5599999999999987, 0.6100000000000012, 0.75, 0.8499999999999996, 0.9499999999999993, 1.0999999999999996, 1.3000000000000007, 1.4500000000000028, 1.6500000000000021, 
    1.9499999999999993, 3.8500000000000014, 4.899999999999999, 6.149999999999999, 7.75, 9.75, 12.100000000000009, 15.5 };
  double p8800_d40x1y1_xerrplus[] = { 0.03649999999999998, 0.04200000000000004, 0.04799999999999993, 0.05499999999999994, 0.06300000000000006, 0.07050000000000001, 0.08000000000000007, 0.09499999999999997, 0.1100000000000001, 
    0.125, 0.14500000000000002, 0.1649999999999996, 0.18500000000000005, 0.21499999999999986, 0.25, 0.27999999999999936, 0.3250000000000002, 0.375, 0.4249999999999998, 
    0.4900000000000002, 0.5600000000000005, 0.6099999999999994, 0.75, 0.8499999999999996, 0.9500000000000011, 1.1000000000000014, 1.2999999999999972, 1.4499999999999993, 1.6499999999999986, 
    1.9499999999999993, 3.8500000000000014, 4.899999999999999, 6.149999999999999, 7.75, 9.75, 12.099999999999994, 15.5 };
  double p8800_d40x1y1_yval[] = { 0.271, 0.289, 0.311, 0.336, 0.361, 0.388, 0.413, 0.434, 0.452, 
    0.452, 0.451, 0.427, 0.388, 0.336, 0.277, 0.221, 0.18, 0.154, 0.142, 
    0.141, 0.15, 0.164, 0.179, 0.201, 0.225, 0.245, 0.276, 0.322, 0.348, 
    0.351, 0.43, 0.51, 0.545, 0.562, 0.609, 0.614, 0.645 };
  double p8800_d40x1y1_yerrminus[] = { 0.013550227909083301, 0.013872270934854176, 0.014306324856094243, 0.01545642218553828, 0.01696752574610983, 0.01823669384637468, 0.019411929667099048, 0.019965226955774882, 0.02124573099330781, 
    0.02079447651772941, 0.02074962533106562, 0.019647615175639, 0.017467284258292703, 0.015130785486550263, 0.012479887553980606, 0.00996562982505371, 0.008309705409940835, 0.006824980471766933, 0.006475822777068562, 
    0.0063590999999999995, 0.006738323233564861, 0.007493920736170086, 0.00786989536143906, 0.00896423945463306, 0.009902556488099425, 0.011100960544025009, 0.012590562179664576, 0.016180299750004633, 0.017744587907302888, 
    0.017897558492710672, 0.022757524030526697, 0.029549208111216793, 0.039004641390480696, 0.045396939324143876, 0.057652714784301354, 0.0751115513087035, 0.11254279685968356 };
  double p8800_d40x1y1_yerrplus[] = { 0.013550227909083301, 0.013872270934854176, 0.014306324856094243, 0.01545642218553828, 0.01696752574610983, 0.01823669384637468, 0.019411929667099048, 0.019965226955774882, 0.02124573099330781, 
    0.02079447651772941, 0.02074962533106562, 0.019647615175639, 0.017467284258292703, 0.015130785486550263, 0.012479887553980606, 0.00996562982505371, 0.008309705409940835, 0.006824980471766933, 0.006475822777068562, 
    0.0063590999999999995, 0.006738323233564861, 0.007493920736170086, 0.00786989536143906, 0.00896423945463306, 0.009902556488099425, 0.011100960544025009, 0.012590562179664576, 0.016180299750004633, 0.017744587907302888, 
    0.017897558492710672, 0.022757524030526697, 0.029549208111216793, 0.039004641390480696, 0.045396939324143876, 0.057652714784301354, 0.0751115513087035, 0.11254279685968356 };
  double p8800_d40x1y1_ystatminus[] = { 7.859E-5, 8.669999999999999E-5, 9.641E-5, 1.1424000000000002E-4, 1.3356999999999998E-4, 1.5908E-4, 1.8998E-4, 2.2134E-4, 2.712E-4, 
    3.2091999999999993E-4, 3.8785999999999996E-4, 4.696999999999999E-4, 5.044E-4, 5.712000000000001E-4, 6.094000000000001E-4, 6.409E-4, 7.02E-4, 8.162E-4, 0.0010508, 
    0.0013958999999999998, 0.0019500000000000001, 0.002952, 0.002327, 0.0030150000000000003, 0.0036000000000000008, 0.004164999999999999, 0.005520000000000001, 0.008372000000000001, 0.007656, 
    0.007722, 0.0086, 0.01071, 0.010355, 0.016860000000000003, 0.029841000000000003, 0.052804, 0.094815 };
  double p8800_d40x1y1_ystatplus[] = { 7.859E-5, 8.669999999999999E-5, 9.641E-5, 1.1424000000000002E-4, 1.3356999999999998E-4, 1.5908E-4, 1.8998E-4, 2.2134E-4, 2.712E-4, 
    3.2091999999999993E-4, 3.8785999999999996E-4, 4.696999999999999E-4, 5.044E-4, 5.712000000000001E-4, 6.094000000000001E-4, 6.409E-4, 7.02E-4, 8.162E-4, 0.0010508, 
    0.0013958999999999998, 0.0019500000000000001, 0.002952, 0.002327, 0.0030150000000000003, 0.0036000000000000008, 0.004164999999999999, 0.005520000000000001, 0.008372000000000001, 0.007656, 
    0.007722, 0.0086, 0.01071, 0.010355, 0.016860000000000003, 0.029841000000000003, 0.052804, 0.094815 };
  int p8800_d40x1y1_numpoints = 37;
  //TGraphAsymmErrors p8800_d40x1y1 = TGraphAsymmErrors(p8800_d40x1y1_numpoints, p8800_d40x1y1_xval, p8800_d40x1y1_yval, p8800_d40x1y1_xerrminus, p8800_d40x1y1_xerrplus, p8800_d40x1y1_yerrminus, p8800_d40x1y1_yerrplus);
  TGraphAsymmErrors p8800_d40x1y1 = TGraphAsymmErrors(p8800_d40x1y1_numpoints, p8800_d40x1y1_xval, p8800_d40x1y1_yval,0,0, p8800_d40x1y1_yerrminus, p8800_d40x1y1_yerrplus);
  p8800_d40x1y1.SetName("ATLAS_0_5");
  p8800_d40x1y1.SetTitle("ATLAS_0_5");
  p8800_d40x1y1.SetMarkerColor(kBlue);
  p8800_d40x1y1.SetMarkerColor(kBlue);
  p8800_d40x1y1.SetLineColor(kBlue);
//end ATLAS data

  TH1D * p;
  p = new TH1D("raa276",";p_{T};R_{AA}",27,raaaxis);
  for(int i = 1; i<p->GetSize()-1; i++){
    p->SetBinContent(i,raaval[tempCentralityBin][i-1]);
    p->SetBinError(i,raavalstat[tempCentralityBin][i-1]); 
  }
  p->SetMarkerColor(kRed);
  p->SetLineColor(kRed);
  p->Draw("same");
  TBox *bp[27];
  for(int i = 0; i<27; i++) bp[i] = new TBox(0.1,0.1,0.2,0.2);
  for(int i = 1; i<p->GetSize()-1; i++){
    bp[i-1]->SetFillStyle(0);
    bp[i-1]->SetLineColor(kRed);
    bp[i-1]->SetLineWidth(1);
    bp[i-1]->SetX1(p->GetXaxis()->GetBinLowEdge(i));
    bp[i-1]->SetX2(p->GetXaxis()->GetBinUpEdge(i));
    bp[i-1]->SetY1((p->GetBinContent(i))*(1-raavalsyst[tempCentralityBin][i-1]));
    bp[i-1]->SetY2((p->GetBinContent(i)*(1+raavalsyst[tempCentralityBin][i-1])));
    bp[i-1]->Draw("same");
  }
  TLegend * legRaa276 = new TLegend(0.53,0.725,0.9,0.91);
  TH1D * dummy = new TH1D("dummy","dummy",10,0,10);
  dummy->SetMarkerColor(kBlack); dummy->SetMarkerStyle(8); dummy->SetFillColor(kOrange);
  legRaa276->AddEntry(dummy,"CMS 5.02 TeV R_{AA}","pf");
  legRaa276->AddEntry(p,"CMS 2.76 TeV R_{AA}","p");

  gStyle->SetErrorX(0.);
  if(centralityBin==0) legRaa276->AddEntry(&p8800_d40x1y1,"ATLAS 2.76 TeV R_{AA}","p");
  legRaa276->SetTextFont(62);
  legRaa276->Draw("same");
  if(centralityBin==0) p8800_d40x1y1.Draw("P same");
  c276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.C",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin]));
  c276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.png",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin]));
  c276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.pdf",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin]));
 
  if(centralityBin==0){
    TH1D * dummy2 = new TH1D("dummy2","dummy2",10,0,10);
    dummy2->SetFillStyle(3002);dummy2->SetFillColor(kRed);dummy2->SetLineWidth(0);
    legRaa276->AddEntry(dummy2,"Y. Chien et al. 0-10%","f");
    legRaa276->Draw("same");
    gettheoryRAA(c276,s,centralityBin,"With276");
    delete dummy2;
  }
  delete legRaa276;
  delete dummy;
  for(int i = 0; i<27; i++) delete bp[i];
  return;
}

void gettheoryRAA(TCanvas * c_th, Settings s, int centralityBin, std::string saveString = ""){
  //Vitev**********************************************************************************************************
  float temp_x;
  float temp_y;
  vector<float> x;
  vector<float> y_d;
  vector<float> y_u;
  int theoryCent_Low, theoryCent_High;
  if(centralityBin==0 || centralityBin==1 || centralityBin==31){theoryCent_Low=0; theoryCent_High=10;}
  if(centralityBin==24){theoryCent_Low=30; theoryCent_High=50;}
  ifstream input_file_d(Form("theoryPredictions/IvanVitev/R-%d%d.5100GeVch.dn",theoryCent_Low,theoryCent_High));
  ifstream input_file_u(Form("theoryPredictions/IvanVitev/R-%d%d.5100GeVch.up",theoryCent_Low,theoryCent_High));
  //get datai
  std::cout << "reading theory prediction data" << std::endl;
  while(!input_file_d.eof()){ 
    input_file_d>>temp_x;
    input_file_d>>temp_y;
    //std::cout << temp_x << " " << temp_y << std::endl;
    x.push_back(temp_x);
    y_d.push_back(temp_y);
  }
  while(!input_file_u.eof()){ 
    input_file_u>>temp_x;
    input_file_u>>temp_y;
    //std::cout << temp_y << std::endl;
    y_u.push_back(temp_y);
  }
  std::cout << "done reading " << x.size() << " Points"  << std::endl;
  
  //put data in histograms
  const int graphPts = 1952;
  TGraph * vitev = new TGraph(2*graphPts);
  for (int i=0;i<graphPts;i++) {
    //std::cout << x[i] << " " << y_d[i] << " " << y_u[i] << std::endl;
    vitev->SetPoint(i,x[i],y_d[i]);
    vitev->SetPoint(graphPts+i,x[graphPts-i-1],y_u[graphPts-i-1]);
  }
  vitev->SetFillStyle(3002);
  vitev->SetFillColor(kRed);
  vitev->SetLineWidth(0);
  vitev->Draw("same f");

  //Jiechen Xu ************************************************************************************************
  vector<float> x2;
  vector<float> y2;
  int theoryCent_Low2, theoryCent_High2;
  if(centralityBin==0 || centralityBin==1 || centralityBin==31){theoryCent_Low2=0; theoryCent_High2=10;}
  if(centralityBin==24){theoryCent_Low2=30; theoryCent_High2=50;}
  std::cout << "Reading Jiechen Xu points" << std::endl;
  ifstream input_file_JX(Form("theoryPredictions/JiechenXu/CUJET3_RAA-pT_%d-%d.dat",theoryCent_Low2,theoryCent_High2)); 
  while(!input_file_JX.eof()){
      input_file_JX>>temp_x;
      input_file_JX>>temp_y;
      x2.push_back(temp_x);
      y2.push_back(temp_y);
      for(int i = 0; i<5; i++) input_file_JX>>temp_y;//skip to next line
    }         
  std::cout << "Done reading" << std::endl;
  //put data in histograms
  const int graphPts2 = 371;
  TGraph * JiechenXu = new TGraph(graphPts2); 
  for (int i=0;i<graphPts2;i++) {
    //std::cout << x2[i] << " " << y2[i]  << std::endl;
    JiechenXu->SetPoint(i,x2[i],y2[i]);
  }
  JiechenXu->SetLineWidth(3);
  JiechenXu->SetLineColor(kBlue+1);
  JiechenXu->Draw("same C");

  TLegend * leg_th = new TLegend(0.5,0.75,0.9,0.85);
  leg_th->AddEntry(vitev,Form("Y. Chien et al. %d-%d%s",theoryCent_Low,theoryCent_High,"%"),"f");
  leg_th->AddEntry(JiechenXu,Form("J. Xu et al. %d-%d%s (h^{#pm}+#pi^{0})",theoryCent_Low2,theoryCent_High2,"%"),"l");
  if(saveString=="") leg_th->Draw("same");
  c_th->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory%s.C",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin],saveString.c_str()));
  c_th->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory%s.png",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin],saveString.c_str()));
  c_th->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory%s.pdf",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin],saveString.c_str()));
  delete leg_th;
  return;
}
#endif
