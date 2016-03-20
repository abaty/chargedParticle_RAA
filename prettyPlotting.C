#ifndef PRETTYPLOT
#define PRETTYPLOT

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
void gettheoryRAA(TCanvas * c_th, Settings s, int centralityBin);

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
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),0.01));//for nonclosure pp 
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_HInormSyst[c]->GetBinContent(i)));//add in PbPb normalization uncert
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),s.h_normSyst->GetBinContent(i)));//add in pp normalization uncert
      
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),hyperonPbPb->GetBinContent(i)-1));//PbPb hyperon study
      s.RAA_totSyst[c]->SetBinContent(i,Quad(s.RAA_totSyst[c]->GetBinContent(i),hyperonpp->GetBinContent(i)-1));//pp hyperon study
    }
    s.RAA_totSyst[c]->Write();

    //plotting
    canv->Clear();
    h[c]->GetXaxis()->SetRangeUser(0.7,350);
    h[c]->GetYaxis()->SetRangeUser(0,1.45);
    h[c]->Draw();
  
    float lumiUncert = Quad(0.05,s.TAAuncert[c]/100.0);//10% for pp lumi, added in quad with TAA
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

    if(c==0 || c==1 || c==23 || c==24 || c==25 || c==30){
      TCanvas * canv_276 = (TCanvas*)canv->Clone("canv_276");
      get276RAA(canv_276,s,c);
    }
    if(c==0 || c==1  || c==24){
      TCanvas * canv_th = (TCanvas*)canv->Clone("canv_th");
      gettheoryRAA(canv_th,s,c);
    }
    delete line1;
  }
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
  TLegend * legRaa276 = new TLegend(0.5,0.75,0.9,0.9);
  legRaa276->AddEntry(p,"CMS 2.76 TeV R_{AA}","p");
  legRaa276->Draw("same");
  c276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.png",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin]));
  c276->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_Compare276.pdf",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin]));
  delete legRaa276;
  for(int i = 0; i<27; i++) delete bp[i];
  return;
}

void gettheoryRAA(TCanvas * c_th, Settings s, int centralityBin){
  float temp_x;
  float temp_y;
  vector<float> x;
  vector<float> y_d;
  vector<float> y_u;
  int theoryCent_Low, theoryCent_High;
  if(centralityBin==0 || centralityBin==1){theoryCent_Low=0; theoryCent_High=10;}
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
  TLegend * leg_th = new TLegend(0.5,0.75,0.9,0.85);
  leg_th->AddEntry(vitev,Form("Y. Chien et al. %d-%d%s (arXiv:1509.02936)",theoryCent_Low,theoryCent_High,"%"),"f");
  leg_th->Draw("same");
  c_th->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory.png",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin]));
  c_th->SaveAs(Form("plots/prettyPlots/RAA_%d_%d_CompareTheory.pdf",5*s.lowCentBin[centralityBin],5*s.highCentBin[centralityBin]));
  delete leg_th;
  return;
}
#endif
