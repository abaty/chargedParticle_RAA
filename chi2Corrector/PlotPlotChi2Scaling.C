#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLegend.h>
#include <TLatex.h>

TH1D* divideHistosDiffBins(TH1D* h_Num, TH1D* h_Den);
TH1D* divideHistosDiffBinsBothErrors(TH1D* h_Num, TH1D* h_Den);
void normalizeByBinWidth(TH1D *histo);
void normalizeByBinWidth(TH2D *histo);
void setHistoOutlook(TH1D *histo, int color, int markerStyle);
TH1D* From_xT_to_pT_histo(TH1D *h_input, TH1D *h_binning, double sqrts, double exp);
void convertToYield(TH1D* h_ppRef_LogLogpt, float sigma_inel);
void convertToYieldWithErrors(TH1D* h_ppRef_LogLogpt, float sigma_inel);
void divideKK(TH1D* h1, TH1D* h2);
//------------------------

void PlotPlotChi2Scaling()
{
   gROOT->Reset();
   gROOT->ProcessLine(".x rootlogonChristof.C");
   gROOT->ForceStyle();
   gStyle->SetPalette(1);
   gStyle->SetTitleYOffset(1.25);
   gStyle->SetTitleXOffset(1.03);
//   gStyle->SetPadRightMargin(0.16);
   gStyle->SetOptFit(0);

   TFile *f_out = new TFile("PlotPlotChi2Scaling_PbPb.root","recreate");

   bool doSave = false;
//   bool doSave = true;

   TFile *f_0_5 = new TFile("PlotPlotTrackProperties_ForChi2Correction_0_5.root");
   TFile *f_5_10 = new TFile("PlotPlotTrackProperties_ForChi2Correction_5_10.root");
   TFile *f_10_30 = new TFile("PlotPlotTrackProperties_ForChi2Correction_10_30.root");
   TFile *f_30_50 = new TFile("PlotPlotTrackProperties_ForChi2Correction_30_50.root");
   TFile *f_50_70 = new TFile("PlotPlotTrackProperties_ForChi2Correction_50_70.root");
   TFile *f_70_90 = new TFile("PlotPlotTrackProperties_ForChi2Correction_70_90.root");
   TFile *f_30_100 = new TFile("PlotPlotTrackProperties_ForChi2Correction_30_100.root");

   TH1D *hScaling_0_5 = (TH1D*)f_0_5->Get("hScaling");
   TH1D *hScaling_5_10 = (TH1D*)f_5_10->Get("hScaling");
   TH1D *hScaling_10_30 = (TH1D*)f_10_30->Get("hScaling");
   TH1D *hScaling_30_50 = (TH1D*)f_30_50->Get("hScaling");
   TH1D *hScaling_50_70 = (TH1D*)f_50_70->Get("hScaling");
   TH1D *hScaling_70_90 = (TH1D*)f_70_90->Get("hScaling");
   TH1D *hScaling_30_100 = (TH1D*)f_30_100->Get("hScaling");

   hScaling_0_5->SetName("hScaling_0_5");
   hScaling_5_10->SetName("hScaling_5_10");
   hScaling_10_30->SetName("hScaling_10_30");
   hScaling_30_50->SetName("hScaling_30_50");
   hScaling_50_70->SetName("hScaling_50_70");
   hScaling_70_90->SetName("hScaling_70_90");

   TCanvas *c0 = new TCanvas("c0","c0");
   c0->cd();
   c0->SetLogx();
   
   hScaling_0_5->GetXaxis()->CenterTitle();
   hScaling_0_5->GetYaxis()->CenterTitle();
   hScaling_0_5->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hScaling_0_5->GetYaxis()->SetTitle("#chi^{2} scaling");
   hScaling_0_5->SetMarkerStyle(20);
   hScaling_0_5->SetMinimum(0.);
   hScaling_0_5->SetMaximum(2.);
   hScaling_0_5->Draw("p");

   hScaling_5_10->SetMarkerStyle(20);
   hScaling_5_10->SetMarkerColor(4);
   hScaling_5_10->SetLineColor(4);
   hScaling_5_10->Draw("psame");

   hScaling_10_30->SetMarkerStyle(20);
   hScaling_10_30->SetMarkerColor(6);
   hScaling_10_30->SetLineColor(6);
   hScaling_10_30->Draw("psame");

   hScaling_30_50->SetMarkerStyle(20);
   hScaling_30_50->SetMarkerColor(8);
   hScaling_30_50->SetLineColor(8);
   hScaling_30_50->Draw("psame");

   hScaling_50_70->SetMarkerStyle(20);
   hScaling_50_70->SetMarkerColor(46);
   hScaling_50_70->SetLineColor(46);
   hScaling_50_70->Draw("psame");

   hScaling_70_90->SetMarkerStyle(20);
   hScaling_70_90->SetMarkerColor(42);
   hScaling_70_90->SetLineColor(42);
   hScaling_70_90->Draw("psame");

   hScaling_30_100->SetMarkerStyle(20);
   hScaling_30_100->SetMarkerColor(2);
   hScaling_30_100->SetLineColor(2);
   hScaling_30_100->Draw("psame");

   TF1* fit_0_5 = new TF1("fit_0_5","[0]",10.,40.);
   TF1* fit_5_10 = new TF1("fit_5_10","[0]",20.,60.);
   TF1* fit_10_30 = new TF1("fit_10_30","[0]",20.,60.);
   TF1* fit_30_50 = new TF1("fit_30_50","[0]",20.,60.);
   TF1* fit_50_70 = new TF1("fit_50_70","[0]",20.,60.);
   TF1* fit_70_90 = new TF1("fit_70_90","[0]",20.,60.);

   hScaling_0_5->Fit("fit_0_5","R0","same");
   hScaling_5_10->Fit("fit_5_10","R0","same");
   hScaling_10_30->Fit("fit_10_30","R0","same");
   hScaling_30_50->Fit("fit_30_50","R0","same");
   hScaling_50_70->Fit("fit_50_70","R0","same");
   hScaling_70_90->Fit("fit_70_90","R0","same");

   //to plot
   TF1* fit_0_5_lemma = new TF1("fit_0_5_lemma","[0]",30.,400.);
   TF1* fit_5_10_lemma = new TF1("fit_5_10_lemma","[0]",30.,400.);
   TF1* fit_10_30_lemma = new TF1("fit_10_30_lemma","[0]",30.,400.);
   TF1* fit_30_50_lemma = new TF1("fit_30_50_lemma","[0]",30.,400.);
   TF1* fit_50_70_lemma = new TF1("fit_50_70_lemma","[0]",30.,400.);
   TF1* fit_70_90_lemma = new TF1("fit_70_90_lemma","[0]",30.,400.);

   fit_0_5_lemma->SetParameter(0,fit_0_5->GetParameter(0));
   fit_5_10_lemma->SetParameter(0,fit_5_10->GetParameter(0));
   fit_10_30_lemma->SetParameter(0,fit_10_30->GetParameter(0));
   fit_30_50_lemma->SetParameter(0,fit_30_50->GetParameter(0));
   fit_50_70_lemma->SetParameter(0,fit_50_70->GetParameter(0));
   fit_70_90_lemma->SetParameter(0,fit_70_90->GetParameter(0));

   fit_0_5_lemma->SetLineColor(1);
   fit_5_10_lemma->SetLineColor(4);
   fit_10_30_lemma->SetLineColor(6);
   fit_30_50_lemma->SetLineColor(8);
   fit_50_70_lemma->SetLineColor(46);
   fit_70_90_lemma->SetLineColor(42);

   fit_0_5_lemma->Draw("same");
   fit_5_10_lemma->Draw("same");
   fit_10_30_lemma->Draw("same");
   fit_30_50_lemma->Draw("same");
   fit_50_70_lemma->Draw("same");
   fit_70_90_lemma->Draw("same");

   std::cerr<<" factor, 0-5%  : " << fit_0_5->GetParameter(0) << std::endl;
   std::cerr<<" factor, 5-10% : " << fit_5_10->GetParameter(0) << std::endl;
   std::cerr<<" factor, 10-30%: " << fit_10_30->GetParameter(0) << std::endl;
   std::cerr<<" factor, 30-50%: " << fit_30_50->GetParameter(0) << std::endl;
   std::cerr<<" factor, 50-70%: " << fit_50_70->GetParameter(0) << std::endl;
   std::cerr<<" factor, 70-90%: " << fit_70_90->GetParameter(0) << std::endl;

   f_out->cd();
   hScaling_0_5->Write();
   hScaling_5_10->Write();
   hScaling_10_30->Write();
   hScaling_30_50->Write();
   hScaling_50_70->Write();
   hScaling_70_90->Write();
   fit_0_5_lemma->Write();
   fit_5_10_lemma->Write();
   fit_10_30_lemma->Write();
   fit_30_50_lemma->Write();
   fit_50_70_lemma->Write();
   fit_70_90_lemma->Write();
   f_out->Close();
/*
   gStyle->SetOptFit(0);
   TF1 *f_fit1 = new TF1("f_fit1","[0]+[1]*log10(x)",0.7,20.);
//   TF1 *f_fit2 = new TF1("f_fit2","[0]+[1]*log(x)",20.,160.);
   TF1 *f_fit2 = new TF1("f_fit2","[0]",20.,400.);
//   TF1 *f_fit = new TF1("f_fit","f_fit1(0)+f_fit2(2)",0.7,400.);

//   hScaling_0_5->Fit("pol3","","",0.7,140.);
   Double_t par[4];
   hScaling->Fit(f_fit1,"R");
   double param = f_fit1->Eval(20.);
   std::cerr<<"param: " << param << std::endl;
   f_fit2->FixParameter(0,param);
//   hScaling->Fit(f_fit2,"R+");
   f_fit2->Draw("same");   

   f_fit2->GetParameters(&par[2]);
//   f_fit->SetParameters(par);
//   hScaling->Fit(f_fit,"R+");
*/
   TLegend *leg3 = new TLegend(0.2,0.75,0.74,0.95,"","brNDC");
   leg3->AddEntry(hScaling_0_5,"PbPb, 0-5% ","pl");
   leg3->AddEntry(hScaling_5_10,"PbPb, 5-10% ","pl");
   leg3->AddEntry(hScaling_10_30,"PbPb, 10-30% ","pl");
   leg3->AddEntry(hScaling_30_50,"PbPb, 30-50% ","pl");
   leg3->AddEntry(hScaling_50_70,"PbPb, 50-70% ","pl");
   leg3->AddEntry(hScaling_70_90,"PbPb, 70-90% ","pl");
   leg3->AddEntry(hScaling_30_100,"PbPb, 30-100%","pl");
   leg3->SetFillStyle(0);
   leg3->SetFillColor(0);
   leg3->SetBorderSize(0);
   leg3->Draw();

   if(doSave) {
      c0->SaveAs("Figs/PlotPlotChi2Scaling_c0.gif");
      c0->SaveAs("Figs/PlotPlotChi2Scaling_c0.eps");
      c0->SaveAs("Figs/PlotPlotChi2Scaling_c0.C");
   }
}


TH1D* divideHistosDiffBins(TH1D* h_Num, TH1D* h_Den) {
   TH1D *h_ratio = (TH1D*)h_Num->Clone("h_ratio");
   h_ratio->Reset();
   for(int i = 1; i <= h_Num->GetNbinsX(); i++) {
      float content = h_Num->GetBinContent(i);
      float error = h_Num->GetBinError(i);
      float center = h_Num->GetBinCenter(i);
      int which_bin_in_h_Den = h_Den->FindBin(center);
      float content_h_Den = h_Den->GetBinContent(which_bin_in_h_Den);

      if(content_h_Den==0)
         continue;

      h_ratio->SetBinContent(i,content/content_h_Den);      
      h_ratio->SetBinError(i,error/content_h_Den);      
   }
   return h_ratio;
}


TH1D* divideHistosDiffBinsBothErrors(TH1D* h_Num, TH1D* h_Den) {
   TH1D *h_ratio = (TH1D*)h_Num->Clone("h_ratio");
   h_ratio->Reset();
   for(int i = 1; i <= h_Num->GetNbinsX(); i++) {
      float content = h_Num->GetBinContent(i);
      float error = h_Num->GetBinError(i);
      float center = h_Num->GetBinCenter(i);
      int which_bin_in_h_Den = h_Den->FindBin(center);
      float content_h_Den = h_Den->GetBinContent(which_bin_in_h_Den);
      float error_h_Den = h_Den->GetBinError(which_bin_in_h_Den);

      if(content_h_Den==0)
         continue;

      h_ratio->SetBinContent(i,content/content_h_Den);      
//      std::cerr<<"error/content: " << error/content << std::endl;
//      std::cerr<<"error/content: " << error/content << std::endl;
      h_ratio->SetBinError(i,(content/content_h_Den)*TMath::Sqrt((error/content)*(error/content) + (error_h_Den/content_h_Den)*(error_h_Den/content_h_Den)));     
   }
   return h_ratio;
}


void normalizeByBinWidth(TH2D *histo) {
   for(int i = 1; i <= histo->GetNbinsX(); i++) {
      float content = histo->GetBinContent(i);
      float error = histo->GetBinError(i);
      histo->SetBinContent(i,content/histo->GetBinWidth(i));
      histo->SetBinError(i,error/histo->GetBinWidth(i));
   }
}


void normalizeByBinWidth(TH1D *histo) {
   for(int i = 1; i <= histo->GetNbinsX(); i++) {
      float content = histo->GetBinContent(i);
      float error = histo->GetBinError(i);
      histo->SetBinContent(i,content/histo->GetBinWidth(i));
      histo->SetBinError(i,error/histo->GetBinWidth(i));
   }
}


void setHistoOutlook(TH1D *hist, int color, int markerStyle) {
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(1.2);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
}


TH1D* From_xT_to_pT_histo(TH1D *h_input, TH1D *h_binning, double sqrts, double exp) {

  TH1D *h_pt = (TH1D*)h_binning->Clone("h_pt");
  h_pt->Reset();

  for(int i=1; i<=h_input->GetNbinsX(); i++) {
     double xT = h_input->GetBinCenter(i);
     double content = h_input->GetBinContent(i)/TMath::Power(sqrts,exp);
     double pT = 5020.*xT/2.;
     if(pT>120) break;
     h_pt->SetBinContent(h_pt->FindBin(pT),content);
  }

  return h_pt;
}


void convertToYield(TH1D* h_ppRef_LogLogpt, float sigma_inel) {

   for(int i = 1; i <= h_ppRef_LogLogpt->GetNbinsX(); i++) {
      float content_Ed3sigmadp3 = h_ppRef_LogLogpt->GetBinContent(i);
      float bincenter = h_ppRef_LogLogpt->GetBinCenter(i);
      //multiply by 2*pi*pT
      float content_d2sigmadetadpT = content_Ed3sigmadp3*2.*TMath::Pi()*bincenter;
      //multiply by 2: reference pp is for the average of the + and - particles
      content_d2sigmadetadpT *= 2.;
      //convert d2sigmadetadpT to d2NdetadpT
      float content_d2NdetadpT = content_d2sigmadetadpT/sigma_inel;
      //convert d2NdetadpT to dNdpT in |eta|<1
      content_d2NdetadpT *= 2.;
      h_ppRef_LogLogpt->SetBinContent(i,content_d2NdetadpT);
   }
}


void convertToYieldWithErrors(TH1D* h_ppRef_LogLogpt, float sigma_inel) {

   for(int i = 1; i <= h_ppRef_LogLogpt->GetNbinsX(); i++) {
      float content_Ed3sigmadp3 = h_ppRef_LogLogpt->GetBinContent(i);
      float error_Ed3sigmadp3 = h_ppRef_LogLogpt->GetBinError(i);
      float bincenter = h_ppRef_LogLogpt->GetBinCenter(i);
      //multiply by 2*pi*pT
      float content_d2sigmadetadpT = content_Ed3sigmadp3*2.*TMath::Pi()*bincenter;
      float error_d2sigmadetadpT = error_Ed3sigmadp3*2.*TMath::Pi()*bincenter;
      //multiply by 2: reference pp is for the average of the + and - particles
      content_d2sigmadetadpT *= 2.;
      error_d2sigmadetadpT *= 2.;
      //convert d2sigmadetadpT to d2NdetadpT
      float content_d2NdetadpT = content_d2sigmadetadpT/sigma_inel;
      float error_d2NdetadpT = error_d2sigmadetadpT/sigma_inel;
      //convert d2NdetadpT to dNdpT in |eta|<1
      content_d2NdetadpT *= 2.;
      error_d2NdetadpT *= 2.;
      h_ppRef_LogLogpt->SetBinContent(i,content_d2NdetadpT);
      h_ppRef_LogLogpt->SetBinError(i,error_d2NdetadpT);
   }
}

void divideKK(TH1D* h1, TH1D* h2)
{
   h1->Divide(h2);
   for(int i = 1; i<h1->GetNbinsX(); i++) {
      double epsilon = h1->GetBinContent(i);
      double sample = h2->GetBinContent(i);
      double binomError = TMath::Sqrt(epsilon*(1-epsilon)/sample);
      h1->SetBinError(i,binomError);
   }
}
