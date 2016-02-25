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

void makePlotsPbPb(Settings s)
{
  for(int c = 0; c<s.nCentBins; c++){
    for(int j = 0; j<s.HInTriggers; j++){
      s.HIByTrigger[j][c]->SetLineColor(j+1);
      s.HIByTrigger[j][c]->SetLineWidth(1);
      s.HIByTrigger[j][c]->SetMarkerColor(j+1);
      s.HIByTrigger[j][c]->SetFillColor(j+1);
      s.HIUsedByTrigger[j][c]->SetLineColor(kBlack);
      s.HIUsedByTrigger[j][c]->SetLineWidth(2);
      s.HIUsedByTrigger[j][c]->SetMarkerColor(j+1);
      s.HIUsedByTrigger[j][c]->SetMarkerSize(0);
      s.HIUsedByTrigger[j][c]->SetFillColor(j+1);
      s.HIJetsByTrigger[j][c]->SetLineColor(j+1);
      s.HIJetsByTrigger[j][c]->SetLineWidth(1);
      s.HIJetsByTrigger[j][c]->SetMarkerColor(j+1);
      s.HIJetsByTrigger[j][c]->SetMarkerSize(0.8);
      s.HIJetsByTrigger[j][c]->SetFillColor(j+1);
      if(j==4){
        s.HIJetsByTrigger[j][c]->SetMarkerColor(kOrange-3);
        s.HIJetsByTrigger[j][c]->SetLineColor(kOrange-3);
      }
    }
    for(int j = 0; j<s.HInTriggers_trk; j++){
      s.HIByTrigger_trk[j][c]->SetLineColor(j+1);
      s.HIByTrigger_trk[j][c]->SetLineWidth(1);
      s.HIByTrigger_trk[j][c]->SetMarkerColor(j+1);
      s.HIByTrigger_trk[j][c]->SetFillColor(j+1);
      s.HIUsedByTrigger_trk[j][c]->SetLineColor(kBlack);
      s.HIUsedByTrigger_trk[j][c]->SetLineWidth(2);
      s.HIUsedByTrigger_trk[j][c]->SetMarkerColor(j+1);
      s.HIUsedByTrigger_trk[j][c]->SetMarkerSize(0);
      s.HIUsedByTrigger_trk[j][c]->SetFillColor(j+1);
      s.HIMaxtrkByTrigger[j][c]->SetLineColor(j+1);
      s.HIMaxtrkByTrigger[j][c]->SetLineWidth(1);
      s.HIMaxtrkByTrigger[j][c]->SetMarkerColor(j+1);
      s.HIMaxtrkByTrigger[j][c]->SetMarkerSize(0.8);
      s.HIMaxtrkByTrigger[j][c]->SetFillColor(j+1);
    }
  }
  
//******************************************************************************************************************
//************************************************JET TRIGGER PLOTS**********************************************************
//******************************************************************************************************************
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  TLine * l[8];
  TLatex * lat = new TLatex(1,1,"test");
  TLegend * leg;
  c1->SetLogy();
  for(int c = 0; c<s.nCentBins; c++){ 
    leg = new TLegend(0.5,0.6,0.9,0.9);
    s.HIJets[c]->Scale(100);
    s.HIJets[c]->SetMarkerSize(0);
    float Ymin = 0.00000000001;
    float Ymax = 10;
    c1->Clear();
    s.HIJets[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIJets[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers; i++) s.HIJetsByTrigger[i][c]->Draw("same");
    leg->Clear();
    leg->AddEntry(s.HIJets[c],"Jet Spectrum (x100)","l");
    leg->AddEntry(s.HIJetsByTrigger[0][c],"MB (I)","p");
    leg->AddEntry(s.HIJetsByTrigger[1][c],"jet 40 (II)","p");
    leg->AddEntry(s.HIJetsByTrigger[2][c],"jet 60 (III)","p");
    leg->AddEntry(s.HIJetsByTrigger[3][c],"jet 80 (IV)","p");
    leg->AddEntry(s.HIJetsByTrigger[4][c],"jet 100 (V)","p");
    leg->AddEntry((TObject*)0,"akVs4Calo Jets, |#eta|<2","");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c1->SaveAs(Form("plots/png/HIJets_FullSpectrum_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c1->SaveAs(Form("plots/pdf/HIJets_FullSpectrum_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c1->Clear();
    
    leg->SetX1NDC(0.5);
    leg->SetX2NDC(0.8);
    s.HIJets[c]->GetXaxis()->SetRangeUser(20,200);
    Ymin = 0.00000001;
    Ymax = 100;
    s.HIJets[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIJets[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers; i++) s.HIJetsByTrigger[i][c]->Draw("same");
    leg->Draw("same"); 
    for(int i = 0; i<4; i++) 
    {
      l[i] = new TLine(s.HItriggerBins[i+1],Ymin,s.HItriggerBins[i+1],Ymax); 
      l[i]->SetLineWidth(2);
      l[i]->SetLineStyle(2);
      l[i]->SetLineColor(1);
      l[i]->Draw("same");
    }
    lat->DrawLatex(45,Ymin*3,"I");
    lat->DrawLatex(65,Ymin*3,"II");
    lat->DrawLatex(85,Ymin*3,"III");
    lat->DrawLatex(105,Ymin*3,"IV");
    lat->DrawLatex(125,Ymin*3,"V");
    c1->SaveAs(Form("plots/png/HIJets_FullSpectrum_XZoom_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c1->SaveAs(Form("plots/pdf/HIJets_FullSpectrum_XZoom_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->Clear();
    for(int i = 0; i<4; i++) delete l[i]; 
    delete leg;
  }
  for(int c = 0; c<s.nCentBins;c++){
    c1->SetLogy(0);
    for(int i = 0; i<s.HInTriggers-1; i++) s.HIJetsByTrigger[s.HInTriggers-1-i][c]->Divide(s.HIJetsByTrigger[s.HInTriggers-2-i][c]);
    s.HIJetsByTrigger[1][c]->GetYaxis()->SetRangeUser(0,2);
    s.HIJetsByTrigger[1][c]->GetXaxis()->SetRangeUser(20,140);
    s.HIJetsByTrigger[1][c]->GetXaxis()->SetTitle("Leading jet p_{T}");
    s.HIJetsByTrigger[1][c]->Draw();
    s.HIJetsByTrigger[2][c]->Draw("same");
    s.HIJetsByTrigger[3][c]->Draw("same");
    s.HIJetsByTrigger[4][c]->Draw("same");
    leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->AddEntry(s.HIJetsByTrigger[1][c],"Jet40/MB","p");
    leg->AddEntry(s.HIJetsByTrigger[2][c],"Jet60/Jet40","p");
    leg->AddEntry(s.HIJetsByTrigger[3][c],"Jet80/Jet60","p");
    leg->AddEntry(s.HIJetsByTrigger[4][c],"Jet100/Jet80","p");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c1->SaveAs(Form("plots/png/HI_JetRelativeTurnOnes_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/HI_JetRelativeTurnOnes_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->Clear();
    delete leg;
  }
  
  c1->SetLogy();
  c1->SetLogx();
  for(int c = 0; c<s.nCentBins;c++){
    s.HI[c]->SetMarkerSize(0.8);
    s.HI[c]->GetYaxis()->SetRangeUser(TMath::Max(s.HI[20]->GetMinimum()/200.0,1e-15),s.HI[20]->GetMaximum()*10);
    s.HI[c]->Draw();
    s.HIUsedByTrigger[0][c]->SetFillColor(kGray);
    s.HIUsedByTrigger[4][c]->SetFillColor(kCyan+2);
    s.HIUsedByTrigger[3][c]->Add(s.HIUsedByTrigger[4][c]);
    s.HIUsedByTrigger[2][c]->Add(s.HIUsedByTrigger[3][c]);
    s.HIUsedByTrigger[1][c]->Add(s.HIUsedByTrigger[2][c]);
    s.HIUsedByTrigger[0][c]->Add(s.HIUsedByTrigger[1][c]);
    for(int i = 0; i<s.HInTriggers; i++) s.HIUsedByTrigger[i][c]->Draw("HIST same");
    s.HI[c]->Draw("sameaxis");
    s.HI[c]->Draw("same");
    leg = new TLegend(0.47,0.62,0.87,0.92);
    leg->AddEntry(s.HI[c],"PbPb track Spectrum","p");
    leg->AddEntry(s.HIUsedByTrigger[0][c],"MB trigger","f");
    leg->AddEntry(s.HIUsedByTrigger[1][c],"Jet40 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger[2][c],"Jet60 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger[3][c],"Jet80 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger[4][c],"Jet100 trigger","f");
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
     
    c1->SaveAs(Form("plots/png/HITrack_FullSpectrum_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/HITrack_FullSpectrum_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
  }

  c1->SetLogy(0);
  for(int c = 0; c<s.nCentBins;c++){
    c1->Clear();
    s.RAA[c] = (TH1D*)s.HI[c]->Clone(Form("RAA_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    s.RAA[c]->Scale(1/s.nColl[c]);
    s.RAA[c]->Divide(s.pp); 
    s.RAA[c]->Write();

    s.RAA[c]->GetYaxis()->SetRangeUser(0,1.2); 
    s.RAA[c]->GetYaxis()->SetTitle("R_{AA}"); 
    s.RAA[c]->Draw();
    leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->AddEntry((TObject*)0,"#sqrt{s_{NN}} = 5.02 TeV","");
    leg->Draw("same");
    TLine *line = new TLine(0.5,1,400,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    c1->SaveAs(Form("plots/png/RAA_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/RAA_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
    delete line;
  }
 
  TFile * KKRAA = TFile::Open("Krajczar_RAA_160216.root","read");
  TH1D * KKRAA_h[21];
  for(int c = 0; c<21; c++){
    if(c<20) KKRAA_h[c] = (TH1D*)KKRAA->Get(Form("hForPlotting_%d_%d",c*5,5+c*5));
    else     KKRAA_h[c] = (TH1D*)KKRAA->Get("hTrackPt_trkCorr_PbPb_copy1");
    c1->Clear();
    s.RAA[c]->Draw();
    KKRAA_h[c]->SetMarkerStyle(25); 
    KKRAA_h[c]->Draw("same"); 
    leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->AddEntry(s.RAA[c],"AB's RAA","p");
    leg->AddEntry(KKRAA_h[c],"KK's RAA","p");
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->AddEntry((TObject*)0,"#sqrt{s_{NN}} = 5.02 TeV","");
    leg->Draw("same");
    TLine *line = new TLine(0.5,1,400,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    c1->SaveAs(Form("plots/png/RAA_Comparison_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/RAA_Comparison_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
    delete line;
  }
  KKRAA->Close();

 
  delete lat;
  delete c1;
  
//************************************************************************************************************
//***********************************************TRACK TRIGGER PLOTS******************************************
//************************************************************************************************************
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  for(int c = 0; c<s.nCentBins; c++){
    c2->SetLogy();
    c2->Clear();
    s.HIMaxtrk[c]->SetMarkerSize(0);
    s.HIMaxtrkByTrigger[4][c]->SetMarkerColor(kOrange);
    s.HIMaxtrkByTrigger[4][c]->SetLineColor(kOrange);
    s.HIMaxtrk[c]->Scale(100);
    float Ymin = 0.000000001;
    float Ymax = 1000;
    s.HIMaxtrk[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIMaxtrk[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIMaxtrkByTrigger[i][c]->Draw("same");
    leg = new TLegend(0.5,0.6,0.8,0.9);
    leg->AddEntry(s.HIMaxtrk[c],"Leading Trk p_{T}Spectrum (x100)","l");
    leg->AddEntry(s.HIMaxtrkByTrigger[0][c],"MB (I)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[1][c],"Track 12 (II)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[2][c],"Track 18 (III)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[3][c],"Track 24 (IV)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[4][c],"Track 34 (V)","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[5][c],"Track 45 (VI)","p");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c2->SaveAs(Form("plots/png/HIMaxtrk_FullSpectrum_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c2->SaveAs(Form("plots/pdf/HIMaxtrk_FullSpectrum_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c2->Clear();
 
    s.HIMaxtrk[c]->GetXaxis()->SetRangeUser(0,70);
    Ymin = 0.00000001;
    Ymax = 1000;
    s.HIMaxtrk[c]->GetYaxis()->SetRangeUser(Ymin,Ymax);
    s.HIMaxtrk[c]->Draw("h");
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIMaxtrkByTrigger[i][c]->Draw("same");
    leg->Draw("same"); 
    for(int i = 0; i<5; i++) 
    {
      l[i] = new TLine(s.HItriggerBins_trk[i+1],Ymin,s.HItriggerBins_trk[i+1],Ymax); 
      l[i]->SetLineWidth(2);
      l[i]->SetLineStyle(2);
      l[i]->SetLineColor(1);
      l[i]->Draw("same");
    }
    lat->DrawLatex(8,Ymin*3,"I");
    lat->DrawLatex(16,Ymin*3,"II");
    lat->DrawLatex(22,Ymin*3,"III");
    lat->DrawLatex(28,Ymin*3,"IV");
    lat->DrawLatex(38,Ymin*3,"V");
    lat->DrawLatex(50,Ymin*3,"VI");
    c2->SaveAs(Form("plots/png/HIMaxtrk_FullSpectrum_XZoom_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    c2->SaveAs(Form("plots/pdf/HIMaxtrk_FullSpectrum_XZoom_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
    for(int i = 0; i<5; i++) delete l[i]; 
    delete leg;
  
    c2->SetLogy(0);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIMaxtrkByTrigger[i][c]->Rebin(2);
    for(int i = 0; i<s.HInTriggers_trk-1; i++) s.HIMaxtrkByTrigger[s.nTriggers_trk-1-i][c]->Divide(s.HIMaxtrkByTrigger[s.nTriggers_trk-2-i][c]);
    s.HIMaxtrkByTrigger[1][c]->GetYaxis()->SetRangeUser(0,2);
    s.HIMaxtrkByTrigger[1][c]->GetXaxis()->SetRangeUser(10,90);
    s.HIMaxtrkByTrigger[1][c]->GetXaxis()->SetTitle("Leading Track p_{T}");
    s.HIMaxtrkByTrigger[1][c]->Draw();
    s.HIMaxtrkByTrigger[2][c]->Draw("same");
    s.HIMaxtrkByTrigger[3][c]->Draw("same");
    s.HIMaxtrkByTrigger[4][c]->Draw("same");
    s.HIMaxtrkByTrigger[5][c]->Draw("same");
    leg = new TLegend(0.2,0.6,0.5,0.9);
    leg->AddEntry(s.HIMaxtrkByTrigger[1][c],"Trk12/MB","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[2][c],"Trk18/Trk12","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[3][c],"Trk24/Trk18","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[4][c],"Trk34/Trk24","p");
    leg->AddEntry(s.HIMaxtrkByTrigger[5][c],"Trk45/Trk34","p");
    leg->AddEntry((TObject*)0,Form("%d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
    c2->SaveAs(Form("plots/png/HITrackRelativeTurnOnes_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HITrackRelativeTurnOnes_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
  }

  c2->SetLogy();
  c2->SetLogx();
  for(int c = 0; c<s.nCentBins; c++){
    c2->Clear();
    s.HI_trk[c]->SetMarkerSize(0.8);
    s.HI_trk[c]->GetYaxis()->SetRangeUser(TMath::Max(s.HI_trk[20]->GetMinimum()/200.0,1e-15),s.HI_trk[20]->GetMaximum()*10);
    s.HI_trk[c]->Draw();
    s.HIUsedByTrigger_trk[0][c]->SetFillColor(kGray);
    s.HIUsedByTrigger_trk[3][c]->SetFillColor(kBlue);
    s.HIUsedByTrigger_trk[4][c]->SetFillColor(kViolet);
    s.HIUsedByTrigger_trk[5][c]->SetFillColor(kCyan+2);
    s.HIUsedByTrigger_trk[4][c]->Add(s.HIUsedByTrigger_trk[5][c]);
    s.HIUsedByTrigger_trk[3][c]->Add(s.HIUsedByTrigger_trk[4][c]);
    s.HIUsedByTrigger_trk[2][c]->Add(s.HIUsedByTrigger_trk[3][c]);
    s.HIUsedByTrigger_trk[1][c]->Add(s.HIUsedByTrigger_trk[2][c]);
    s.HIUsedByTrigger_trk[0][c]->Add(s.HIUsedByTrigger_trk[1][c]);
    for(int i = 0; i<s.HInTriggers_trk; i++) s.HIUsedByTrigger_trk[i][c]->Draw("HIST same");
    leg = new TLegend(0.5,0.6,0.8,0.9);
    s.HI_trk[c]->Draw("sameaxis");
    s.HI_trk[c]->Draw("same");
    leg->Clear();
    leg->AddEntry(s.HI_trk[c],"PbPb track Spectrum","p");
    leg->AddEntry(s.HIUsedByTrigger_trk[0][c],"MB trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[1][c],"Track12 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[2][c],"Track18 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[3][c],"Track24 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[4][c],"Track34 trigger","f");
    leg->AddEntry(s.HIUsedByTrigger_trk[5][c],"Track45 trigger","f");
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->Draw("same");
     
    c2->SaveAs(Form("plots/png/HITrack_FullSpectrum_trk_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HITrack_FullSpectrum_trk_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
  }
  
  c2->SetLogy(0);
  for(int c = 0; c<s.nCentBins;c++){
    c2->Clear();
    s.RAA_trk[c] = (TH1D*)s.HI_trk[c]->Clone(Form("RAA_trk_%d_%d",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    s.RAA_trk[c]->Scale(1/s.nColl[c]);
    s.RAA_trk[c]->Divide(s.pp_trk); 
    //s.RAA_trk[c]->Write();

    s.RAA_trk[c]->GetYaxis()->SetRangeUser(0,1.2); 
    s.RAA_trk[c]->GetYaxis()->SetTitle("R_{AA}"); 
    s.RAA_trk[c]->SetMarkerStyle(25); 
    s.RAA_trk[c]->Draw();
    s.RAA[c]->Draw("same");
    leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->AddEntry((TObject*)0,Form("|#eta|<1   %d-%d%s",s.lowCentBin[c]*5,s.highCentBin[c]*5,"%"),"");
    leg->AddEntry((TObject*)0,"#sqrt{s_{NN}} = 5.02 TeV","");
    leg->AddEntry(s.RAA[c],"Jet Triggers","p");
    leg->AddEntry(s.RAA_trk[c],"Track Triggers","p");
    leg->Draw("same");
    TLine *line = new TLine(0.5,1,400,1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("same");
    c1->SaveAs(Form("plots/png/RAA_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c1->SaveAs(Form("plots/pdf/RAA_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    delete leg;
    delete line;
  }

  c2->SetLogy(0);
  TH1D * HIjtVsTrk[s.nCentBins];
  for(int c = 0; c<s.nCentBins;c++){
    HIjtVsTrk[c] = (TH1D*)s.HI[c]->Clone(Form("HIJetVsTrk%d",c));
    HIjtVsTrk[c]->Divide(s.HI_trk[c]);
    HIjtVsTrk[c]->GetYaxis()->SetTitle("PbPb Jet triggers/PbPb Track triggers");
    HIjtVsTrk[c]->GetYaxis()->SetRangeUser(0.5,1.5);
    HIjtVsTrk[c]->Draw();
    c2->SaveAs(Form("plots/png/HI_jetVsTrkTriggers_%d_%d.png",s.lowCentBin[c]*5,s.highCentBin[c]*5));
    c2->SaveAs(Form("plots/pdf/HI_JetVsTrkTriggers_%d_%d.pdf",s.lowCentBin[c]*5,s.highCentBin[c]*5)); 
  }
  delete c2;

  return;
}
