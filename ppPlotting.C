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

void makePlotsPP(Settings s)
{
//******************************************************************************************************************
//************************************************JET TRIGGER PLOTS**********************************************************
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

  c1->SetLogy(0);
  for(int i = 0; i<s.nTriggers-1; i++) ppJetsByTrigger[s.nTriggers-1-i]->Divide(ppJetsByTrigger[s.nTriggers-2-i]);
  ppJetsByTrigger[1]->GetYaxis()->SetRangeUser(0,2);
  ppJetsByTrigger[1]->GetXaxis()->SetRangeUser(20,140);
  ppJetsByTrigger[1]->GetXaxis()->SetTitle("Leading jet p_{T}");
  ppJetsByTrigger[1]->Draw();
  ppJetsByTrigger[2]->Draw("same");
  ppJetsByTrigger[3]->Draw("same");
  TLegend * leg2 = new TLegend(0.2,0.6,0.5,0.9);
  leg2->AddEntry(ppJetsByTrigger[1],"Jet40/MB","p");
  leg2->AddEntry(ppJetsByTrigger[2],"Jet60/Jet40","p");
  leg2->AddEntry(ppJetsByTrigger[3],"Jet80/Jet60","p");
  leg2->Draw("same");
  c1->SaveAs("plots/JetRelativeTurnOnes.png");
  c1->SaveAs("plots/JetRelativeTurnOnes.pdf");
  
  c1->SetLogy();
  c1->SetLogx();
  pp->SetMarkerSize(0.8);
  pp->Draw();
  ppUsedByTrigger[0]->SetFillColor(kGray);
  ppUsedByTrigger[3]->SetFillColor(kCyan+2);
  ppUsedByTrigger[2]->Add(ppUsedByTrigger[3]);
  ppUsedByTrigger[1]->Add(ppUsedByTrigger[2]);
  ppUsedByTrigger[0]->Add(ppUsedByTrigger[1]);
  for(int i = 0; i<s.nTriggers; i++) ppUsedByTrigger[i]->Draw("HIST same");
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
  TH1D * jtVsTrk = (TH1D*)pp->Clone("jtVsTrkRatio");
  //RpPb bin shift correction
  
  //interpolation points
  float ppintyvals[s.ntrkBins+1] = {4.1202/0.5946*0.999, 2.8116/0.6211*1.003, 1.9778/0.6552*1.002, 1.4330/0.6984*0.994 , 1.0405/0.7219*1.002, 0.7719/0.7515*1.002, 0.5816/0.7809*1.002, 0.3893/0.825*1.001, 0.23381/0.866*1.000, 0.14395/0.901*1.004, 0.09103/0.925*1.005, 0.05937/0.965*0.997, 0.03906/0.984*0.998, 0.014787/1.023*1.062, 0.003806/1.052*1.042, 0.001181/1.056*1.033, 0.0004290/1.048*1.024, 0.0001787/1.054*1.018, 0.00008152/1.031*1.015, 0.00002216/1.023* 1.109, 0.000004653/1.036*1.061, 0.000001402/1.054*1.040, 0.0000003180/1.072*1.111, 0.00000006850/1.142*1.065, 0.00000001971/1.189*1.044, 0.000000006013/1.259*1.051, 0.000000001910/1.308*1.033, 0.0000000007181/1.342*1.024, 0.0000000002083/1.382*1.078,0.00000000005311/1.407*1.05,0.00000000001639/1.363*1.034386,0.000000000005354/1.381*1.047316,0.000000000001709/1.316*1.032760,0,0,0,0,0,0,0,0,0 };
  float ppintyvals_pPb[s.ntrkBins+1] = {4.1202*0.999, 2.8116*1.003, 1.9778*1.002, 1.4330*0.994 , 1.0405*1.002, 0.7719*1.002, 0.5816*1.002, 0.3893*1.001, 0.23381*1.000, 0.14395*1.004, 0.09103*1.005, 0.05937*0.997, 0.03906*0.998, 0.014787*1.062, 0.003806*1.042, 0.001181*1.033, 0.0004290*1.024, 0.0001787*1.018, 0.00008152*1.015, 0.00002216* 1.109, 0.000004653*1.061, 0.000001402*1.040, 0.0000003180*1.111, 0.00000006850*1.065, 0.00000001971*1.044, 0.000000006013*1.051, 0.000000001910*1.033, 0.0000000007181*1.024, 0.0000000002083*1.078,0.00000000005311*1.05,0.00000000001639*1.034386,0.000000000005354*1.047316,0.000000000001709*1.032760,0,0,0,0,0,0,0,0,0 };
  float ppintyvals_RpPb[s.ntrkBins+1] = {0.5946, 0.6211, 0.6552, 0.6984,0.7219,0.7515, 0.7809, 0.825, 0.866, 0.901, 0.925, 0.965, 0.984, 1.023, 1.052,1.056, 1.048, 1.054, 1.031, 1.023, 1.036, 1.054, 1.072, 1.142, 1.189, 1.259, 1.308, 1.342,1.382,1.407,1.363,1.381,1.316,0,0,0,0,0,0,0,0,0 };
  float ppintyerrRpPb[s.ntrkBins+1] = {0.0027/0.5946,0.0028/0.6211,0.0030/0.6552,0.0032/0.6984,0.0033/0.7219,0.0034/0.7515,0.0036/0.7809,0.003/0.825,0.003/0.866,0.003/0.901,0.003/0.925,0.003/0.965,0.003/0.984,0.002/1.023,0.002/1.052,0.002/1.056,0.002/1.048,0.003/1.054,0.003/1.031,0.003/1.023,0.005/1.036,0.008/1.054,0.001/1.072,0.002/1.142,0.001/1.189,0.002/1.259,0.002/1.308,0.004/1.342,0.004/1.382,0.008/1.407,0.013/1.363,0.019/1.381,0.030/1.316,0,0,0,0,0,0,0,0,0 };
  float ppintyerr_pPb[s.ntrkBins+1] = {0.0188/4.1202,0.0128/2.8116,0.0090/1.9778,0.0065/1.4330,0.0047/1.0405,0.0035/0.7719,0.0026/0.5816,0.0013/0.3893,0.00076/0.23381,0.00047/0.14395,0.00029/0.09103,0.00019/0.05937,0.00013/0.03906,0.000026/0.014787,0.000007/0.003806,0.000002/0.001181,0.0000009/0.0004290,0.0000004/0.0001787,0.00000025/0.00008152,0.00000025/0.00002216,0.00000006/0.000004653,0.000000023/0.000001402,0.000000011/0.0000003180,0.0000000004/0.00000006850,0.00000000011/0.00000001971,0.00000000002/0.000000006013,0.000000000009/0.000000001910,0.000000000004/0.0000000007181,0.0000000000020/0.0000000002083,0.0000000000007/0.00000000005311,0.00000000000031/0.00000000001639,0.00000000000016/0.000000000005354,0.000000000000073/0.000000000001709,0,0,0,0,0,0,0,0,0 };
  TH1D * RpPb = new TH1D("pp_interpolation",";pt;E#frac{d^{3}N}{dPtdyd#phi}",s.ntrkBins,s.xtrkbins);
  TH1D * RpPbSpectrum = new TH1D("RpPb_interp",";pt;RpPb",s.ntrkBins,s.xtrkbins);
  TH1D * pPbSpectrum = new TH1D("pPb_data",";pt;E#frac{d^{3}N}{dPtdyd#phi}",s.ntrkBins,s.xtrkbins);
  for(int i=1; i<s.ntrkBins+1; i++) RpPb->SetBinContent(i,ppintyvals[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) RpPb->SetBinError(i,TMath::Power(TMath::Abs(TMath::Power(ppintyerrRpPb[i-1],2)-TMath::Power(ppintyerr_pPb[i-1],2)),0.5)*ppintyvals[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) RpPbSpectrum->SetBinContent(i,ppintyvals_RpPb[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) RpPbSpectrum->SetBinError(i,ppintyvals_RpPb[i-1]*ppintyerrRpPb[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) pPbSpectrum->SetBinContent(i,ppintyvals_pPb[i-1]);
  for(int i=1; i<s.ntrkBins+1; i++) pPbSpectrum->SetBinError(i,ppintyvals_pPb[i-1]*ppintyerr_pPb[i-1]);
  
  float ppScale = RpPb->GetBinContent(17)/pp->GetBinContent(17);
  pp->Scale(ppScale); 

  TFile * binShift = TFile::Open("BinningAndResolutionCorrection_TrackTrigger.root","read");
  TH1D *hBinningAndResol_extendedPt = (TH1D*)binShift->Get("hPt_pseudo2_copy1");
  hBinningAndResol_extendedPt->SetDirectory(0);
  binShift->Close();
  TH1D *pp_BinShift = (TH1D*)pp->Clone("pp_BinShift"); 
  for(int i = 2; i<hBinningAndResol_extendedPt->GetSize()+1; i++)
  {
    pp_BinShift->SetBinContent(i-1,pp_BinShift->GetBinContent(i-1)/hBinningAndResol_extendedPt->GetBinContent(i));
    pp_BinShift->SetBinError(i-1,pp_BinShift->GetBinError(i-1)/hBinningAndResol_extendedPt->GetBinError(i));
  }
  TFile * outBinF = TFile::Open("ppSpectrum.root","update");
  pp_BinShift->Write();
  outBinF->Close();
  
  c1->SetLogy(0);
  pPbSpectrum->Divide(pp);
  pPbSpectrum->GetYaxis()->SetTitle("RpPb");
  pPbSpectrum->GetYaxis()->SetRangeUser(0.5,1.5);
  pPbSpectrum->Draw();
  c1->SaveAs("plots/pp_RpPb.png");
  
  pp->Divide(RpPb);
  pp->GetYaxis()->SetTitle("data/interp");
  pp->Draw(); 
  c1->SaveAs("plots/pp_dataVsInterp.png");
  
//************************************************************************************************************
//***********************************************TRACK TRIGGER PLOTS******************************************
//************************************************************************************************************
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  ppMaxtrkByTrigger[4]->SetMarkerColor(kOrange);
  ppMaxtrkByTrigger[4]->SetLineColor(kOrange);
  c2->SetLogy();
  ppMaxtrk->Scale(100);
  Ymin = 0.00000000001;
  Ymax = 10;
  ppMaxtrk->GetYaxis()->SetRangeUser(Ymin,Ymax);
  ppMaxtrk->Draw("h");
  for(int i = 0; i<s.nTriggers_trk; i++) ppMaxtrkByTrigger[i]->Draw("same");
  TLegend * leg_trk = new TLegend(0.6,0.6,0.9,0.9);
  leg_trk->AddEntry(ppMaxtrk,"Leading Trk p_{T}Spectrum (x100)","l");
  leg_trk->AddEntry(ppMaxtrkByTrigger[0],"MB (I)","p");
  leg_trk->AddEntry(ppMaxtrkByTrigger[1],"Track 18 (II)","p");
  leg_trk->AddEntry(ppMaxtrkByTrigger[2],"Track 24 (III)","p");
  leg_trk->AddEntry(ppMaxtrkByTrigger[3],"Track 34 (IV)","p");
  leg_trk->AddEntry(ppMaxtrkByTrigger[4],"Track 45 (V)","p");
  leg_trk->AddEntry(ppMaxtrkByTrigger[5],"Track 53 (VI)","p");
  leg_trk->Draw("same");
  c2->SaveAs("plots/ppMaxtrk_FullSpectrum.png"); 
  c2->SaveAs("plots/ppMaxtrk_FullSpectrum.pdf"); 

  ppMaxtrk->GetXaxis()->SetRangeUser(10,90);
  Ymin = 0.0000000001;
  Ymax = 1;
  ppMaxtrk->GetYaxis()->SetRangeUser(Ymin,Ymax);
  ppMaxtrk->Draw("h");
  for(int i = 0; i<s.nTriggers_trk; i++) ppMaxtrkByTrigger[i]->Draw("same");
  leg_trk->Draw("same"); 
  TLine * l_trk[5];
  for(int i = 0; i<5; i++) 
  {
    l_trk[i] = new TLine(s.triggerBins_trk[i+1],Ymin,s.triggerBins_trk[i+1],Ymax); 
    l_trk[i]->SetLineWidth(2);
    l_trk[i]->SetLineStyle(2);
    l_trk[i]->SetLineColor(1);
    l_trk[i]->Draw("same");
  }
  lat->DrawLatex(16,Ymin*3,"I");
  lat->DrawLatex(22,Ymin*3,"II");
  lat->DrawLatex(28,Ymin*3,"III");
  lat->DrawLatex(40,Ymin*3,"IV");
  lat->DrawLatex(49,Ymin*3,"V");
  lat->DrawLatex(57,Ymin*3,"VI");
  c2->SaveAs("plots/ppMaxtrk_FullSpectrum_XZoom.png"); 
  c2->SaveAs("plots/ppMaxtrk_FullSpectrum_XZoom.pdf"); 

  c2->SetLogy(0);
  for(int i = 0; i<s.nTriggers_trk; i++) ppMaxtrkByTrigger[i]->Rebin(2);
  for(int i = 0; i<s.nTriggers_trk-1; i++) ppMaxtrkByTrigger[s.nTriggers_trk-1-i]->Divide(ppMaxtrkByTrigger[s.nTriggers_trk-2-i]);
  ppMaxtrkByTrigger[1]->GetYaxis()->SetRangeUser(0,2);
  ppMaxtrkByTrigger[1]->GetXaxis()->SetRangeUser(10,90);
  ppMaxtrkByTrigger[1]->GetXaxis()->SetTitle("Leading Track p_{T}");
  ppMaxtrkByTrigger[1]->Draw();
  ppMaxtrkByTrigger[2]->Draw("same");
  ppMaxtrkByTrigger[3]->Draw("same");
  ppMaxtrkByTrigger[4]->Draw("same");
  ppMaxtrkByTrigger[5]->Draw("same");
  TLegend * leg2_trk = new TLegend(0.2,0.6,0.5,0.9);
  leg2_trk->AddEntry(ppMaxtrkByTrigger[1],"Trk18/MB","p");
  leg2_trk->AddEntry(ppMaxtrkByTrigger[2],"Trk24/Trk18","p");
  leg2_trk->AddEntry(ppMaxtrkByTrigger[3],"Trk34/Trk24","p");
  leg2_trk->AddEntry(ppMaxtrkByTrigger[4],"Trk45/Trk34","p");
  leg2_trk->AddEntry(ppMaxtrkByTrigger[5],"Trk53/Trk45","p");
  leg2_trk->Draw("same");
  c2->SaveAs("plots/TrackRelativeTurnOnes.png");
  c2->SaveAs("plots/TrackRelativeTurnOnes.pdf");
  
  c2->SetLogy();
  c2->SetLogx();
  pp_trk->SetMarkerSize(0.8);
  pp_trk->Draw();
  ppUsedByTrigger_trk[0]->SetFillColor(kGray);
  ppUsedByTrigger_trk[3]->SetFillColor(kBlue);
  ppUsedByTrigger_trk[4]->SetFillColor(kViolet);
  ppUsedByTrigger_trk[5]->SetFillColor(kCyan+2);
  ppUsedByTrigger_trk[4]->Add(ppUsedByTrigger_trk[5]);
  ppUsedByTrigger_trk[3]->Add(ppUsedByTrigger_trk[4]);
  ppUsedByTrigger_trk[2]->Add(ppUsedByTrigger_trk[3]);
  ppUsedByTrigger_trk[1]->Add(ppUsedByTrigger_trk[2]);
  ppUsedByTrigger_trk[0]->Add(ppUsedByTrigger_trk[1]);
  for(int i = 0; i<s.nTriggers_trk; i++) ppUsedByTrigger_trk[i]->Draw("HIST same");
  pp_trk->Draw("sameaxis");
  pp_trk->Draw("same");
  leg_trk->Clear();
  leg_trk->AddEntry(pp_trk,"pp track Spectrum","p");
  leg_trk->AddEntry(ppUsedByTrigger_trk[0],"MB trigger","f");
  leg_trk->AddEntry(ppUsedByTrigger_trk[1],"Track18 trigger","f");
  leg_trk->AddEntry(ppUsedByTrigger_trk[2],"Track24 trigger","f");
  leg_trk->AddEntry(ppUsedByTrigger_trk[3],"Track34 trigger","f");
  leg_trk->AddEntry(ppUsedByTrigger_trk[4],"Track45 trigger","f");
  leg_trk->AddEntry(ppUsedByTrigger_trk[5],"Track53 trigger","f");
  leg_trk->AddEntry((TObject*)0,"|#eta|<1","");
  leg_trk->Draw("same");
   
  c2->SaveAs("plots/ppTrack_FullSpectrum_trk.png");
  c2->SaveAs("plots/ppTrack_FullSpectrum_trk.pdf");

  c2->SetLogy(0);
  jtVsTrk->Divide(pp_trk);
  jtVsTrk->GetYaxis()->SetTitle("pp Jet triggers/pp Track triggers");
  jtVsTrk->GetYaxis()->SetRangeUser(0.5,1.5);
  jtVsTrk->Draw();
  c2->SaveAs("plots/pp_jetVsTrkTriggers.png");
  c2->SaveAs("plots/pp_JetVsTrkTriggers.pdf");

  delete c1;
  delete c2;
  delete lat;
  delete leg;
  delete leg2;
  delete leg_trk; 
  delete leg2_trk; 
  delete RpPb;
  delete RpPbSpectrum;
  delete pPbSpectrum;
  delete l;
  delete l_trk;
  return;
}
