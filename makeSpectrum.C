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

  //RpPb bin shift correction
/*  TFile * binShift = TFile::Open("BinningAndResolutionCorrection_TrackTrigger.root","read");
  TH1D *hBinningAndResol_extendedPt = (TH1D*)f_binningAndResol_extendedPt->Get("hPt_pseudo2_copy1");
  hBinningAndResol_extendedPt->SetDirectory(0);
  binShift->Close();
  TH1D *pp_BinShift = (TH1D*)pp->Clone("pp_BinShift"); 
  for(int i = 2; i<hBinningAndResol_extendedPt->GetSize()+1; i++)
  {
    pp_BinShift->SetBinContent(i-1,hBinningAndResol_extendedPt->GetBinContent(i));
    pp_BinShift->SetBinError(i-1,hBinningAndResol_extendedPt->GetBinError(i));
  }
  for(int i = 1; i<pp_binShift->GetSize()+1; i++)
  {
    if(pp_binShift
  }*/

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
  
  pp->Divide(RpPb);
  float ppScale = 1/pp->GetBinContent(1);
  pp->Scale(ppScale);
  c1->SetLogy(0);
  pp->GetYaxis()->SetTitle("data/interp");
  pp->Draw(); 
  c1->SaveAs("plots/pp_dataVsInterp.png");
}
