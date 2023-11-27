#include <TROOT.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TBox.h>
#include <TDatabasePDG.h>
#include <TWbox.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TMarker.h>
#include <stdio.h>
#include <iostream>

void BeautifyPad();
void BeautifyHisto(TH1* histogram);
void BeautifyHistoCosTheta(TH1* histogram, int i);
void BeautifyHistoPhi(TH1* histogram, int i);
void BeautifyHistoCosThetaRun3(TH1* histogram);
void BeautifyHistoPhiRun3(TH1* histogram);
void BeautifyHistoCosThetaRun2(TH1* histogram);
void BeautifyHistoPhiRun2(TH1* histogram);



//_______________________________________
void ComputeAxErun2(){

  TFile* file   = new TFile("run2-MC/CohHe-2.root", "READ");
  TH1F* MRecCosTheta[16];
  TH1F* MGenCosTheta[16];
  TH1F* MRecPhi[12];
  TH1F* MGenPhi[12];
  for (Int_t i = 0; i < 16; i++)
  {
    MRecCosTheta[i] = (TH1F*) file->Get(Form("MRecCosTheta_%i", i));
    MGenCosTheta[i] = (TH1F*) file->Get(Form("MGenCosTheta_%i", i));
    MRecCosTheta[i]->Sumw2();
    MGenCosTheta[i]->Sumw2();
    MRecCosTheta[i]->Divide(MGenCosTheta[i]);
    MRecCosTheta[i]->SetName(Form("AxECosTheta_%i", i));
  }
    for (Int_t i = 0; i < 12; i++)
  {
    MRecPhi[i] = (TH1F*) file->Get(Form("MRecPhi_%i", i));
    MGenPhi[i] = (TH1F*) file->Get(Form("MGenPhi_%i", i));
    MRecPhi[i]->Sumw2();
    MGenPhi[i]->Sumw2();
    MRecPhi[i]->Divide(MGenPhi[i]);
    MRecPhi[i]->SetName(Form("AxEPhi_%i", i));
  }

  
  TH1F* MRecCosThetaIntegrated;
  TH1F* MGenCosThetaIntegrated;
  TH1F* MRecPhiIntegrated;
  TH1F* MGenPhiIntegrated;
  MRecCosThetaIntegrated = (TH1F*) file->Get("CosThetaRecH");
  MGenCosThetaIntegrated = (TH1F*) file->Get("CosThetaGen2H");
  MRecCosThetaIntegrated->Sumw2();
  MGenCosThetaIntegrated->Sumw2();
  MRecCosThetaIntegrated->Rebin(10);
  MGenCosThetaIntegrated->Rebin(10);
  MRecCosThetaIntegrated->Divide(MGenCosThetaIntegrated);
  MRecCosThetaIntegrated->SetName("AxECosTheta");

  MRecPhiIntegrated = (TH1F*) file->Get("PhiRecH");
  MGenPhiIntegrated = (TH1F*) file->Get("PhiGen2H");
  MRecPhiIntegrated->Sumw2();
  MGenPhiIntegrated->Sumw2();
  MRecPhiIntegrated->Rebin(2);
  MGenPhiIntegrated->Rebin(2);
  MRecPhiIntegrated->Divide(MGenPhiIntegrated);
  MRecPhiIntegrated->SetName("AxEPhi");
  











  TFile *SavingFile = new TFile("run3-MC/AxE/AxE-run2.root", "RECREATE");
  MRecCosThetaIntegrated->Write();
  MRecPhiIntegrated->Write();
  for (Int_t i = 0; i < 16; i++)
  {
    MRecCosTheta[i]->Write();
  }
    for (Int_t i = 0; i < 12; i++)
  {
    MRecPhi[i]->Write();
  }
  SavingFile->Close();


  for (Int_t i = 0; i < 16; i++)
  {
    BeautifyHistoCosTheta(MRecCosTheta[i], i);
  }
    for (Int_t i = 0; i < 12; i++)
  {
    BeautifyHistoPhi(MRecPhi[i], i);
  }





  BeautifyHistoCosThetaRun2(MRecCosThetaIntegrated);
  BeautifyHistoPhiRun2(MRecPhiIntegrated);




}
//_______________________________________
void ComputeAxErun3(){

  TFile* fileSaved = new TFile("run3-MC/AnalysisResultsStarlightRun3.root", "READ");
  TH1F* MRecCosThetaIntegrated;
  TH1F* MGenCosThetaIntegrated;
  TH1F* MRecPhiIntegrated;
  TH1F* MGenPhiIntegrated;
  MRecCosThetaIntegrated = (TH1F*) fileSaved->Get("process-m-c/costhetaRhoReconstructed");
  MGenCosThetaIntegrated = (TH1F*) fileSaved->Get("process-m-c/costhetaRhoGenerated");
  MRecCosThetaIntegrated->Sumw2();
  MGenCosThetaIntegrated->Sumw2();
  MRecCosThetaIntegrated->Rebin(10);
  MGenCosThetaIntegrated->Rebin(10);
  MRecCosThetaIntegrated->Divide(MGenCosThetaIntegrated);
  MRecCosThetaIntegrated->SetName("AxECosTheta");

  MRecPhiIntegrated = (TH1F*) fileSaved->Get("process-m-c/phiRhoReconstructed");
  MGenPhiIntegrated = (TH1F*) fileSaved->Get("process-m-c/phiRhoGenerated");
  MRecPhiIntegrated->Sumw2();
  MGenPhiIntegrated->Sumw2();
  MRecPhiIntegrated->Rebin(30);
  MGenPhiIntegrated->Rebin(30);
  MRecPhiIntegrated->Divide(MGenPhiIntegrated);
  MRecPhiIntegrated->SetName("AxEPhi");
  






  TFile *SavingFile = new TFile("run3-MC/AxE/AxE-run3.root", "RECREATE");
  MRecCosThetaIntegrated->Write();
  MRecPhiIntegrated->Write();
  SavingFile->Close();
  BeautifyHistoCosThetaRun3(MRecCosThetaIntegrated);
  BeautifyHistoPhiRun3(MRecPhiIntegrated);



}
//______________________________________________
void BeautifyPad(){
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
}
//______________________________________________
void BeautifyHisto(TH1* histogram){
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.45);
  histogram->GetXaxis()->SetTitleSize(0.045);
  histogram->GetYaxis()->SetTitleSize(0.045);
  histogram->GetXaxis()->SetLabelSize(0.045);
  histogram->GetYaxis()->SetLabelSize(0.045);
  histogram->GetXaxis()->SetTitleFont(42);
  histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelFont(42);
  histogram->GetYaxis()->SetLabelFont(42);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");
}
//______________________________________________
void BeautifyHistoCosTheta(TH1* histogram, int i){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("M_{#mu#mu} [GeV/#it{c}^{2}]");
  histogram->GetYaxis()->SetTitle("AxE corrected invariant mass");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  histogram->GetXaxis()->SetRangeUser(0., 2.);
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 2 AxE" );
  if(i == 0)  latex10->DrawLatex(0.2,0.72, "-1    < cos#theta < -0.7");
  if(i == 1)  latex10->DrawLatex(0.2,0.72, "-0.7  < cos#theta < -0.55");
  if(i == 2)  latex10->DrawLatex(0.2,0.72, "-0.55 < cos#theta < -0.45");
  if(i == 3)  latex10->DrawLatex(0.2,0.72, "-0.45 < cos#theta < -0.35");
  if(i == 4)  latex10->DrawLatex(0.2,0.72, "-0.35 < cos#theta < -0.25");
  if(i == 5)  latex10->DrawLatex(0.2,0.72, "-0.25 < cos#theta < -0.15");
  if(i == 6)  latex10->DrawLatex(0.2,0.72, "-0.15 < cos#theta < -0.05");
  if(i == 7)  latex10->DrawLatex(0.2,0.72, "-0.05 < cos#theta <  0.");
  if(i == 8)  latex10->DrawLatex(0.2,0.72, " 0.   < cos#theta <  0.05");
  if(i == 9)  latex10->DrawLatex(0.2,0.72, " 0.05 < cos#theta <  0.15");
  if(i == 10) latex10->DrawLatex(0.2,0.72, " 0.15 < cos#theta <  0.25");
  if(i == 11) latex10->DrawLatex(0.2,0.72, " 0.25 < cos#theta <  0.35");
  if(i == 12) latex10->DrawLatex(0.2,0.72, " 0.35 < cos#theta <  0.45");
  if(i == 13) latex10->DrawLatex(0.2,0.72, " 0.45 < cos#theta <  0.55");
  if(i == 14) latex10->DrawLatex(0.2,0.72, " 0.55 < cos#theta <  0.7");
  if(i == 15) latex10->DrawLatex(0.2,0.72, " 0.7  < cos#theta <  1");
  gPad->SaveAs(Form("run3-MC/AxE/AxE-run2-costheta-%i.pdf", i), "recreate");


}
//______________________________________________
void BeautifyHistoPhi(TH1* histogram, int i){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("M_{#mu#mu} [GeV/#it{c}^{2}]");
  histogram->GetYaxis()->SetTitle("AxE corrected invariant mass");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  histogram->GetXaxis()->SetRangeUser(0., 2.);
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 2 AxE" );
  if(i == 0)  latex10->DrawLatex(0.2,0.72, "0    < #varphi < 1#times#frac{2#pi}{12}");
  if(i == 1)  latex10->DrawLatex(0.2,0.72, "1#times#frac{2#pi}{12}  < #varphi < 2#times#frac{2#pi}{12}");
  if(i == 2)  latex10->DrawLatex(0.2,0.72, "2#times#frac{2#pi}{12}  < #varphi < 3#times#frac{2#pi}{12}");
  if(i == 3)  latex10->DrawLatex(0.2,0.72, "3#times#frac{2#pi}{12}  < #varphi < 4#times#frac{2#pi}{12}");
  if(i == 4)  latex10->DrawLatex(0.2,0.72, "4#times#frac{2#pi}{12}  < #varphi < 5#times#frac{2#pi}{12}");
  if(i == 5)  latex10->DrawLatex(0.2,0.72, "5#times#frac{2#pi}{12}  < #varphi < 6#times#frac{2#pi}{12}");
  if(i == 6)  latex10->DrawLatex(0.2,0.72, "6#times#frac{2#pi}{12}  < #varphi < 7#times#frac{2#pi}{12}");
  if(i == 7)  latex10->DrawLatex(0.2,0.72, "7#times#frac{2#pi}{12}  < #varphi < 8#times#frac{2#pi}{12}");
  if(i == 8)  latex10->DrawLatex(0.2,0.72, "8#times#frac{2#pi}{12}  < #varphi < 9#times#frac{2#pi}{12}");
  if(i == 9)  latex10->DrawLatex(0.2,0.72, "9#times#frac{2#pi}{12}  < #varphi < 10#times#frac{2#pi}{12}");
  if(i == 10) latex10->DrawLatex(0.2,0.72, "10#times#frac{2#pi}{12} < #varphi < 11#times#frac{2#pi}{12}");
  if(i == 11) latex10->DrawLatex(0.2,0.72, "11#times#frac{2#pi}{12} < #varphi < 12#times#frac{2#pi}{12}");
  gPad->SaveAs(Form("run3-MC/AxE/AxE-run2-phi-%i.pdf", i), "recreate");


}
//______________________________________________
void BeautifyHistoCosThetaRun3(TH1* histogram){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("cos#it{#theta}");
  histogram->GetYaxis()->SetTitle("AxE for Run 3 MC");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  // histogram->GetXaxis()->SetRangeUser(0., 2.);
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 3 AxE" );
  gPad->SaveAs("run3-MC/AxE/AxE-run3-costheta.pdf", "recreate");
}
//______________________________________________
void BeautifyHistoPhiRun3(TH1* histogram){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("#varphi");
  histogram->GetYaxis()->SetTitle("AxE for Run 3 MC");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  // histogram->GetXaxis()->SetRangeUser(0., 2.);
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 3 AxE" );
  gPad->SaveAs("run3-MC/AxE/AxE-run3-phi.pdf", "recreate");
}
//______________________________________________
void BeautifyHistoCosThetaRun2(TH1* histogram){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("cos#it{#theta}");
  histogram->GetYaxis()->SetTitle("AxE for Run 2 MC");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  // histogram->GetXaxis()->SetRangeUser(0., 2.);
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 2 AxE" );
  gPad->SaveAs("run3-MC/AxE/AxE-run2-costheta.pdf", "recreate");
}
//______________________________________________
void BeautifyHistoPhiRun2(TH1* histogram){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("#varphi");
  histogram->GetYaxis()->SetTitle("AxE for Run 2 MC");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  histogram->GetXaxis()->SetRangeUser(0., 2.*TMath::Pi());
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 2 AxE" );
  gPad->SaveAs("run3-MC/AxE/AxE-run2-phi.pdf", "recreate");
}
