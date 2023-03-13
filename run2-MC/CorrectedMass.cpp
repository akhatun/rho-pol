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



//_______________________________________
void CorrectedMass(){

  TFile* file   = new TFile("run2-MC/AxE/AxE-run2.root", "READ");
  TFile* file2  = new TFile("AnalysisResults-run3-data-pass1.root", "READ");
  TH1F* AxECosTheta[16];
  TH1F* AxEPhi[12];
  TH1F* MassCosTheta[16];
  TH1F* MassPhi[12];
  for (Int_t i = 0; i < 16; i++)
  {
    MassCosTheta[i] = (TH1F*) file2->Get(Form("u-p-c-rho0-pol/hMassCosTheta_%i", i));
    AxECosTheta[i]  = (TH1F*) file ->Get(Form("AxECosTheta_%i", i));
    MassCosTheta[i]->Sumw2();
    AxECosTheta[i] ->Sumw2();
    MassCosTheta[i]->Divide(AxECosTheta[i]);
    AxEPhi[i]->SetName(Form("CorrCosTheta_%i", i));
  }
    for (Int_t i = 0; i < 12; i++)
  {
    MassPhi[i] = (TH1F*) file2->Get(Form("u-p-c-rho0-pol/hMassPhi_%i", i));
    AxEPhi[i]  = (TH1F*) file ->Get(Form("AxEPhi_%i", i));
    AxEPhi[i]->Sumw2();
    MassPhi[i]->Divide(AxEPhi[i]);
    AxEPhi[i] ->SetName(Form("CorrPhi_%i", i));
  }

  






  TFile *SavingFile = new TFile("run2-MC/AxE/Corrected-run2.root", "RECREATE");
  SavingFile->Close();
  for (Int_t i = 0; i < 16; i++)
  {
    MassCosTheta[i]->Write();
  }
    for (Int_t i = 0; i < 12; i++)
  {
    MassPhi[i]->Write();
  }


  for (Int_t i = 0; i < 16; i++)
  {
    BeautifyHistoCosTheta(MassCosTheta[i], i);
  }
    for (Int_t i = 0; i < 12; i++)
  {
    BeautifyHistoPhi(MassPhi[i], i);
  }









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
  latex10->DrawLatex(0.2,0.80, "Run 3 events corrected with Run 2 AxE" );
  // if(i == 0) latex10->DrawLatex(0.2,0.72, "-1 < cos#theta < -0.7");
  // if(i == 1) latex10->DrawLatex(0.2,0.72, "-0 7< cos#theta < -0.55");
  // if(i == 2) latex10->DrawLatex(0.2,0.72, "-0 55< cos#theta < -0.45");
  // if(i == 3) latex10->DrawLatex(0.2,0.72, "-0 55< cos#theta < -0.45");
  gPad->SaveAs(Form("run2-MC/AxE/Corr-costheta-%i.pdf", i), "recreate");


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
  latex10->DrawLatex(0.2,0.80, "Run 3 events corrected with Run 2 AxE" );
  // if(i == 0) latex10->DrawLatex(0.2,0.72, "-1 < cos#theta < -0.7");
  // if(i == 1) latex10->DrawLatex(0.2,0.72, "-0 7< cos#theta < -0.55");
  // if(i == 2) latex10->DrawLatex(0.2,0.72, "-0 55< cos#theta < -0.45");
  // if(i == 3) latex10->DrawLatex(0.2,0.72, "-0 55< cos#theta < -0.45");
  gPad->SaveAs(Form("run2-MC/AxE/Corr-phi-%i.pdf", i), "recreate");


}
