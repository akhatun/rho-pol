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
void BeautifyHistoCosTheta(TH1* histogram);
void BeautifyHistoPhi(TH1* histogram);



//_______________________________________
void TrivialPlotting(){

  TFile* file   = new TFile("run2-MC/AxE/Corrected-run2.root", "READ");
  TH1F* MassCosTheta[16];
  TH1F* MassPhi[12];
  double integralCosTheta[16];
  double integralPhi[12];
  double errorCosTheta[16];
  double errorPhi[12];
  for (Int_t i = 0; i < 16; i++)
  {
    MassCosTheta[i] = (TH1F*) file->Get(Form("CorrCosTheta_%i", i));
    MassCosTheta[i]->Sumw2();
    integralCosTheta[i] = MassCosTheta[i]->Integral();
    errorCosTheta[i]    = TMath::Sqrt(integralCosTheta[i]);
  }
    for (Int_t i = 0; i < 12; i++)
  {
    MassPhi[i] = (TH1F*) file->Get(Form("CorrPhi_%i", i));
    MassPhi[i]->Sumw2();
    integralPhi[i] = MassPhi[i]->Integral();
    errorPhi[i]    = TMath::Sqrt(integralPhi[i]);

  }

  

  Double_t binEdgesCosTheta[] = { -1., -0.7, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0., 
                                  0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.7, 1.};
  Double_t binEdgesPhi[] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
  for(auto & i : binEdgesPhi) { i *= 2.*TMath::Pi()/12.; }       

  TH1F* PolarisationCosTheta = new TH1F("PolarisationCosTheta", "PolarisationCosTheta", 16, binEdgesCosTheta);
  TH1F* PolarisationPhi      = new TH1F("PolarisationPhi",      "PolarisationPhi",      12, binEdgesPhi     );
  for (int i = 0; i < 16; i++)
  {
    PolarisationCosTheta->SetBinContent(i+1, integralCosTheta[i] / (binEdgesCosTheta[i+1] - binEdgesCosTheta[i]));
    PolarisationCosTheta->SetBinError(  i+1, errorCosTheta[i]    / (binEdgesCosTheta[i+1] - binEdgesCosTheta[i]));
  }
  for (int i = 0; i < 12; i++)
  {
    PolarisationPhi->SetBinContent(i+1, integralPhi[i] / (binEdgesPhi[i+1] - binEdgesPhi[i]));
    PolarisationPhi->SetBinError(  i+1, errorPhi[i]    / (binEdgesPhi[i+1] - binEdgesPhi[i]));
  }
  

  BeautifyHistoCosTheta(PolarisationCosTheta);
  BeautifyHistoPhi(PolarisationPhi);

              





  TF1* costhetafit = new TF1("costhetafit", "[0]*((1.-[1]) + (3.*[1]-1.)*x*x )", -1., 1.);
  costhetafit->SetParameter(0, 1000);
  costhetafit->SetParameter(1, 0.);
  new TCanvas;
  BeautifyPad();
  TFitResultPtr r = PolarisationCosTheta->Fit("costhetafit", "RIS");
  TMatrixD cov = r->GetCorrelationMatrix();
  cov.Print();
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.2,0.72, Form("r_{00} = %.3f #pm %.3f (stat.)", costhetafit->GetParameter(1), costhetafit->GetParError(1)) );
  gPad->SaveAs("run2-MC/AxE/trivial-costheta-fit.pdf", "recreate");




  TF1* phifit = new TF1("phifit", "[0]*( 1.-[1]* TMath::Cos(2*x) )", 0., 2.*TMath::Pi());
  phifit->SetParameter(0, 1000);
  phifit->SetParameter(1, 0.);
  new TCanvas;
  BeautifyPad();
  TFitResultPtr t = PolarisationPhi->Fit("phifit", "RIS");
  TMatrixD covt = t->GetCorrelationMatrix();
  covt.Print();
  latex10->DrawLatex(0.2,0.72, Form("r_{1,-1} = %.3f #pm %.3f (stat.)", phifit->GetParameter(1), phifit->GetParError(1)) );
  gPad->SaveAs("run2-MC/AxE/trivial-phi-fit.pdf", "recreate");





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
void BeautifyHistoCosTheta(TH1* histogram){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("cos(#it{#theta})");
  histogram->GetYaxis()->SetTitle("AxE corrected integral of the inv. mass distribution");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  // histogram->GetXaxis()->SetRangeUser(0., 2.);
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 3 events corrected with Run 2 AxE" );
  gPad->SaveAs("run2-MC/AxE/trivial-costheta.pdf", "recreate");
  TFile* file2 = new TFile("run2-MC/CohHe.root", "READ");
  TH1F* MCgen = (TH1F*) file2->Get("CosThetaGen2H");
  MCgen->Draw("same");
  gPad->SaveAs("run2-MC/AxE/trivial-costheta-with-mc.pdf", "recreate");




}
//______________________________________________
void BeautifyHistoPhi(TH1* histogram){
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histogram);
  histogram->GetXaxis()->SetTitle("#varphi");
  histogram->GetYaxis()->SetTitle("AxE corrected integral of the inv. mass distribution");
  histogram->GetYaxis()->SetRangeUser(0., histogram->GetMaximum()*1.5);
  // histogram->GetXaxis()->SetRangeUser(0., 2.);
  histogram->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.36 TeV");
  latex10->DrawLatex(0.2,0.80, "Run 3 events corrected with Run 2 AxE" );
  gPad->SaveAs("run2-MC/AxE/trivial-phi.pdf", "recreate");
  TFile* file2 = new TFile("run2-MC/CohHe.root", "READ");
  TH1F* MCgen = (TH1F*) file2->Get("PhiGen2H");
  MCgen->Draw("same");
  gPad->SaveAs("run2-MC/AxE/trivial-phi-with-mc.pdf", "recreate");


}
