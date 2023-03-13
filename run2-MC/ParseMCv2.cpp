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






//_____________________________________________________________________________
Double_t CosThetaHelicityFrame(TLorentzVector muonPositive, TLorentzVector muonNegative, TLorentzVector possibleJPsi )
{
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );

  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  
  //Translate the dimuon parameters in the dimuon rest frame
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);

  // Axes
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  Double_t CosThetaHE = zaxis.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaHE;

}
//_____________________________________________________________________________
Double_t PhiHelicityFrame(TLorentzVector muonPositive, TLorentzVector muonNegative, TLorentzVector possibleJPsi )
{

  //Half of the energy per pair of the colliding nucleons.
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  
  //Translate the dimuon parameters in the dimuon rest frame
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);

  // Axes
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
  return   phi;
}
//_______________________________________
void ParseMC(){

  TFile* file   = new TFile("run2-MC/AnalysisResults-coherentrho-run2-MC.root", "READ");
  TTree* ReconTree = (TTree*)file->Get("Rho0Central/Selected");
  TTree* GenerTree = (TTree*)file->Get("Rho0Central/Generated");



  Float_t pTrec;
  Float_t pTgen;
  ReconTree->SetBranchAddress("Pt_T",  &pTrec);
  // ReconTree->SetBranchAddress("Pt_MC_T",  &pTgen);
  Float_t Mrec;
  Float_t Mgen;
  ReconTree->SetBranchAddress("Mass_T",  &Mrec);
  GenerTree->SetBranchAddress("Mass_MC_T",  &Mgen);
  Float_t Yrec;
  Float_t Ygen;
  ReconTree->SetBranchAddress("Rapidity_T",  &Yrec);
  GenerTree->SetBranchAddress("Rapidity_MC_T",  &Ygen);
  Short_t  Qrec[2];
  ReconTree->SetBranchAddress("TrackQ_T",  &Qrec);
  Int_t  Qgen[2];
  GenerTree->SetBranchAddress("Q1_MC_T",  &Qgen[0]);
  GenerTree->SetBranchAddress("Q2_MC_T",  &Qgen[1]);
  Float_t  SinglePtrec[2];
  Float_t  SinglePXrec[2];
  Float_t  SinglePYrec[2];
  Float_t  SinglePZrec[2];
  Float_t  SinglePtgen[2];
  Float_t  SinglePtgen2[2];
  // ReconTree->SetBranchAddress("TrackPt_T",   &SinglePtrec);
  ReconTree->SetBranchAddress("TrackPx_T",   &SinglePXrec);
  ReconTree->SetBranchAddress("TrackPy_T",   &SinglePYrec);
  ReconTree->SetBranchAddress("TrackPz_T",   &SinglePZrec);
  ReconTree->SetBranchAddress("TrackPtGen_T",&SinglePtgen);
  // ReconTree->SetBranchAddress("TrackPtGen_T",&SinglePtgen2);
  GenerTree->SetBranchAddress("Pt1_MC_T",&SinglePtgen2[0]);
  GenerTree->SetBranchAddress("Pt2_MC_T",&SinglePtgen2[1]);;
  Float_t SingleEtarec[2];
  Float_t SingleEtagen[2];
  Float_t SingleEtagen2[2];
  ReconTree->SetBranchAddress("TrackEta_T",   &SingleEtarec);
  ReconTree->SetBranchAddress("TrackEtaGen_T",&SingleEtagen);
  // ReconTree->SetBranchAddress("TrackEtaGen_T",&SingleEtagen2);
  GenerTree->SetBranchAddress("Eta1_MC_T",&SingleEtagen2[0]);
  GenerTree->SetBranchAddress("Eta2_MC_T",&SingleEtagen2[1]);
  Float_t SinglePhirec[2];
  Float_t SinglePhigen[2];
  Float_t SinglePhigen2[2];
  ReconTree->SetBranchAddress("TrackPhi_T",   &SinglePhirec);
  ReconTree->SetBranchAddress("TrackPhiGen_T",&SinglePhigen);
  // ReconTree->SetBranchAddress("TrackPhiGen_T",&SinglePhigen2);
  GenerTree->SetBranchAddress("Phi1_MC_T",&SinglePhigen2[0]);
  GenerTree->SetBranchAddress("Phi2_MC_T",&SinglePhigen2[1]);

  Double_t CosThetaHErec;
  Double_t CosThetaHEgen;
  Double_t PhiHErec;
  Double_t PhiHEgen;

  Int_t nentriesRec = (Int_t) ReconTree->GetEntries();
  Int_t nentriesGen = (Int_t) GenerTree->GetEntries();


  TH1F *PtRecH    = new TH1F("PtRecH",    "PtRecH",    4000, 0, 20);
  TH1F *PtGenH    = new TH1F("PtGenH",    "PtGenH",    4000, 0, 20);
  TH1F *PtGen2H   = new TH1F("PtGen2H",   "PtGen2H",   4000, 0, 20);
  TH1F *PtGenSH   = new TH1F("PtGenSH",   "PtGenSH",   4000, 0, 20);
  TH1F *PtGenS2H  = new TH1F("PtGenS2H",  "PtGenS2H",  4000, 0, 20);
  TH1F *EtaGenSH  = new TH1F("EtaGenSH",  "EtaGenSH",  4000, 0, 20);
  TH1F *PhiGenSH  = new TH1F("PhiGenSH",  "PhiGenSH",  4000, 0, 20);
  TH1F *MRecH     = new TH1F("MRecH",     "MRecH",     4000, 0, 20);
  TH1F *MRecMineH = new TH1F("MRecMineH", "MRecMineH", 4000, 0, 20);
  TH1F *MGenH     = new TH1F("MGenH",     "MGenH",     4000, 0, 20);
  TH1F *MGenMineH = new TH1F("MGenMineH", "MGenMineH", 4000, 0, 20);
  TH1F *MGenMine2H= new TH1F("MGenMine2H","MGenMine2H",4000, 0, 20);

  TH1F *PhiRecMinusGenH      = new TH1F("PhiRecMinusGenH",      "PhiRecMinusGenH",      80, -2, 2);
  TH1F *CosThetaRecMinusGenH = new TH1F("CosThetaRecMinusGenH", "CosThetaRecMinusGenH", 80, -0.5, 0.5);


  TH1F *CosThetaRecH = new TH1F("CosThetaRecH", "CosThetaRecH", 250, -1., 1.);
  TH1F *CosThetaGenH = new TH1F("CosThetaGenH", "CosThetaGenH", 250, -1., 1.);
  TH1F *CosThetaGen2H= new TH1F("CosThetaGen2H", "CosThetaGen2H",250, -1., 1.);
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      250, -2.*TMath::Pi(), 4.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      250, -2.*TMath::Pi(), 4.*TMath::Pi());
  TH1F *PhiGen2H     = new TH1F("PhiGen2H",     "PhiGen2H",     250, -2.*TMath::Pi(), 4.*TMath::Pi());
  TH1F *PhiRecHH     = new TH1F("PhiRecHH",     "PhiRecHH",     500, -2.*TMath::Pi(), 2.*TMath::Pi());
  TH1F *PhiGenHH     = new TH1F("PhiGenHH",     "PhiGenHH",     500, -2.*TMath::Pi(), 2.*TMath::Pi());
  TH2F *CosThetaMigrationH = new TH2F("CosThetaMigrationH", "CosThetaMigrationH", 250, -1., 1., 250, -1., 1.);
  TH2F *PhiMigrationH      = new TH2F("PhiMigrationH",      "PhiMigrationH",      250, -2.*TMath::Pi(), 2.*TMath::Pi(), 250,-2.*TMath::Pi(), 2.*TMath::Pi());
  


  TH1F *MRecCosTheta_0  = new TH1F("MRecCosTheta_0",  "MRecCosTheta_0", 100, 0., 2.);
  TH1F *MRecCosTheta_1  = new TH1F("MRecCosTheta_1",  "MRecCosTheta_1", 100, 0., 2.);
  TH1F *MRecCosTheta_2  = new TH1F("MRecCosTheta_2",  "MRecCosTheta_2", 100, 0., 2.);
  TH1F *MRecCosTheta_3  = new TH1F("MRecCosTheta_3",  "MRecCosTheta_3", 100, 0., 2.);
  TH1F *MRecCosTheta_4  = new TH1F("MRecCosTheta_4",  "MRecCosTheta_4", 100, 0., 2.);
  TH1F *MRecCosTheta_5  = new TH1F("MRecCosTheta_5",  "MRecCosTheta_5", 100, 0., 2.);
  TH1F *MRecCosTheta_6  = new TH1F("MRecCosTheta_6",  "MRecCosTheta_6", 100, 0., 2.);
  TH1F *MRecCosTheta_7  = new TH1F("MRecCosTheta_7",  "MRecCosTheta_7", 100, 0., 2.);
  TH1F *MRecCosTheta_8  = new TH1F("MRecCosTheta_8",  "MRecCosTheta_8", 100, 0., 2.);
  TH1F *MRecCosTheta_9  = new TH1F("MRecCosTheta_9",  "MRecCosTheta_9", 100, 0., 2.);
  TH1F *MRecCosTheta_10 = new TH1F("MRecCosTheta_10", "MRecCosTheta_10", 100, 0., 2.);
  TH1F *MRecCosTheta_11 = new TH1F("MRecCosTheta_11", "MRecCosTheta_11", 100, 0., 2.);
  TH1F *MRecCosTheta_12 = new TH1F("MRecCosTheta_12", "MRecCosTheta_12", 100, 0., 2.);
  TH1F *MRecCosTheta_13 = new TH1F("MRecCosTheta_13", "MRecCosTheta_13", 100, 0., 2.);
  TH1F *MRecCosTheta_14 = new TH1F("MRecCosTheta_14", "MRecCosTheta_14", 100, 0., 2.);
  TH1F *MRecCosTheta_15 = new TH1F("MRecCosTheta_15", "MRecCosTheta_15", 100, 0., 2.);

  TH1F *MRecPhi_0  = new TH1F("MRecPhi_0",  "MRecPhi_0", 100, 0., 2.);
  TH1F *MRecPhi_1  = new TH1F("MRecPhi_1",  "MRecPhi_1", 100, 0., 2.);
  TH1F *MRecPhi_2  = new TH1F("MRecPhi_2",  "MRecPhi_2", 100, 0., 2.);
  TH1F *MRecPhi_3  = new TH1F("MRecPhi_3",  "MRecPhi_3", 100, 0., 2.);
  TH1F *MRecPhi_4  = new TH1F("MRecPhi_4",  "MRecPhi_4", 100, 0., 2.);
  TH1F *MRecPhi_5  = new TH1F("MRecPhi_5",  "MRecPhi_5", 100, 0., 2.);
  TH1F *MRecPhi_6  = new TH1F("MRecPhi_6",  "MRecPhi_6", 100, 0., 2.);
  TH1F *MRecPhi_7  = new TH1F("MRecPhi_7",  "MRecPhi_7", 100, 0., 2.);
  TH1F *MRecPhi_8  = new TH1F("MRecPhi_8",  "MRecPhi_8", 100, 0., 2.);
  TH1F *MRecPhi_9  = new TH1F("MRecPhi_9",  "MRecPhi_9", 100, 0., 2.);
  TH1F *MRecPhi_10 = new TH1F("MRecPhi_10", "MRecPhi_10", 100, 0., 2.);
  TH1F *MRecPhi_11 = new TH1F("MRecPhi_11", "MRecPhi_11", 100, 0., 2.);
  TH1F *MRecPhi_12 = new TH1F("MRecPhi_12", "MRecPhi_12", 100, 0., 2.);
  TH1F *MRecPhi_13 = new TH1F("MRecPhi_13", "MRecPhi_13", 100, 0., 2.);
  TH1F *MRecPhi_14 = new TH1F("MRecPhi_14", "MRecPhi_14", 100, 0., 2.);
  TH1F *MRecPhi_15 = new TH1F("MRecPhi_15", "MRecPhi_15", 100, 0., 2.);

  TH1F *MGenCosTheta_0  = new TH1F("MGenCosTheta_0",  "MGenCosTheta_0", 100, 0., 2.);
  TH1F *MGenCosTheta_1  = new TH1F("MGenCosTheta_1",  "MGenCosTheta_1", 100, 0., 2.);
  TH1F *MGenCosTheta_2  = new TH1F("MGenCosTheta_2",  "MGenCosTheta_2", 100, 0., 2.);
  TH1F *MGenCosTheta_3  = new TH1F("MGenCosTheta_3",  "MGenCosTheta_3", 100, 0., 2.);
  TH1F *MGenCosTheta_4  = new TH1F("MGenCosTheta_4",  "MGenCosTheta_4", 100, 0., 2.);
  TH1F *MGenCosTheta_5  = new TH1F("MGenCosTheta_5",  "MGenCosTheta_5", 100, 0., 2.);
  TH1F *MGenCosTheta_6  = new TH1F("MGenCosTheta_6",  "MGenCosTheta_6", 100, 0., 2.);
  TH1F *MGenCosTheta_7  = new TH1F("MGenCosTheta_7",  "MGenCosTheta_7", 100, 0., 2.);
  TH1F *MGenCosTheta_8  = new TH1F("MGenCosTheta_8",  "MGenCosTheta_8", 100, 0., 2.);
  TH1F *MGenCosTheta_9  = new TH1F("MGenCosTheta_9",  "MGenCosTheta_9", 100, 0., 2.);
  TH1F *MGenCosTheta_10 = new TH1F("MGenCosTheta_10", "MGenCosTheta_10", 100, 0., 2.);
  TH1F *MGenCosTheta_11 = new TH1F("MGenCosTheta_11", "MGenCosTheta_11", 100, 0., 2.);
  TH1F *MGenCosTheta_12 = new TH1F("MGenCosTheta_12", "MGenCosTheta_12", 100, 0., 2.);
  TH1F *MGenCosTheta_13 = new TH1F("MGenCosTheta_13", "MGenCosTheta_13", 100, 0., 2.);
  TH1F *MGenCosTheta_14 = new TH1F("MGenCosTheta_14", "MGenCosTheta_14", 100, 0., 2.);
  TH1F *MGenCosTheta_15 = new TH1F("MGenCosTheta_15", "MGenCosTheta_15", 100, 0., 2.);

  TH1F *MGenPhi_0  = new TH1F("MGenPhi_0",  "MGenPhi_0", 100, 0., 2.);
  TH1F *MGenPhi_1  = new TH1F("MGenPhi_1",  "MGenPhi_1", 100, 0., 2.);
  TH1F *MGenPhi_2  = new TH1F("MGenPhi_2",  "MGenPhi_2", 100, 0., 2.);
  TH1F *MGenPhi_3  = new TH1F("MGenPhi_3",  "MGenPhi_3", 100, 0., 2.);
  TH1F *MGenPhi_4  = new TH1F("MGenPhi_4",  "MGenPhi_4", 100, 0., 2.);
  TH1F *MGenPhi_5  = new TH1F("MGenPhi_5",  "MGenPhi_5", 100, 0., 2.);
  TH1F *MGenPhi_6  = new TH1F("MGenPhi_6",  "MGenPhi_6", 100, 0., 2.);
  TH1F *MGenPhi_7  = new TH1F("MGenPhi_7",  "MGenPhi_7", 100, 0., 2.);
  TH1F *MGenPhi_8  = new TH1F("MGenPhi_8",  "MGenPhi_8", 100, 0., 2.);
  TH1F *MGenPhi_9  = new TH1F("MGenPhi_9",  "MGenPhi_9", 100, 0., 2.);
  TH1F *MGenPhi_10 = new TH1F("MGenPhi_10", "MGenPhi_10", 100, 0., 2.);
  TH1F *MGenPhi_11 = new TH1F("MGenPhi_11", "MGenPhi_11", 100, 0., 2.);
  TH1F *MGenPhi_12 = new TH1F("MGenPhi_12", "MGenPhi_12", 100, 0., 2.);
  TH1F *MGenPhi_13 = new TH1F("MGenPhi_13", "MGenPhi_13", 100, 0., 2.);
  TH1F *MGenPhi_14 = new TH1F("MGenPhi_14", "MGenPhi_14", 100, 0., 2.);
  TH1F *MGenPhi_15 = new TH1F("MGenPhi_15", "MGenPhi_15", 100, 0., 2.);




  TH1F*** InvMassH;
  InvMassH = new TH1F**[24];
  for (Int_t i = 0; i < 24; i++) {
    InvMassH[i] = new TH1F*[24];
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j] = new TH1F(Form("InvMassH_%d_%d", i, j),Form("InvMassH_%d_%d", i, j),2000, 0, 20);
    }
  }



  TH1F*** InvMassHg;
  InvMassHg = new TH1F**[24];
  for (Int_t i = 0; i < 24; i++) {
    InvMassHg[i] = new TH1F*[24];
    for (Int_t j = 0; j < 24; j++) {
      InvMassHg[i][j] = new TH1F(Form("InvMassHg_%d_%d", i, j),Form("InvMassHg_%d_%d", i, j),2000, 0, 20);
    }
  }

  TH1F* FullInvMassH = new TH1F("FullInvMassH","FullInvMassH",2000, 0, 20);
  TH1F* FullInvMassHgen = new TH1F("FullInvMassHgen","FullInvMassHgen",2000, 0, 20);


  for (Int_t i=0; i<nentriesGen; i++) {
    GenerTree->GetEntry(i);
    if( (Ygen > -0.8 && Ygen < 0.8) ) {
      MGenH->Fill(Mgen);
      TLorentzVector piPlusMC2, piMinusMC2, rhoMC2;
      Double_t EnergyGen2[2];
      Float_t pXgen2[2], pYgen2[2], pZgen2[2];
      pXgen2[0] = SinglePtgen2[0] * TMath::Cos(SinglePhigen2[0]);
      pXgen2[1] = SinglePtgen2[1] * TMath::Cos(SinglePhigen2[1]);
      pYgen2[0] = SinglePtgen2[0] * TMath::Sin(SinglePhigen2[0]);
      pYgen2[1] = SinglePtgen2[1] * TMath::Sin(SinglePhigen2[1]);
      pZgen2[0] = SinglePtgen2[0] * TMath::SinH(SingleEtagen2[0]);
      pZgen2[1] = SinglePtgen2[1] * TMath::SinH(SingleEtagen2[1]);
      EnergyGen2[0] = TMath::Sqrt( pXgen2[0]*pXgen2[0] + pYgen2[0]*pYgen2[0] + pZgen2[0]*pZgen2[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
      EnergyGen2[1] = TMath::Sqrt( pXgen2[1]*pXgen2[1] + pYgen2[1]*pYgen2[1] + pZgen2[1]*pZgen2[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
      if(Qgen[1] < -0.5 && Qgen[0] > 0.5){
        piPlusMC2.SetPxPyPzE(  pXgen2[1],pYgen2[1],pZgen2[1],EnergyGen2[1]);
        piMinusMC2.SetPxPyPzE( pXgen2[0],pYgen2[0],pZgen2[0],EnergyGen2[0]);
      } else if(Qgen[0] < -0.5 && Qgen[1] > 0.5){
        piPlusMC2.SetPxPyPzE(  pXgen2[0],pYgen2[0],pZgen2[0],EnergyGen2[0]);
        piMinusMC2.SetPxPyPzE( pXgen2[1],pYgen2[1],pZgen2[1],EnergyGen2[1]);
      }
      rhoMC2 += piPlusMC2;
      rhoMC2 += piMinusMC2;
      Double_t CosThetaHEgen2= CosThetaHelicityFrame(piPlusMC2,piMinusMC2,rhoMC2);
      Double_t PhiHEgen2     = TMath::Pi() + PhiHelicityFrame(piPlusMC2, piMinusMC2, rhoMC2);
      CosThetaGen2H       ->Fill(CosThetaHEgen2);
      PhiGen2H            ->Fill(PhiHEgen2);

      if (CosThetaHEgen2>-1.    && CosThetaHEgen2<-0.7)  {MGenCosTheta_0->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>-0.7   && CosThetaHEgen2<-0.55) {MGenCosTheta_1->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>-0.55  && CosThetaHEgen2<-0.45) {MGenCosTheta_2->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>-0.45  && CosThetaHEgen2<-0.35) {MGenCosTheta_3->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>-0.35  && CosThetaHEgen2<-0.25) {MGenCosTheta_4->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>-0.25  && CosThetaHEgen2<-0.15) {MGenCosTheta_5->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>-0.15  && CosThetaHEgen2<-0.05) {MGenCosTheta_6->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>-0.05  && CosThetaHEgen2<0.)    {MGenCosTheta_7->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.     && CosThetaHEgen2<0.05)  {MGenCosTheta_8->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.05   && CosThetaHEgen2<0.15)  {MGenCosTheta_9->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.15   && CosThetaHEgen2<0.25)  {MGenCosTheta_10->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.25   && CosThetaHEgen2<0.35)  {MGenCosTheta_11->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.35   && CosThetaHEgen2<0.45)  {MGenCosTheta_12->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.45   && CosThetaHEgen2<0.55)  {MGenCosTheta_13->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.55   && CosThetaHEgen2<0.7)   {MGenCosTheta_14->Fill(rhoMC2.M());}
      if (CosThetaHEgen2>0.7    && CosThetaHEgen2<1.)    {MGenCosTheta_15->Fill(rhoMC2.M());}
              
              
      if (PhiHEgen2>2.*TMath::Pi()/12. * 0.  && PhiHEgen2<2.*TMath::Pi()/12. * 1.)  {MGenPhi_0->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 1.  && PhiHEgen2<2.*TMath::Pi()/12. * 2.)  {MGenPhi_1->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 2.  && PhiHEgen2<2.*TMath::Pi()/12. * 3.)  {MGenPhi_2->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 3.  && PhiHEgen2<2.*TMath::Pi()/12. * 4.)  {MGenPhi_3->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 4.  && PhiHEgen2<2.*TMath::Pi()/12. * 5.)  {MGenPhi_4->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 5.  && PhiHEgen2<2.*TMath::Pi()/12. * 6.)  {MGenPhi_5->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 6.  && PhiHEgen2<2.*TMath::Pi()/12. * 7.)  {MGenPhi_6->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 7.  && PhiHEgen2<2.*TMath::Pi()/12. * 8.)  {MGenPhi_7->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 8.  && PhiHEgen2<2.*TMath::Pi()/12. * 9.)  {MGenPhi_8->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 9.  && PhiHEgen2<2.*TMath::Pi()/12. * 10.) {MGenPhi_9->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 10. && PhiHEgen2<2.*TMath::Pi()/12. * 11.) {MGenPhi_10->Fill(rhoMC2.M());}
      if (PhiHEgen2>2.*TMath::Pi()/12. * 11. && PhiHEgen2<2.*TMath::Pi()/12. * 12.) {MGenPhi_11->Fill(rhoMC2.M());}


    }
  }

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);
    TLorentzVector piPlus, piMinus, rho;
    TLorentzVector piPlusMC, piMinusMC, rhoMC;
    // TLorentzVector piPlusMC2, piMinusMC2, rhoMC2;
    SinglePtrec[0] = TMath::Sqrt(SinglePXrec[0]*SinglePXrec[0] + SinglePYrec[0]*SinglePYrec[0]);
    SinglePtrec[1] = TMath::Sqrt(SinglePXrec[1]*SinglePXrec[1] + SinglePYrec[1]*SinglePYrec[1]);
    Double_t EnergyRec[2], EnergyGen[2];
    // Double_t EnergyRec[2], EnergyGen[2], EnergyGen2[2];
    EnergyRec[0] = TMath::Sqrt( SinglePXrec[0]*SinglePXrec[0] + SinglePYrec[0]*SinglePYrec[0] + SinglePZrec[0]*SinglePZrec[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyRec[1] = TMath::Sqrt( SinglePXrec[1]*SinglePXrec[1] + SinglePYrec[1]*SinglePYrec[1] + SinglePZrec[1]*SinglePZrec[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    Float_t pXgen[2], pYgen[2], pZgen[2];
    pXgen[0] = SinglePtgen[0] * TMath::Cos(SinglePhigen[0]);
    pXgen[1] = SinglePtgen[1] * TMath::Cos(SinglePhigen[1]);
    pYgen[0] = SinglePtgen[0] * TMath::Sin(SinglePhigen[0]);
    pYgen[1] = SinglePtgen[1] * TMath::Sin(SinglePhigen[1]);
    pZgen[0] = SinglePtgen[0] * TMath::SinH(SingleEtagen[0]);
    pZgen[1] = SinglePtgen[1] * TMath::SinH(SingleEtagen[1]);
    // Float_t pXgen2[2], pYgen2[2], pZgen2[2];
    // pXgen2[0] = SinglePtgen2[0] * TMath::Cos(SinglePhigen2[0]);
    // pXgen2[1] = SinglePtgen2[1] * TMath::Cos(SinglePhigen2[1]);
    // pYgen2[0] = SinglePtgen2[0] * TMath::Sin(SinglePhigen2[0]);
    // pYgen2[1] = SinglePtgen2[1] * TMath::Sin(SinglePhigen2[1]);
    // pZgen2[0] = SinglePtgen2[0] * TMath::SinH(SingleEtagen2[0]);
    // pZgen2[1] = SinglePtgen2[1] * TMath::SinH(SingleEtagen2[1]);
    // EnergyGen[0] = TMath::Sqrt( SinglePtgen[0]*SinglePtgen[0] + pZgen[0]*pZgen[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    // EnergyGen[1] = TMath::Sqrt( SinglePtgen[1]*SinglePtgen[1] + pZgen[1]*pZgen[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyGen[0] = TMath::Sqrt( pXgen[0]*pXgen[0] + pYgen[0]*pYgen[0] + pZgen[0]*pZgen[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyGen[1] = TMath::Sqrt( pXgen[1]*pXgen[1] + pYgen[1]*pYgen[1] + pZgen[1]*pZgen[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    // EnergyGen2[0] = TMath::Sqrt( pXgen2[0]*pXgen2[0] + pYgen2[0]*pYgen2[0] + pZgen2[0]*pZgen2[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    // EnergyGen2[1] = TMath::Sqrt( pXgen2[1]*pXgen2[1] + pYgen2[1]*pYgen2[1] + pZgen2[1]*pZgen2[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    if(Qrec[0] < -0.5 && Qrec[1] > 0.5){
      piPlus.SetPxPyPzE(    SinglePXrec[1],SinglePYrec[1],SinglePZrec[1],EnergyRec[1]);
      piMinus.SetPxPyPzE(   SinglePXrec[0],SinglePYrec[0],SinglePZrec[0],EnergyRec[0]);
      piPlusMC.SetPxPyPzE(  pXgen[1],pYgen[1],pZgen[1],EnergyGen[1]);
      piMinusMC.SetPxPyPzE( pXgen[0],pYgen[0],pZgen[0],EnergyGen[0]);
      // piPlusMC2.SetPxPyPzE(  pXgen2[1],pYgen2[1],pZgen2[1],EnergyGen2[1]);
      // piMinusMC2.SetPxPyPzE( pXgen2[0],pYgen2[0],pZgen2[0],EnergyGen2[0]);
    } else if(Qrec[1] < -0.5 && Qrec[0] > 0.5){
      piPlus.SetPxPyPzE(  SinglePXrec[0],SinglePYrec[0],SinglePZrec[0],EnergyRec[0]);
      piMinus.SetPxPyPzE( SinglePXrec[1],SinglePYrec[1],SinglePZrec[1],EnergyRec[1]);
      piPlusMC.SetPxPyPzE(  pXgen[0],pYgen[0],pZgen[0],EnergyGen[0]);
      piMinusMC.SetPxPyPzE( pXgen[1],pYgen[1],pZgen[1],EnergyGen[1]);
      // piPlusMC2.SetPxPyPzE(  pXgen2[0],pYgen2[0],pZgen2[0],EnergyGen2[0]);
      // piMinusMC2.SetPxPyPzE( pXgen2[1],pYgen2[1],pZgen2[1],EnergyGen2[1]);
    }
    rho   += piPlus;
    rho   += piMinus;
    rhoMC += piPlusMC;
    rhoMC += piMinusMC;
    // rhoMC2 += piPlusMC2;
    // rhoMC2 += piMinusMC2;



      if( (pTrec > 0. && pTrec < 0.20) && (Yrec > -0.8 && Yrec < 0.8) ) {

        MRecH->Fill(Mrec);
        MRecMineH->Fill(rho.M());
        MGenMineH->Fill(rhoMC.M());
        PtGenH->Fill(rhoMC.Pt());

        Double_t CosThetaHErec = CosThetaHelicityFrame(piPlus,   piMinus,   rho);
        Double_t CosThetaHEgen = CosThetaHelicityFrame(piPlusMC, piMinusMC, rhoMC);
        // Double_t CosThetaHEgen2= CosThetaHelicityFrame(piPlusMC2,piMinusMC2,rhoMC2);
        Double_t PhiHErec      = TMath::Pi() + PhiHelicityFrame(piPlus,   piMinus,   rho);
        Double_t PhiHEgen      = TMath::Pi() + PhiHelicityFrame(piPlusMC, piMinusMC, rhoMC);
        // Double_t PhiHEgen2     = TMath::Pi() + PhiHelicityFrame(piPlusMC2, piMinusMC2, rhoMC2);
        CosThetaRecH        ->Fill(CosThetaHErec);
        CosThetaGenH        ->Fill(CosThetaHEgen);
        // CosThetaGen2H       ->Fill(CosThetaHEgen2);
        PhiRecH             ->Fill(PhiHErec);
        PhiGenH             ->Fill(PhiHEgen);
        // PhiGen2H            ->Fill(PhiHEgen2);
        PhiMigrationH       ->Fill(PhiHEgen,PhiHErec);
        CosThetaMigrationH  ->Fill(CosThetaHEgen,CosThetaHErec);
        PhiRecMinusGenH     ->Fill(PhiHErec-PhiHEgen);
        CosThetaRecMinusGenH->Fill(CosThetaHErec-CosThetaHEgen);



        if (CosThetaHErec>-1.    && CosThetaHErec<-0.7)  {MRecCosTheta_0->Fill(rho.M());}
        if (CosThetaHErec>-0.7   && CosThetaHErec<-0.55) {MRecCosTheta_1->Fill(rho.M());}
        if (CosThetaHErec>-0.55  && CosThetaHErec<-0.45) {MRecCosTheta_2->Fill(rho.M());}
        if (CosThetaHErec>-0.45  && CosThetaHErec<-0.35) {MRecCosTheta_3->Fill(rho.M());}
        if (CosThetaHErec>-0.35  && CosThetaHErec<-0.25) {MRecCosTheta_4->Fill(rho.M());}
        if (CosThetaHErec>-0.25  && CosThetaHErec<-0.15) {MRecCosTheta_5->Fill(rho.M());}
        if (CosThetaHErec>-0.15  && CosThetaHErec<-0.05) {MRecCosTheta_6->Fill(rho.M());}
        if (CosThetaHErec>-0.05  && CosThetaHErec<0.)    {MRecCosTheta_7->Fill(rho.M());}
        if (CosThetaHErec>0.     && CosThetaHErec<0.05)  {MRecCosTheta_8->Fill(rho.M());}
        if (CosThetaHErec>0.05   && CosThetaHErec<0.15)  {MRecCosTheta_9->Fill(rho.M());}
        if (CosThetaHErec>0.15   && CosThetaHErec<0.25)  {MRecCosTheta_10->Fill(rho.M());}
        if (CosThetaHErec>0.25   && CosThetaHErec<0.35)  {MRecCosTheta_11->Fill(rho.M());}
        if (CosThetaHErec>0.35   && CosThetaHErec<0.45)  {MRecCosTheta_12->Fill(rho.M());}
        if (CosThetaHErec>0.45   && CosThetaHErec<0.55)  {MRecCosTheta_13->Fill(rho.M());}
        if (CosThetaHErec>0.55   && CosThetaHErec<0.7)   {MRecCosTheta_14->Fill(rho.M());}
        if (CosThetaHErec>0.7    && CosThetaHErec<1.)    {MRecCosTheta_15->Fill(rho.M());}
              
              
        if (PhiHErec>2.*TMath::Pi()/12. * 0.  && PhiHErec<2.*TMath::Pi()/12. * 1.)  {MRecPhi_0->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 1.  && PhiHErec<2.*TMath::Pi()/12. * 2.)  {MRecPhi_1->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 2.  && PhiHErec<2.*TMath::Pi()/12. * 3.)  {MRecPhi_2->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 3.  && PhiHErec<2.*TMath::Pi()/12. * 4.)  {MRecPhi_3->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 4.  && PhiHErec<2.*TMath::Pi()/12. * 5.)  {MRecPhi_4->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 5.  && PhiHErec<2.*TMath::Pi()/12. * 6.)  {MRecPhi_5->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 6.  && PhiHErec<2.*TMath::Pi()/12. * 7.)  {MRecPhi_6->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 7.  && PhiHErec<2.*TMath::Pi()/12. * 8.)  {MRecPhi_7->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 8.  && PhiHErec<2.*TMath::Pi()/12. * 9.)  {MRecPhi_8->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 9.  && PhiHErec<2.*TMath::Pi()/12. * 10.) {MRecPhi_9->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 10. && PhiHErec<2.*TMath::Pi()/12. * 11.) {MRecPhi_10->Fill(rho.M());}
        if (PhiHErec>2.*TMath::Pi()/12. * 11. && PhiHErec<2.*TMath::Pi()/12. * 12.) {MRecPhi_11->Fill(rho.M());}













        // if (CosThetaHEgen>-1.    && CosThetaHEgen<-0.7)  {MGenCosTheta_0->Fill(rhoMC.M());}
        // if (CosThetaHEgen>-0.7   && CosThetaHEgen<-0.55) {MGenCosTheta_1->Fill(rhoMC.M());}
        // if (CosThetaHEgen>-0.55  && CosThetaHEgen<-0.45) {MGenCosTheta_2->Fill(rhoMC.M());}
        // if (CosThetaHEgen>-0.45  && CosThetaHEgen<-0.35) {MGenCosTheta_3->Fill(rhoMC.M());}
        // if (CosThetaHEgen>-0.35  && CosThetaHEgen<-0.25) {MGenCosTheta_4->Fill(rhoMC.M());}
        // if (CosThetaHEgen>-0.25  && CosThetaHEgen<-0.15) {MGenCosTheta_5->Fill(rhoMC.M());}
        // if (CosThetaHEgen>-0.15  && CosThetaHEgen<-0.05) {MGenCosTheta_6->Fill(rhoMC.M());}
        // if (CosThetaHEgen>-0.05  && CosThetaHEgen<0.)    {MGenCosTheta_7->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.     && CosThetaHEgen<0.05)  {MGenCosTheta_8->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.05   && CosThetaHEgen<0.15)  {MGenCosTheta_9->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.15   && CosThetaHEgen<0.25)  {MGenCosTheta_10->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.25   && CosThetaHEgen<0.35)  {MGenCosTheta_11->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.35   && CosThetaHEgen<0.45)  {MGenCosTheta_12->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.45   && CosThetaHEgen<0.55)  {MGenCosTheta_13->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.55   && CosThetaHEgen<0.7)   {MGenCosTheta_14->Fill(rhoMC.M());}
        // if (CosThetaHEgen>0.7    && CosThetaHEgen<1.)    {MGenCosTheta_15->Fill(rhoMC.M());}
              
              
        // if (PhiHEgen>2.*TMath::Pi()/12. * 0.  && PhiHEgen<2.*TMath::Pi()/12. * 1.)  {MGenPhi_0->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 1.  && PhiHEgen<2.*TMath::Pi()/12. * 2.)  {MGenPhi_1->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 2.  && PhiHEgen<2.*TMath::Pi()/12. * 3.)  {MGenPhi_2->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 3.  && PhiHEgen<2.*TMath::Pi()/12. * 4.)  {MGenPhi_3->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 4.  && PhiHEgen<2.*TMath::Pi()/12. * 5.)  {MGenPhi_4->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 5.  && PhiHEgen<2.*TMath::Pi()/12. * 6.)  {MGenPhi_5->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 6.  && PhiHEgen<2.*TMath::Pi()/12. * 7.)  {MGenPhi_6->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 7.  && PhiHEgen<2.*TMath::Pi()/12. * 8.)  {MGenPhi_7->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 8.  && PhiHEgen<2.*TMath::Pi()/12. * 9.)  {MGenPhi_8->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 9.  && PhiHEgen<2.*TMath::Pi()/12. * 10.) {MGenPhi_9->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 10. && PhiHEgen<2.*TMath::Pi()/12. * 11.) {MGenPhi_10->Fill(rhoMC.M());}
        // if (PhiHEgen>2.*TMath::Pi()/12. * 11. && PhiHEgen<2.*TMath::Pi()/12. * 12.) {MGenPhi_11->Fill(rhoMC.M());}






        // if (CosThetaHEgen2>-1.    && CosThetaHEgen2<-0.7)  {MGenCosTheta_0->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>-0.7   && CosThetaHEgen2<-0.55) {MGenCosTheta_1->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>-0.55  && CosThetaHEgen2<-0.45) {MGenCosTheta_2->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>-0.45  && CosThetaHEgen2<-0.35) {MGenCosTheta_3->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>-0.35  && CosThetaHEgen2<-0.25) {MGenCosTheta_4->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>-0.25  && CosThetaHEgen2<-0.15) {MGenCosTheta_5->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>-0.15  && CosThetaHEgen2<-0.05) {MGenCosTheta_6->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>-0.05  && CosThetaHEgen2<0.)    {MGenCosTheta_7->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.     && CosThetaHEgen2<0.05)  {MGenCosTheta_8->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.05   && CosThetaHEgen2<0.15)  {MGenCosTheta_9->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.15   && CosThetaHEgen2<0.25)  {MGenCosTheta_10->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.25   && CosThetaHEgen2<0.35)  {MGenCosTheta_11->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.35   && CosThetaHEgen2<0.45)  {MGenCosTheta_12->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.45   && CosThetaHEgen2<0.55)  {MGenCosTheta_13->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.55   && CosThetaHEgen2<0.7)   {MGenCosTheta_14->Fill(rhoMC2.M());}
        // if (CosThetaHEgen2>0.7    && CosThetaHEgen2<1.)    {MGenCosTheta_15->Fill(rhoMC2.M());}
              
              
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 0.  && PhiHEgen2<2.*TMath::Pi()/12. * 1.)  {MGenPhi_0->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 1.  && PhiHEgen2<2.*TMath::Pi()/12. * 2.)  {MGenPhi_1->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 2.  && PhiHEgen2<2.*TMath::Pi()/12. * 3.)  {MGenPhi_2->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 3.  && PhiHEgen2<2.*TMath::Pi()/12. * 4.)  {MGenPhi_3->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 4.  && PhiHEgen2<2.*TMath::Pi()/12. * 5.)  {MGenPhi_4->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 5.  && PhiHEgen2<2.*TMath::Pi()/12. * 6.)  {MGenPhi_5->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 6.  && PhiHEgen2<2.*TMath::Pi()/12. * 7.)  {MGenPhi_6->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 7.  && PhiHEgen2<2.*TMath::Pi()/12. * 8.)  {MGenPhi_7->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 8.  && PhiHEgen2<2.*TMath::Pi()/12. * 9.)  {MGenPhi_8->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 9.  && PhiHEgen2<2.*TMath::Pi()/12. * 10.) {MGenPhi_9->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 10. && PhiHEgen2<2.*TMath::Pi()/12. * 11.) {MGenPhi_10->Fill(rhoMC2.M());}
        // if (PhiHEgen2>2.*TMath::Pi()/12. * 11. && PhiHEgen2<2.*TMath::Pi()/12. * 12.) {MGenPhi_11->Fill(rhoMC2.M());}






    }
  }




  TFile *SavingFile = new TFile("run2-MC/CohHe.root", "RECREATE");
  MRecH->Write();
  MGenH->Write();
  MRecMineH->Write();
  MGenMineH->Write();

  PhiRecH->Write();
  PhiRecHH->Write();
  PhiGenH->Write();
  PhiGenHH->Write();
  CosThetaRecH->Write();
  CosThetaGenH->Write();
  for (Int_t i = 0; i < 24; i++) {
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j]->Write();
      InvMassHg[i][j]->Write();
    }
  }
  FullInvMassH->Write();
  FullInvMassHgen->Write();
  PhiRecMinusGenH->Write();
  CosThetaRecMinusGenH->Write();


  MRecCosTheta_0->Write();
  MRecCosTheta_1->Write();
  MRecCosTheta_2->Write();
  MRecCosTheta_3->Write();
  MRecCosTheta_4->Write();
  MRecCosTheta_5->Write();
  MRecCosTheta_6->Write();
  MRecCosTheta_7->Write();
  MRecCosTheta_8->Write();
  MRecCosTheta_9->Write();
  MRecCosTheta_10->Write();
  MRecCosTheta_11->Write();
  MRecCosTheta_12->Write();
  MRecCosTheta_13->Write();
  MRecCosTheta_14->Write();
  MRecCosTheta_15->Write();

  MGenCosTheta_0->Write();
  MGenCosTheta_1->Write();
  MGenCosTheta_2->Write();
  MGenCosTheta_3->Write();
  MGenCosTheta_4->Write();
  MGenCosTheta_5->Write();
  MGenCosTheta_6->Write();
  MGenCosTheta_7->Write();
  MGenCosTheta_8->Write();
  MGenCosTheta_9->Write();
  MGenCosTheta_10->Write();
  MGenCosTheta_11->Write();
  MGenCosTheta_12->Write();
  MGenCosTheta_13->Write();
  MGenCosTheta_14->Write();
  MGenCosTheta_15->Write();

  MRecPhi_0->Write();
  MRecPhi_1->Write();
  MRecPhi_2->Write();
  MRecPhi_3->Write();
  MRecPhi_4->Write();
  MRecPhi_5->Write();
  MRecPhi_6->Write();
  MRecPhi_7->Write();
  MRecPhi_8->Write();
  MRecPhi_9->Write();
  MRecPhi_10->Write();
  MRecPhi_11->Write();

  MGenPhi_0->Write();
  MGenPhi_1->Write();
  MGenPhi_2->Write();
  MGenPhi_3->Write();
  MGenPhi_4->Write();
  MGenPhi_5->Write();
  MGenPhi_6->Write();
  MGenPhi_7->Write();
  MGenPhi_8->Write();
  MGenPhi_9->Write();
  MGenPhi_10->Write();
  MGenPhi_11->Write();


  SavingFile->Close();

  new TCanvas;
  PhiRecMinusGenH->Draw();










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
void DrawExample(){
  new TCanvas;
  new TCanvas;
  TFile* f = new TFile("run2-MC/CohHe.root");
  TH1F*  h = (TH1F*)f->Get("MRecH");
  BeautifyPad();
  BeautifyHisto(h);
  h->GetXaxis()->SetTitle("M_{#pi#pi} [GeV/#it{c}^{2}]");
  h->GetYaxis()->SetTitle("Reconstructed MC counts");
  h->GetYaxis()->SetRangeUser(0., h->GetMaximum()*1.5);
  h->GetXaxis()->SetRangeUser(0.3, 1.2);
  h->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"STARlight, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, reconstructed events from Run 2" );
  // latex10->DrawLatex(0.2,0.72, "0 < cos#theta < 0.08#bar{3}, #pi < #varphi < 4#pi/3");
  gPad->SaveAs("run2-MC/reconstructed-mass-MC.pdf", "recreate");
}
//______________________________________________
void DrawResolutionPhi(){
  new TCanvas;
  new TCanvas;
  TFile* f = new TFile("run2-MC/CohHe.root");
  TH1F*  h = (TH1F*)f->Get("PhiRecMinusGenH");
  BeautifyPad();
  BeautifyHisto(h);
  h->GetXaxis()->SetTitle("#varphi_{rec} - #varphi_{gen}");
  h->GetYaxis()->SetTitle("Counts [a.u.]");
  h->GetYaxis()->SetRangeUser(0., h->GetMaximum()*1.5);
  h->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
  h->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"STARlight, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, events from Run 2" );
  gPad->SaveAs("run2-MC/resolution-phi.pdf", "recreate");
}
//______________________________________________
void DrawResolutionCosTheta(){
  new TCanvas;
  new TCanvas;
  TFile* f = new TFile("run2-MC/CohHe.root");
  TH1F*  h = (TH1F*)f->Get("CosThetaRecMinusGenH");
  BeautifyPad();
  BeautifyHisto(h);
  h->GetXaxis()->SetTitle("cos#theta_{rec} - cos#theta_{gen}");
  h->GetYaxis()->SetTitle("Counts [a.u.]");
  h->GetYaxis()->SetRangeUser(0., h->GetMaximum()*1.5);
  h->GetXaxis()->SetRangeUser(-1.,1.);
  h->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"STARlight, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, events from Run 2" );
  gPad->SaveAs("run2-MC/resolution-costheta.pdf", "recreate");
}
//______________________________________________
void DrawAxE(){
  new TCanvas;
  new TCanvas;
  TFile* f = new TFile("run2-MC/CohHe.root");
  TH1F*  h  = (TH1F*)f->Get("CosThetaRecH");
  TH1F*  h2 = (TH1F*)f->Get("CosThetaGenH");
  BeautifyPad();
  BeautifyHisto(h);
  h->GetXaxis()->SetTitle("cos#theta_{rec}");
  h->GetYaxis()->SetTitle("Counts [a.u.]");
  h->GetYaxis()->SetRangeUser(0., h->GetMaximum()*1.5);
  h->GetXaxis()->SetRangeUser(-1, 1);
  h->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"STARlight, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, reconstructed events from Run 2" );
  gPad->SaveAs("run2-MC/costheta.pdf", "recreate");




  // Double_t Numerator   = 0;
  // Double_t Denominator = 0;
  // Int_t startbin = h->GetXaxis()->FindBin(-0.5);
  // Int_t endbin   = h->GetXaxis()->FindBin(0.5);
  // for (size_t i = startbin; i < endbin; i++) {
  //   Numerator   += h ->GetBinContent(i);
  //   Denominator += h2->GetBinContent(i);
  //   cout << "Rec i " << h->GetBinContent(i) << endl;
  //   // cout << "Gen i " << h2->GetBinContent(i) << endl;
  // }
  // cout << "Numerator   = " << Numerator   << endl;
  // cout << "Denominator = " << Denominator << endl;
  // cout << "AxE         = " << Numerator/Denominator << endl;
}
//______________________________________________
void DrawAxEmass(){
  new TCanvas;
  new TCanvas;
  TFile* f = new TFile("run2-MC/CohHe.root");
  TH1F*  h  = (TH1F*)f->Get("MRecH");
  TH1F*  h2 = (TH1F*)f->Get("MGenH");
  h->Sumw2();
  h2->Sumw2();
  h->Divide(h2);
  BeautifyPad();
  BeautifyHisto(h);
  h->GetXaxis()->SetTitle("M_{#pi#pi} [GeV/#it{c}^{2}]");
  h->GetYaxis()->SetTitle("AxE");
  h->GetYaxis()->SetRangeUser(0., h->GetMaximum()*1.5);
  h->GetXaxis()->SetRangeUser(0, 2);
  h->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"STARlight, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, reconstructed events from Run 2" );
  gPad->SaveAs("run2-MC/mass-axe.pdf", "recreate");




  // Double_t Numerator   = 0;
  // Double_t Denominator = 0;
  // Int_t startbin = h->GetXaxis()->FindBin(-0.5);
  // Int_t endbin   = h->GetXaxis()->FindBin(0.5);
  // for (size_t i = startbin; i < endbin; i++) {
  //   Numerator   += h ->GetBinContent(i);
  //   Denominator += h2->GetBinContent(i);
  //   cout << "Rec i " << h->GetBinContent(i) << endl;
  //   // cout << "Gen i " << h2->GetBinContent(i) << endl;
  // }
  // cout << "Numerator   = " << Numerator   << endl;
  // cout << "Denominator = " << Denominator << endl;
  // cout << "AxE         = " << Numerator/Denominator << endl;
}

