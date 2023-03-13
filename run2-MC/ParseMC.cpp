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
  Float_t Mrec;
  Float_t Mgen;
  ReconTree->SetBranchAddress("Mass_T",  &Mrec);
  GenerTree->SetBranchAddress("Mass_MC_T",  &Mgen);
  Float_t Yrec;
  Float_t Ygen;
  ReconTree->SetBranchAddress("Rapidity_T",  &Yrec);
  Short_t  Qrec[2];
  ReconTree->SetBranchAddress("TrackQ_T",  &Qrec);
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
  ReconTree->SetBranchAddress("TrackPtGen_T",&SinglePtgen2);
  GenerTree->SetBranchAddress("Pt1_MC_T",&SinglePtgen[0]);
  GenerTree->SetBranchAddress("Pt2_MC_T",&SinglePtgen[1]);;
  Float_t SingleEtarec[2];
  Float_t SingleEtagen[2];
  Float_t SingleEtagen2[2];
  ReconTree->SetBranchAddress("TrackEta_T",   &SingleEtarec);
  ReconTree->SetBranchAddress("TrackEtaGen_T",&SingleEtagen2);
  GenerTree->SetBranchAddress("Eta1_MC_T",&SingleEtagen[0]);
  GenerTree->SetBranchAddress("Eta2_MC_T",&SingleEtagen[1]);
  Float_t SinglePhirec[2];
  Float_t SinglePhigen[2];
  Float_t SinglePhigen2[2];
  ReconTree->SetBranchAddress("TrackPhi_T",   &SinglePhirec);
  ReconTree->SetBranchAddress("TrackPhiGen_T",&SinglePhigen2);
  GenerTree->SetBranchAddress("Phi1_MC_T",&SinglePhigen[0]);
  GenerTree->SetBranchAddress("Phi2_MC_T",&SinglePhigen[1]);

  Double_t CosThetaHErec;
  Double_t CosThetaHEgen;
  Double_t PhiHErec;
  Double_t PhiHEgen;

  Int_t nentriesRec = (Int_t) ReconTree->GetEntries();
  // Int_t nentriesGen = (Int_t) GenerTree->GetEntries();


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

  TH1F *PhiRecMinusGenH      = new TH1F("PhiRecMinusGenH",      "PhiRecMinusGenH",      50, -2.*TMath::Pi(), 2.*TMath::Pi());
  TH1F *CosThetaRecMinusGenH = new TH1F("CosThetaRecMinusGenH", "CosThetaRecMinusGenH", 50, -2., 2.);


  TH1F *CosThetaRecH = new TH1F("CosThetaRecH", "CosThetaRecH", 250, -1., 1.);
  TH1F *CosThetaGenH = new TH1F("CosThetaGenH", "CosThetaGenH", 250, -1., 1.);
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      250, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      250, 0., 2.*TMath::Pi());
  TH1F *PhiRecHH     = new TH1F("PhiRecHH",     "PhiRecHH",     500, -2.*TMath::Pi(), 2.*TMath::Pi());
  TH1F *PhiGenHH     = new TH1F("PhiGenHH",     "PhiGenHH",     500, -2.*TMath::Pi(), 2.*TMath::Pi());



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




  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);
    TLorentzVector piPlus, piMinus, rho;
    TLorentzVector piPlusMC, piMinusMC, rhoMC;
    TLorentzVector piPlusMC2, piMinusMC2, rhoMC2;
    SinglePtrec[0] = TMath::Sqrt(SinglePXrec[0]*SinglePXrec[0] + SinglePYrec[0]*SinglePYrec[0]);
    SinglePtrec[1] = TMath::Sqrt(SinglePXrec[1]*SinglePXrec[1] + SinglePYrec[1]*SinglePYrec[1]);
    Double_t EnergyRec[2], EnergyGen[2], EnergyGen2[2];
    EnergyRec[0] = TMath::Sqrt( SinglePXrec[0]*SinglePXrec[0] + SinglePYrec[0]*SinglePYrec[0] + SinglePZrec[0]*SinglePZrec[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyRec[1] = TMath::Sqrt( SinglePXrec[1]*SinglePXrec[1] + SinglePYrec[1]*SinglePYrec[1] + SinglePZrec[1]*SinglePZrec[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    Float_t pXgen[2], pYgen[2], pZgen[2];
    pXgen[0] = SinglePtgen[0] * TMath::Cos(SinglePhigen[0]);
    pXgen[1] = SinglePtgen[1] * TMath::Cos(SinglePhigen[1]);
    pYgen[0] = SinglePtgen[0] * TMath::Sin(SinglePhigen[0]);
    pYgen[1] = SinglePtgen[1] * TMath::Sin(SinglePhigen[1]);
    pZgen[0] = SinglePtgen[0] * TMath::SinH(SingleEtagen[0]);
    pZgen[1] = SinglePtgen[1] * TMath::SinH(SingleEtagen[1]);
    Float_t pXgen2[2], pYgen2[2], pZgen2[2];
    // SinglePtgen2[0] -= 0.0435;
    // SinglePtgen2[1] -= 0.0435;
    pXgen2[0] = SinglePtgen2[0] * TMath::Cos(SinglePhigen2[0]);
    pXgen2[1] = SinglePtgen2[1] * TMath::Cos(SinglePhigen2[1]);
    pYgen2[0] = SinglePtgen2[0] * TMath::Sin(SinglePhigen2[0]);
    pYgen2[1] = SinglePtgen2[1] * TMath::Sin(SinglePhigen2[1]);
    pZgen2[0] = SinglePtgen2[0] * TMath::SinH(SingleEtagen2[0]);
    pZgen2[1] = SinglePtgen2[1] * TMath::SinH(SingleEtagen2[1]);
    // EnergyGen[0] = TMath::Sqrt( SinglePtgen[0]*SinglePtgen[0] + pZgen[0]*pZgen[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    // EnergyGen[1] = TMath::Sqrt( SinglePtgen[1]*SinglePtgen[1] + pZgen[1]*pZgen[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyGen[0] = TMath::Sqrt( pXgen[0]*pXgen[0] + pYgen[0]*pYgen[0] + pZgen[0]*pZgen[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyGen[1] = TMath::Sqrt( pXgen[1]*pXgen[1] + pYgen[1]*pYgen[1] + pZgen[1]*pZgen[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyGen2[0] = TMath::Sqrt( pXgen2[0]*pXgen2[0] + pYgen2[0]*pYgen2[0] + pZgen2[0]*pZgen2[0] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    EnergyGen2[1] = TMath::Sqrt( pXgen2[1]*pXgen2[1] + pYgen2[1]*pYgen2[1] + pZgen2[1]*pZgen2[1] + TDatabasePDG::Instance()->GetParticle(211)->Mass()*TDatabasePDG::Instance()->GetParticle(211)->Mass());
    if(Qrec[0] < -0.5 && Qrec[1] > 0.5){
      piPlus.SetPxPyPzE(    SinglePXrec[1],SinglePYrec[1],SinglePZrec[1],EnergyRec[1]);
      piMinus.SetPxPyPzE(   SinglePXrec[0],SinglePYrec[0],SinglePZrec[0],EnergyRec[0]);
      piPlusMC.SetPxPyPzE(  pXgen[1],pYgen[1],pZgen[1],EnergyGen[1]);
      piMinusMC.SetPxPyPzE( pXgen[0],pYgen[0],pZgen[0],EnergyGen[0]);
      piPlusMC2.SetPxPyPzE(  pXgen2[1],pYgen2[1],pZgen2[1],EnergyGen2[1]);
      piMinusMC2.SetPxPyPzE( pXgen2[0],pYgen2[0],pZgen2[0],EnergyGen2[0]);
      // piPlus.SetPtEtaPhiM( SinglePtrec[1],SingleEtarec[1],SinglePhirec[1],TDatabasePDG::Instance()->GetParticle(211)->Mass());
      // piMinus.SetPtEtaPhiM(SinglePtrec[0],SingleEtarec[0],SinglePhirec[0],TDatabasePDG::Instance()->GetParticle(211)->Mass());
      // piPlusMC.SetPtEtaPhiM( SinglePtgen[1],SingleEtagen[1],SinglePhigen[1],TDatabasePDG::Instance()->GetParticle(211)->Mass());
      // piMinusMC.SetPtEtaPhiM(SinglePtgen[0],SingleEtagen[0],SinglePhigen[0],TDatabasePDG::Instance()->GetParticle(211)->Mass());
    } else if(Qrec[1] < -0.5 && Qrec[0] > 0.5){
      // piPlus.SetPtEtaPhiM( SinglePtrec[0],SingleEtarec[0],SinglePhirec[0],TDatabasePDG::Instance()->GetParticle(211)->Mass());
      // piMinus.SetPtEtaPhiM(SinglePtrec[1],SingleEtarec[1],SinglePhirec[1],TDatabasePDG::Instance()->GetParticle(211)->Mass());
      piPlus.SetPxPyPzE(  SinglePXrec[0],SinglePYrec[0],SinglePZrec[0],EnergyRec[0]);
      piMinus.SetPxPyPzE( SinglePXrec[1],SinglePYrec[1],SinglePZrec[1],EnergyRec[1]);
      piPlusMC.SetPxPyPzE(  pXgen[0],pYgen[0],pZgen[0],EnergyGen[0]);
      piMinusMC.SetPxPyPzE( pXgen[1],pYgen[1],pZgen[1],EnergyGen[1]);
      piPlusMC2.SetPxPyPzE(  pXgen2[0],pYgen2[0],pZgen2[0],EnergyGen2[0]);
      piMinusMC2.SetPxPyPzE( pXgen2[1],pYgen2[1],pZgen2[1],EnergyGen2[1]);
      // piPlusMC.SetPtEtaPhiM( SinglePtgen[0],SingleEtagen[0],SinglePhigen[0],TDatabasePDG::Instance()->GetParticle(211)->Mass());
      // piMinusMC.SetPtEtaPhiM(SinglePtgen[1],SingleEtagen[1],SinglePhigen[1],TDatabasePDG::Instance()->GetParticle(211)->Mass());
    }
    rho   += piPlus;
    rho   += piMinus;
    rhoMC += piPlusMC;
    rhoMC += piMinusMC;
    rhoMC2 += piPlusMC2;
    rhoMC2 += piMinusMC2;


      if( (rhoMC.Pt() > 0. && rhoMC.Pt() < 0.20) && (rhoMC.Y() > -0.8 && rhoMC.Y() < 0.8) ) {
        // MGenMineH->Fill(rhoMC.M());
      }  

      if( (pTrec > 0. && pTrec < 0.20) && (Yrec > -0.8 && Yrec < 0.8) ) {
      // if( true ) {

        MRecH->Fill(Mrec);
        MRecMineH->Fill(rho.M());
        MGenMineH->Fill(rhoMC.M());
        MGenMine2H->Fill(rhoMC2.M());
        MGenH->Fill(Mgen);
        PtGenH->Fill(rhoMC.Pt());
        PtGen2H->Fill(rhoMC2.Pt());
        PtGenS2H->Fill(SinglePtgen2[0]);
        PtGenS2H->Fill(SinglePtgen2[1]);
        PtGenSH->Fill(SinglePtgen[0]);
        PtGenSH->Fill(SinglePtgen[1]);
        EtaGenSH->Fill(SingleEtagen2[0]);
        EtaGenSH->Fill(SingleEtagen2[1]);
        PhiGenSH->Fill(SinglePhigen2[0]);
        PhiGenSH->Fill(SinglePhigen2[1]);



        // CosThetaHEgen2 = -1.*CosThetaHEgen;
        // if((CosThetaHEgen2 < 1.0) && (CosThetaHEgen2 > -1.0)){
        // PhiHEgen2 = PhiHEgen+TMath::Pi();
        // PhiGenH->Fill( (PhiHEgen+TMath::Pi()) );
        // if((CosThetaHEgen2 < 0.5) && (CosThetaHEgen2 > -0.5)){
        //   PhiGenHH->Fill( (PhiHEgen+TMath::Pi()) );
        // }
        // CosThetaGenH->Fill( CosThetaHEgen2 );
        // //================================
        // // Translating Phi and CosTheta
        // // into bin numbers
        // //--------------------------------
        // Double_t TraslatedCosThetaGeng = 0.5*(CosThetaHEgen2 + 1.)*24.;
        // Double_t iCosThetaBins2g = -1;
        // Double_t RemainderCosThetag = 0.;
        // RemainderCosThetag = modf(TraslatedCosThetaGeng, &iCosThetaBins2g);
        // Int_t iCosThetaBinsg = (Int_t)  iCosThetaBins2g;
        // //--------------------------------
        // // Binning in phi depends
        // // on CosTheta
        // //--------------------------------
        // Double_t Mg = 6.;
        // // Double_t TraslatedPhiGen = (PhiHEgen2+TMath::Pi())*M/(2.*TMath::Pi());
        // Double_t TraslatedPhiGeng = (PhiHEgen2)*Mg/(2.*TMath::Pi());
        // Double_t iPhiBins2g = -1;
        // Double_t RemainderPhig = 0.;
        // RemainderPhig = modf(TraslatedPhiGeng, &iPhiBins2g);
        // Int_t iPhiBinsg = (Int_t)  iPhiBins2g;
        // //--------------------------------------
        // InvMassHg[iCosThetaBinsg][iPhiBinsg]->Fill(Mrec);

        // if(iCosThetaBinsg > 5 && iCosThetaBinsg < 18) FullInvMassHgen->Fill(1);


          // if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {


          //         // PhiRecMinusGenH->Fill(PhiHErec-PhiHEgen2);
          //         if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
          //         if(TMath::Abs(PhiHErec - PhiHEgen2) > TMath::Pi()){
          //           PhiRecMinusGenH->Fill(PhiHErec-PhiHEgen2+TMath::Pi());
          //         }  else {
          //           PhiRecMinusGenH->Fill(PhiHErec-PhiHEgen2);
          //         }
          //         PhiRecH->Fill( PhiHErec );
          //         if((CosThetaHEgen2 < 0.5) && (CosThetaHEgen2 > -0.5)){
          //           PhiRecHH->Fill( PhiHErec );
          //         }
          //         CosThetaRecH->Fill( CosThetaHErec );

          //         //================================
          //         // Translating Phi and CosTheta
          //         // into bin numbers
          //         //--------------------------------
          //         Double_t TraslatedCosThetaGen = 0.5*(CosThetaHEgen2 + 1.)*24.;
          //         Double_t iCosThetaBins2 = -1;
          //         Double_t RemainderCosTheta = 0.;
          //         RemainderCosTheta = modf(TraslatedCosThetaGen, &iCosThetaBins2);
          //         Int_t iCosThetaBins = (Int_t)  iCosThetaBins2;
          //         //--------------------------------
          //         // Binning in phi depends
          //         // on CosTheta
          //         //--------------------------------
          //         Double_t M = 6.;
          //         // Double_t TraslatedPhiGen = (PhiHEgen2+TMath::Pi())*M/(2.*TMath::Pi());
          //         Double_t TraslatedPhiGen = (PhiHErec)*M/(2.*TMath::Pi());
          //         Double_t iPhiBins2 = -1;
          //         Double_t RemainderPhi = 0.;
          //         RemainderPhi = modf(TraslatedPhiGen, &iPhiBins2);
          //         Int_t iPhiBins = (Int_t)  iPhiBins2;
          //         //--------------------------------------
          //         InvMassH[iCosThetaBins][iPhiBins]->Fill(Mrec);




          //         if(iCosThetaBins > 5 && iCosThetaBins < 18) FullInvMassH->Fill(Mrec);


          // }
      // }
    }
  }




  TFile *SavingFile = new TFile("run2-MC/CohHe.root", "RECREATE");
  MRecH->Write();
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
  TFile* f = new TFile("rho-polarisation/run2-MC/CohHe.root");
  TH1F*  h = (TH1F*)f->Get("MrecH");
  BeautifyPad();
  BeautifyHisto(h);
  h->GetXaxis()->SetTitle("M_{#mu#mu} [GeV/#it{c}^{2}]");
  h->GetYaxis()->SetTitle("Reconstructed MC counts");
  h->GetYaxis()->SetRangeUser(0., h->GetMaximum()*1.5);
  h->GetXaxis()->SetRangeUser(2.2, 4.0);
  h->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, reconstructed events" );
  latex10->DrawLatex(0.2,0.72, "0 < cos#theta < 0.08#bar{3}, #pi < #varphi < 4#pi/3");
  gPad->SaveAs("SignalExtractionCoarse/ClosureYields/reconstructed-MC.pdf", "recreate");
}
//______________________________________________
void DrawResolution(){
  new TCanvas;
  new TCanvas;
  TFile* f = new TFile("SignalExtractionCoarse/ClosureYields/CohHe.root");
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
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, LHC18l7" );
  gPad->SaveAs("SignalExtractionCoarse/ClosureYields/resolution-phi.pdf", "recreate");
}
//______________________________________________
void DrawAxE(){
  new TCanvas;
  new TCanvas;
  TFile* f = new TFile("SignalExtractionCoarse/ClosureYields/CohHe.root");
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
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80, "ALICE Simulation, LHC18l7" );
  gPad->SaveAs("SignalExtractionCoarse/ClosureYields/costheta.pdf", "recreate");




  Double_t Numerator   = 0;
  Double_t Denominator = 0;
  Int_t startbin = h->GetXaxis()->FindBin(-0.5);
  Int_t endbin   = h->GetXaxis()->FindBin(0.5);
  for (size_t i = startbin; i < endbin; i++) {
    Numerator   += h ->GetBinContent(i);
    Denominator += h2->GetBinContent(i);
    cout << "Rec i " << h->GetBinContent(i) << endl;
    // cout << "Gen i " << h2->GetBinContent(i) << endl;
  }
  cout << "Numerator   = " << Numerator   << endl;
  cout << "Denominator = " << Denominator << endl;
  cout << "AxE         = " << Numerator/Denominator << endl;
}
