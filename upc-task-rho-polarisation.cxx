// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
/*
I ran this code using following:
o2-analysis-event-selection|o2-analysis-timestamp| o2-analysis-upc-forward --aod-file <path to ao2d.txt> [--isPbPb] -b
for now AO2D.root I am using is
alien:///alice/data/2015/LHC15o/000246392/pass5_lowIR/PWGZZ/Run3_Conversion/148_20210304-0829_child_1/AOD/001/AO2D.root
*/
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "iostream"
#include "PWGUD/DataModel/UDTables.h"
//#include "Common/DataModel/PIDResponse.h"
//#include "PWGUD/Core/RLhelper.h"
#include <TString.h>
#include "TLorentzVector.h"
using namespace std;
using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
#define mpion 0.1396 // mass of pion
//#define mmuon 0.1057 // mass of muon

struct UPCRho0Pol {
  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
    
//_____________________________________________________________________________
Double_t CosThetaHelicityFrame(TLorentzVector pionPositive,
                                TLorentzVector pionNegative,
                                    TLorentzVector possibleRhoZero)
{
/* - This function computes the Helicity cos(theta) for the

         - helicity of the RhoZero.

         - The idea should be to get back to a reference frame where it

         - is easier to compute and to define the proper z-axis.
       */

      /* - Half of the energy per pair of the colliding nucleons.*/

      Double_t HalfSqrtSnn   = 2680.;

      Double_t MassOfLead208 = 193.6823;

      Double_t MomentumBeam  = TMath::Sqrt(HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208);

      /* - Fill the Lorentz vector for projectile and target.

         - For the moment we do not consider the crossing angle.

         - Projectile runs towards the MUON arm.
       */

      TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile

      TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target

      /* - Translate the dimuon parameters in the dimuon rest frame

         -

       */

      TVector3  beta = ( -1./possibleRhoZero.E() ) * possibleRhoZero.Vect();

      TLorentzVector pPi1Dipion  = pionPositive;

      TLorentzVector pPi2Dipion  = pionNegative;

      TLorentzVector pProjDipion = pProjCM;

      TLorentzVector pTargDipion = pTargCM;

      pPi1Dipion.Boost(beta);

      pPi2Dipion.Boost(beta);

      pProjDipion.Boost(beta);

      pTargDipion.Boost(beta);

      // --- Determine the z axis for the calculation of the polarization angle

      // (i.e. the direction of the dipion in the CM system)

      TVector3 zaxis = (possibleRhoZero.Vect()).Unit();

      /* - Determine the He angle (angle between pi+ and the z axis defined above)

         -

       */

      Double_t CosThetaHE = zaxis.Dot((pPi1Dipion.Vect()).Unit());

      return   CosThetaHE;

     }
//------------------------------------------------------------------------------------------------------
   
Double_t PhiHelicityFrame(TLorentzVector muonPositive, TLorentzVector muonNegative, TLorentzVector possibleJPsi)
    {

          //Half of the energy per pair of the colliding nucleons.
          Double_t HalfSqrtSnn   = 2680.;
          Double_t MassOfLead208 = 193.6823;
          Double_t MomentumBeam  = TMath::Sqrt(HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208);

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
//-----------------------------------------------------------------------------------------------------------------------
  
    void init(o2::framework::InitContext&)
  {

    auto hSelectionCounter = registry.add<TH1>("hSelectionCounter", "hSelectionCounter;;NEvents", HistType::kTH1I, {{10, 0., 10.}});

//THIS IS THE SELECTION CUTS FOR NOW. We will add further cuts
//TString SelectionCuts[8] = {"NoSelection", "CMup11and10Trigger", "V0Selection", "FDSelection", "twotracks", "oppositecharge", "-2.5<Eta<-4", "Pt<1"};
    TString SelectionCuts[7] = {"NoSelection", "2pioncandidates","tpcNClsCrossedRows","opposite_charge","|nsigmapi|<3", "nsigmael","track_momenta>0.3GeV/c"};
//TString SelectionCuts[7] = {"NoSelection", "2pioncandidates","eta_extended","opposite_charge","|nsigmapi|<3", "nsigmael","track_momenta>0.3GeV/c"};
//TString SelectionCuts[6] = {"NoSelection", "2pioncandidates","opposite_charge","|nsigmapi|<3", "nsigmael<-7","track_momenta>0.3GeV/c"};
// now we can set BinLabel in histogram Registry

    for (int i = 0; i < 7; i++) {
      hSelectionCounter->GetXaxis()->SetBinLabel(i + 1, SelectionCuts[i].Data());
    }

    registry.add("hNsigEvsPi1", "NSigmaPi(t1) vs NSigmapi (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.},{100,-15.,15}});
    registry.add("hNsigEvsPi2", "NSigmaEl(t1) vs NSigmaEl (t2);n#sigma_{1};n#sigma_{2}", kTH2F, {{100, -15., 15.},{100,-15.,15}});
    registry.add("hNsigPivsPt1", "Pt vs NSigmaPi (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5},{100,-15.,15}});
    registry.add("hNsigPivsPt2", "Pt vs NSigmaPi (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5},{100,-15.,15}});
      
    registry.add("hNsigElvsPt1", "Pt vs NSigmaEl (t1);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5},{100,-15.,15}});
    registry.add("hNsigElvsPt2", "Pt vs NSigmaEl (t2);#it{p_{t}}, GeV/c;n#sigma_{#e}", kTH2F, {{100, 0., 2.5},{100,-15.,15}});
      
    registry.add("hNsigMuvsPt1", "Pt vs NSigmaMu (t1);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5},{100,-15.,15}});
    registry.add("hNsigMuvsPt2", "Pt vs NSigmaMu (t2);#it{p_{t}}, GeV/c;n#sigma_{#pi}", kTH2F, {{100, 0., 2.5},{100,-15.,15}});
      
    registry.add("hCharge", "Charge;#it{charge};", kTH1F, {{500, -10., 10.}});
   
    registry.add("hPtsingle_track1", "Pt t1;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});
    registry.add("hPtsingle_track2", "Pt t2;#it{p_{t}}, GeV/c;", kTH1F, {{600, 0., 3.}});
   
    registry.add("hP1", "P vs TPC signal;#it{P_{track}}, GeV/c; signal_{TPC} t1", kTH2F, {{100, 0., 2.},{300,0,150}});
    registry.add("hTPCsig", "TPC signal;signal_{TPC} t2; signal_{TPC} t2", kTH2F, {{300, 0., 150.},{300,0,150}});
    registry.add("hMPt1", "Inv.M vs track Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 2.},{100, 0., 2.}});
    
    registry.add("hPhi1", "Phi (t1);#it{#phi};", kTH1F, {{120, 0., 6.28}});
    registry.add("hPhi2", "Phi (t1);#it{#phi};", kTH1F, {{120, 0., 6.28}});
   
    registry.add("hEta_t1", "Eta of t1;#it{#eta};", kTH1F, {{100, -5., 5.}});
    registry.add("hEta_t2", "Eta of t2;#it{#eta};", kTH1F, {{100, -5., 5.}});
      
//        registry.add("hCostheta_Phi", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0., 6.28},{120,-1,1}});
      registry.add("hCostheta_Phi", "Phi vs Costheta;#it{#phi};#it{Cos#Theta};", kTH2F, {{100, 0.*TMath::Pi(), 2.*TMath::Pi()},{120,-1,1}});

      registry.add("hMass", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
      registry.add("hMass1", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
      registry.add("hMass2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
      
        registry.add("hMassCosTheta_0", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_1", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_3", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_4", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_5", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_6", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_7", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_8", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_9", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_10", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_11", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_12", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_13", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_14", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_15", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassCosTheta_16", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        
        
        registry.add("hMassPhi_0", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_1", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_3", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_4", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_5", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_6", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_7", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_8", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_9", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_10", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_11", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        registry.add("hMassPhi_12", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
        
        
      registry.add("hMasseta2", "Raw Inv.M;#it{m_{#mu#mu}}, GeV/c^{2};", kTH1D, {{100, 0., 2.}});
      registry.add("hPt", "Pt;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
      registry.add("hPt1", "Pt1;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
      registry.add("hPt2", "Pt2;#it{p_{t}}, GeV/c;", kTH1D, {{500, 0., 5.}});
      registry.add("hEta", "Eta;#it{#eta};", kTH1F, {{500, -10., 10.}});
      registry.add("hRap", "Rapidity;#it{y};", kTH1F, {{500, -10., 10.}});
      registry.add("hPhi", "Phi;#it{#Phi};", kTH1F, {{120, 0, 2.*TMath::Pi()}});
      registry.add("hCostheta", "Costheta;#it{Cos#Theta};", kTH1F, {{100, -1., 1.}});
      registry.add("hMPt", "Inv.M vs Pt;M, GeV/c^{2};#it{P_{t}}, GeV/c;", kTH2F, {{100, 0., 2.},{100, 0., 2.}});
  }
  
  //Configurable<float> etalow{"etalow", -1.0f, ""}; //
  //Configurable<float> etahigh{"etahigh", 1.0f, ""}; //
  
  using udtracks = soa::Join<aod::UDTracks,aod::UDTracksExtra,aod::UDTracksPID>;
  //using udtracks = soa::Join<aod::UDTracks,aod::UDTracksExtra,aod::UDTracksPID,aod::UDTracksDCA>;

  void process(UDCollisions::iterator const& collision, udtracks const& tracks)
  {
      registry.fill(HIST("hSelectionCounter"), 0);

      TLorentzVector t1, t2, t3, t4, p; // lorentz vectors of tracks and the mother
      //cout <<"###################################################################there are two tracks events"<< tracks.size()<<endl;
      //if(tracks.tpcNClsCrossedRows()<70) continue; //skip bad tracks from TPC
      //if(fabs(track.dcaXY()) > 0.2) continue; //skip tracks that does not point to primary vertex
      
      if (tracks.size()!=2)return;

      registry.fill(HIST("hSelectionCounter"), 1);

      //cout << "number of track is " << tracks.size()<< endl;

      auto track1 = tracks.begin();
      auto track2 = track1+1;
    
      //if(track1.tpcNClsCrossedRows()<70) return;
      //if(track2.tpcNClsCrossedRows()<70) return;
      
      t1.SetXYZM(track1.px(), track1.py(), track1.pz(),mpion);
      t2.SetXYZM(track2.px(), track2.py(), track2.pz(),mpion);
      p=t1+t2;
      
    
      registry.fill(HIST("hSelectionCounter"), 2);
      
     // registry.fill(HIST("hSelectionCounter"), 3);
      
          
      float nsigmapit1 = track1.tpcNSigmaPi();
      float nsigmapit2 = track2.tpcNSigmaPi();
      
      if(track1.sign()==track2.sign()) return;
      
      TLorentzVector tplus,tminus;
      if(track1.sign()>0)tplus=t1;
      else tminus = t1;
      if(track2.sign()<0)tminus=t2;
      else tplus = t2;
      
      registry.fill(HIST("hSelectionCounter"), 3);
      //if(TMath::Abs(nsigmapit1)>2 || TMath::Abs(nsigmapit2) >2) return;
      if((nsigmapit1 * nsigmapit1 + nsigmapit2 * nsigmapit2) > 25.0) return;
     // if(nsigmapit1<-3 || nsigmapit1>3)  return;
     // if(nsigmapit2<-3 || nsigmapit2>3)  return;
      
      registry.fill(HIST("hSelectionCounter"), 4);
      if(track1.tpcNSigmaEl()>-7.0 || track2.tpcNSigmaEl()>-7.0)  return;
      
      
      registry.fill(HIST("hSelectionCounter"), 5);
      if(t1.P()<0.3 || t2.P()<0.3) return;
      //if(p.Pt()>0.2) return;
            
      registry.fill(HIST("hSelectionCounter"), 6);
      //cout<<"pidt1==="<<nsigmapit1<<"    "<<"pidt2====="<<nsigmapit2<<endl;
      if(p.Rapidity()<-0.8 || p.Rapidity() > 0.8) return;
      
      
      registry.fill(HIST("hNsigEvsPi1"), track1.tpcNSigmaPi(),track2.tpcNSigmaPi());
      registry.fill(HIST("hNsigEvsPi2"), track1.tpcNSigmaEl(),track2.tpcNSigmaEl());
      registry.fill(HIST("hNsigPivsPt1"), t1.Pt(),track1.tpcNSigmaPi());
      registry.fill(HIST("hNsigPivsPt2"), t2.Pt(),track2.tpcNSigmaPi());
      registry.fill(HIST("hNsigElvsPt1"), t1.Pt(),track1.tpcNSigmaEl());
      registry.fill(HIST("hNsigElvsPt2"), t2.Pt(),track2.tpcNSigmaEl());
      registry.fill(HIST("hNsigMuvsPt1"), t1.Pt(),track1.tpcNSigmaMu());
      registry.fill(HIST("hNsigMuvsPt2"), t2.Pt(),track2.tpcNSigmaMu());

      registry.fill(HIST("hCharge"),  track1.sign());
      registry.fill(HIST("hPtsingle_track1"), t1.Pt());
      registry.fill(HIST("hPtsingle_track2"), t2.Pt());
      registry.fill(HIST("hP1"), t1.P(),track1.tpcSignal());
      registry.fill(HIST("hTPCsig"),track1.tpcSignal(),track2.tpcSignal());
      registry.fill(HIST("hMPt1"), p.M(),t1.Pt());
      
      
      //if(p.M()<0.85 && p.M()>0.65 && p.Pt() <0.2) {
          
          if(p.Pt() <0.2) {
          
              auto costheta = CosThetaHelicityFrame(tplus,tminus,p);
              registry.fill(HIST("hCostheta"), costheta);
              auto phihel = 1.*TMath::Pi() + PhiHelicityFrame(tplus,tminus,p);
              registry.fill(HIST("hPhi"), phihel);
              registry.fill(HIST("hCostheta_Phi"), phihel,costheta);
              
              if (costheta>-1.    && costheta<-0.7)  {registry.fill(HIST("hMassCosTheta_0"), p.M());};
              if (costheta>-0.7   && costheta<-0.55) {registry.fill(HIST("hMassCosTheta_1"), p.M());};
              if (costheta>-0.55  && costheta<-0.45) {registry.fill(HIST("hMassCosTheta_2"), p.M());};
              if (costheta>-0.45  && costheta<-0.35) {registry.fill(HIST("hMassCosTheta_3"), p.M());};
              if (costheta>-0.35  && costheta<-0.25) {registry.fill(HIST("hMassCosTheta_4"), p.M());};
              if (costheta>-0.25  && costheta<-0.15) {registry.fill(HIST("hMassCosTheta_5"), p.M());};
              if (costheta>-0.15  && costheta<-0.05) {registry.fill(HIST("hMassCosTheta_6"), p.M());};
              if (costheta>-0.05  && costheta<0.)    {registry.fill(HIST("hMassCosTheta_7"), p.M());};
              if (costheta>0.     && costheta<0.05)  {registry.fill(HIST("hMassCosTheta_8"), p.M());};
              if (costheta>0.05   && costheta<0.15)  {registry.fill(HIST("hMassCosTheta_9"), p.M());};
              if (costheta>0.15   && costheta<0.25)  {registry.fill(HIST("hMassCosTheta_10"), p.M());};
              if (costheta>0.25   && costheta<0.35)  {registry.fill(HIST("hMassCosTheta_11"), p.M());};
              if (costheta>0.35   && costheta<0.45)  {registry.fill(HIST("hMassCosTheta_12"), p.M());};
              if (costheta>0.45   && costheta<0.55)  {registry.fill(HIST("hMassCosTheta_13"), p.M());};
              if (costheta>0.55   && costheta<0.7)   {registry.fill(HIST("hMassCosTheta_14"), p.M());};
              if (costheta>0.7    && costheta<1.)    {registry.fill(HIST("hMassCosTheta_15"), p.M());};
              
              
              if (phihel>2.*TMath::Pi()/12. * 0.  && phihel<2.*TMath::Pi()/12. * 1.)  {registry.fill(HIST("hMassPhi_0"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 1.  && phihel<2.*TMath::Pi()/12. * 2.)  {registry.fill(HIST("hMassPhi_1"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 2.  && phihel<2.*TMath::Pi()/12. * 3.)  {registry.fill(HIST("hMassPhi_2"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 3.  && phihel<2.*TMath::Pi()/12. * 4.)  {registry.fill(HIST("hMassPhi_3"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 4.  && phihel<2.*TMath::Pi()/12. * 5.)  {registry.fill(HIST("hMassPhi_4"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 5.  && phihel<2.*TMath::Pi()/12. * 6.)  {registry.fill(HIST("hMassPhi_5"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 6.  && phihel<2.*TMath::Pi()/12. * 7.)  {registry.fill(HIST("hMassPhi_6"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 7.  && phihel<2.*TMath::Pi()/12. * 8.)  {registry.fill(HIST("hMassPhi_7"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 8.  && phihel<2.*TMath::Pi()/12. * 9.)  {registry.fill(HIST("hMassPhi_8"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 9.  && phihel<2.*TMath::Pi()/12. * 10.) {registry.fill(HIST("hMassPhi_9"),  p.M());};
              if (phihel>2.*TMath::Pi()/12. * 10. && phihel<2.*TMath::Pi()/12. * 11.) {registry.fill(HIST("hMassPhi_10"), p.M());};
              if (phihel>2.*TMath::Pi()/12. * 11. && phihel<2.*TMath::Pi()/12. * 12.) {registry.fill(HIST("hMassPhi_11"), p.M());};

      }
      
        
        
     // if((t1.Eta()>1. && t1.Eta()<1.5) && (t2.Eta()>1. && t2.Eta()<1.5)) {registry.fill(HIST("hMasseta1"), p.M());};
      
      
      
      
      
      if(p.Pt() <0.2) {registry.fill(HIST("hMass1"),   p.M());};
      if(p.Pt() >0.2) {registry.fill(HIST("hMass2"),   p.M());};
            
      
      if(track1.sign()>0&&track2.sign()>0) {registry.fill(HIST("hPt1"),   p.Pt());};
      if(track1.sign()<0&&track2.sign()<0) {registry.fill(HIST("hPt2"),   p.Pt());};
     
      registry.fill(HIST("hPt"),   p.Pt());
      registry.fill(HIST("hMass"), p.M());
      registry.fill(HIST("hEta"),  p.Eta());
      registry.fill(HIST("hRap"),  p.Rapidity());
      //registry.fill(HIST("hPhi"),  p.Phi());
      registry.fill(HIST("hMPt"),  p.M(),p.Pt());
      
      //cout << " Pt of di pion is "<< p.Pt() << endl; // pion hypothesis no pid selection for now
      //if (ntracks!=2) return; // selecting two tracks event */

  } // end of process

}; // end of struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UPCRho0Pol>(cfgc)};
}


