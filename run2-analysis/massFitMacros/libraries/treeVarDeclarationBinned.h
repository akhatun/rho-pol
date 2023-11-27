#ifndef TREEVARDECLARATIONBINNED_H
#define TREEVARDECLARATIONBINNED_H

    void setVal(bool settingVal ,bool &valToSet){
        valToSet = settingVal;
    }

    //TTree globl pionter
    TTree *globalTree = NULL;

    //flag for using new data
    extern bool newData;
    //falg for use deltaPhi
    extern bool useDeltaPhi;

    // Declaration of leaf types of REAL DATA EVENTS 
    Int_t           RunNum_T;
    UInt_t          PeriodNumber_T;
    UInt_t          OrbitNumber_T;
    UShort_t        BunchCrossNumber_T;
    Bool_t          LikeSign_T;
    Float_t         Mass_T;
    Float_t         Pt_T;
    Float_t         Rapidity_T;
    Float_t         Phi_T;
    Float_t         ZNAenergy_T;
    Float_t         ZNCenergy_T;
    Float_t         ZPAenergy_T;
    Float_t         ZPCenergy_T;
    Float_t         ZDCAtime_T[4];
    Float_t         ZDCCtime_T[4];
    Float_t         PIDTPCPion_T[2];
    Float_t         PIDTPCElectron_T[2];
    Int_t           TPCsignal_T[2];
    Float_t         TrackP_T[2];
    Float_t         TrackEta_T[2];
    Float_t         TrackPhi_T[2];
    Float_t         TrackPx_T[2];
    Float_t         TrackPy_T[2];
    Float_t         TrackPz_T[2];
    Short_t         TrackQ_T[2];
    Float_t         VtxX_T;
    Float_t         VtxY_T;
    Float_t         VtxZ_T;
    Int_t           VtxContrib_T;
    Float_t         VtxChi2_T;
    Float_t         VtxNDF_T;
    Float_t         SpdVtxX_T;
    Float_t         SpdVtxY_T;
    Float_t         SpdVtxZ_T;
    Int_t           SpdVtxContrib_T;
    Int_t           V0Adecision_T;
    Int_t           V0Cdecision_T;
    Int_t           ADAdecision_T;
    Int_t           ADCdecision_T;
    Bool_t          UBAfired_T;
    Bool_t          UBCfired_T;
    Bool_t          VBAfired_T;
    Bool_t          VBCfired_T;
    Int_t           Ntracklets_T;
    Bool_t          ChipCut_T;
    Bool_t          TriggerSPD_T;
    Bool_t          TriggerTOF_T;
    
    Float_t            deltaPhi; //used also in generated events

    Float_t         meanDeltaPhi;
    Float_t         meanPt;
    
    void createDataTree(TFile *globalFile, string treePath){
        globalTree = (TTree*)globalFile->Get(treePath.c_str());

        globalTree->SetBranchAddress("RunNum_T", &RunNum_T);
        globalTree->SetBranchAddress("PeriodNumber_T", &PeriodNumber_T);
        globalTree->SetBranchAddress("OrbitNumber_T", &OrbitNumber_T);
        globalTree->SetBranchAddress("BunchCrossNumber_T", &BunchCrossNumber_T);
        globalTree->SetBranchAddress("LikeSign_T", &LikeSign_T);
        globalTree->SetBranchAddress("Mass_T", &Mass_T);
        globalTree->SetBranchAddress("Pt_T", &Pt_T);
        globalTree->SetBranchAddress("Rapidity_T", &Rapidity_T);
        globalTree->SetBranchAddress("Phi_T", &Phi_T);
        globalTree->SetBranchAddress("ZNAenergy_T", &ZNAenergy_T);
        globalTree->SetBranchAddress("ZNCenergy_T", &ZNCenergy_T);
        globalTree->SetBranchAddress("ZPAenergy_T", &ZPAenergy_T);
        globalTree->SetBranchAddress("ZPCenergy_T", &ZPCenergy_T);
        globalTree->SetBranchAddress("ZDCAtime_T", ZDCAtime_T);
        globalTree->SetBranchAddress("ZDCCtime_T", ZDCCtime_T);
        globalTree->SetBranchAddress("PIDTPCPion_T", PIDTPCPion_T);
        globalTree->SetBranchAddress("PIDTPCElectron_T", PIDTPCElectron_T);
        globalTree->SetBranchAddress("TPCsignal_T", TPCsignal_T);
        globalTree->SetBranchAddress("TrackP_T", TrackP_T);
        globalTree->SetBranchAddress("TrackEta_T", TrackEta_T);
        globalTree->SetBranchAddress("TrackPhi_T", TrackPhi_T);
        globalTree->SetBranchAddress("TrackPx_T", TrackPx_T);
        globalTree->SetBranchAddress("TrackPy_T", TrackPy_T);
        globalTree->SetBranchAddress("TrackPz_T", TrackPz_T);
        if(newData) globalTree->SetBranchAddress("TrackQ_T", TrackQ_T);
        globalTree->SetBranchAddress("VtxX_T", &VtxX_T);
        globalTree->SetBranchAddress("VtxY_T", &VtxY_T);
        globalTree->SetBranchAddress("VtxZ_T", &VtxZ_T);
        globalTree->SetBranchAddress("VtxContrib_T", &VtxContrib_T);
        globalTree->SetBranchAddress("VtxChi2_T", &VtxChi2_T);
        globalTree->SetBranchAddress("VtxNDF_T", &VtxNDF_T);
        globalTree->SetBranchAddress("SpdVtxX_T", &SpdVtxX_T);
        globalTree->SetBranchAddress("SpdVtxY_T", &SpdVtxY_T);
        globalTree->SetBranchAddress("SpdVtxZ_T", &SpdVtxZ_T);
        globalTree->SetBranchAddress("SpdVtxContrib_T", &SpdVtxContrib_T);
        globalTree->SetBranchAddress("V0Adecision_T", &V0Adecision_T);
        globalTree->SetBranchAddress("V0Cdecision_T", &V0Cdecision_T);
        globalTree->SetBranchAddress("ADAdecision_T", &ADAdecision_T);
        globalTree->SetBranchAddress("ADCdecision_T", &ADCdecision_T);
        globalTree->SetBranchAddress("UBAfired_T", &UBAfired_T);
        globalTree->SetBranchAddress("UBCfired_T", &UBCfired_T);
        globalTree->SetBranchAddress("VBAfired_T", &VBAfired_T);
        globalTree->SetBranchAddress("VBCfired_T", &VBCfired_T);
        globalTree->SetBranchAddress("Ntracklets_T", &Ntracklets_T);
        globalTree->SetBranchAddress("ChipCut_T", &ChipCut_T);
        globalTree->SetBranchAddress("TriggerSPD_T", &TriggerSPD_T);
        globalTree->SetBranchAddress("TriggerTOF_T", &TriggerTOF_T);   
        
        if(useDeltaPhi) globalTree->SetBranchAddress("deltaPhi", &deltaPhi);   

        globalTree->SetBranchAddress("meanPt",&meanPt);
        globalTree->SetBranchAddress("meanDeltaPhi",&meanDeltaPhi);
    }
    
    // Declaration of leaf types of RECONSTRUCTED EVENTS
    Float_t         TrackEtaGen_T[2];
    Float_t         TrackPhiGen_T[2];
    Float_t         TrackPtGen_T[2];

    void createMCRecoTree(TFile *globalFile, string treePath){
        globalTree = (TTree*)globalFile->Get(treePath.c_str());

        globalTree->SetBranchAddress("RunNum_T", &RunNum_T);
        globalTree->SetBranchAddress("PeriodNumber_T", &PeriodNumber_T);
        globalTree->SetBranchAddress("OrbitNumber_T", &OrbitNumber_T);
        globalTree->SetBranchAddress("BunchCrossNumber_T", &BunchCrossNumber_T);
        globalTree->SetBranchAddress("LikeSign_T", &LikeSign_T);
        globalTree->SetBranchAddress("Mass_T", &Mass_T);
        globalTree->SetBranchAddress("Pt_T", &Pt_T);
        globalTree->SetBranchAddress("Rapidity_T", &Rapidity_T);
        globalTree->SetBranchAddress("Phi_T", &Phi_T);
        globalTree->SetBranchAddress("ZNAenergy_T", &ZNAenergy_T);
        globalTree->SetBranchAddress("ZNCenergy_T", &ZNCenergy_T);
        globalTree->SetBranchAddress("ZPAenergy_T", &ZPAenergy_T);
        globalTree->SetBranchAddress("ZPCenergy_T", &ZPCenergy_T);
        globalTree->SetBranchAddress("ZDCAtime_T", ZDCAtime_T);
        globalTree->SetBranchAddress("ZDCCtime_T", ZDCCtime_T);
        globalTree->SetBranchAddress("PIDTPCPion_T", PIDTPCPion_T);
        globalTree->SetBranchAddress("PIDTPCElectron_T", PIDTPCElectron_T);
        globalTree->SetBranchAddress("TPCsignal_T", TPCsignal_T);
        globalTree->SetBranchAddress("TrackP_T", TrackP_T);
        globalTree->SetBranchAddress("TrackEta_T", TrackEta_T);
        globalTree->SetBranchAddress("TrackPhi_T", TrackPhi_T);
        globalTree->SetBranchAddress("TrackPx_T", TrackPx_T);
        globalTree->SetBranchAddress("TrackPy_T", TrackPy_T);
        globalTree->SetBranchAddress("TrackPz_T", TrackPz_T);
        if(newData) globalTree->SetBranchAddress("TrackQ_T", TrackQ_T);
        globalTree->SetBranchAddress("VtxX_T", &VtxX_T);
        globalTree->SetBranchAddress("VtxY_T", &VtxY_T);
        globalTree->SetBranchAddress("VtxZ_T", &VtxZ_T);
        globalTree->SetBranchAddress("VtxContrib_T", &VtxContrib_T);
        globalTree->SetBranchAddress("VtxChi2_T", &VtxChi2_T);
        globalTree->SetBranchAddress("VtxNDF_T", &VtxNDF_T);
        globalTree->SetBranchAddress("SpdVtxX_T", &SpdVtxX_T);
        globalTree->SetBranchAddress("SpdVtxY_T", &SpdVtxY_T);
        globalTree->SetBranchAddress("SpdVtxZ_T", &SpdVtxZ_T);
        globalTree->SetBranchAddress("SpdVtxContrib_T", &SpdVtxContrib_T);
        globalTree->SetBranchAddress("V0Adecision_T", &V0Adecision_T);
        globalTree->SetBranchAddress("V0Cdecision_T", &V0Cdecision_T);
        globalTree->SetBranchAddress("ADAdecision_T", &ADAdecision_T);
        globalTree->SetBranchAddress("ADCdecision_T", &ADCdecision_T);
        globalTree->SetBranchAddress("UBAfired_T", &UBAfired_T);
        globalTree->SetBranchAddress("UBCfired_T", &UBCfired_T);
        globalTree->SetBranchAddress("VBAfired_T", &VBAfired_T);
        globalTree->SetBranchAddress("VBCfired_T", &VBCfired_T);
        globalTree->SetBranchAddress("Ntracklets_T", &Ntracklets_T);
        globalTree->SetBranchAddress("ChipCut_T", &ChipCut_T);
        globalTree->SetBranchAddress("TriggerSPD_T", &TriggerSPD_T);
        globalTree->SetBranchAddress("TriggerTOF_T", &TriggerTOF_T);   
        
        if(useDeltaPhi) globalTree->SetBranchAddress("deltaPhi", &deltaPhi);  

        globalTree->SetBranchAddress("TrackEtaGen_T", &TrackEtaGen_T);  
        globalTree->SetBranchAddress("TrackPhiGen_T", &TrackPhiGen_T);  
        globalTree->SetBranchAddress("TrackPtGen_T", &TrackPtGen_T);  

        globalTree->SetBranchAddress("meanPt",&meanPt);
        globalTree->SetBranchAddress("meanDeltaPhi",&meanDeltaPhi);
    }
    // Declaration of leaf types of GENERATED EVENTS
    Int_t           RunNum_MC_T;
    Float_t         Mass_MC_T;
    Float_t         Pt_MC_T;
    Float_t         Rapidity_MC_T;
    Float_t         Phi_MC_T;
    Int_t           N_MC_T;
    Float_t         Pt1_MC_T;
    Float_t         Phi1_MC_T;
    Int_t           Q1_MC_T;
    Float_t         Eta1_MC_T;
    Float_t         Pt2_MC_T;
    Float_t         Phi2_MC_T;
    Int_t           Q2_MC_T;
    Float_t         Eta2_MC_T;

    void createMCGenTree(TFile *globalFile, string treePath, TTree* = globalTree){
        char ctreePath[100];
        strcpy(ctreePath, treePath.c_str());
        globalTree = (TTree*)globalFile->Get(ctreePath);
        //setting address ONLY of Runs, Mass, Rapidity and Pt to the same address of reald data address
        //in this way I can use the same code to make the cuts for reconstructed and generated events
        globalTree->SetBranchAddress("RunNum_MC_T", &RunNum_T); //same address of real data!
        globalTree->SetBranchAddress("Mass_MC_T", &Mass_T); //same address of real data!
        globalTree->SetBranchAddress("Pt_MC_T", &Pt_T); //same address of real data!
        globalTree->SetBranchAddress("Rapidity_MC_T", &Rapidity_T); //same address of real data!
        globalTree->SetBranchAddress("Phi_MC_T", &Phi_MC_T);
        globalTree->SetBranchAddress("N_MC_T", &N_MC_T);
        globalTree->SetBranchAddress("Pt1_MC_T", &Pt1_MC_T);
        globalTree->SetBranchAddress("Phi1_MC_T", &Phi1_MC_T);
        globalTree->SetBranchAddress("Q1_MC_T", &Q1_MC_T);
        globalTree->SetBranchAddress("Eta1_MC_T", &Eta1_MC_T);
        globalTree->SetBranchAddress("Pt2_MC_T", &Pt2_MC_T);
        globalTree->SetBranchAddress("Phi2_MC_T", &Phi2_MC_T);
        globalTree->SetBranchAddress("Q2_MC_T", &Q2_MC_T);
        globalTree->SetBranchAddress("Eta2_MC_T", &Eta2_MC_T);

        if(useDeltaPhi) globalTree->SetBranchAddress("deltaPhi", &deltaPhi);   

        globalTree->SetBranchAddress("meanPt",&meanPt);
        globalTree->SetBranchAddress("meanDeltaPhi",&meanDeltaPhi);
    }

    
#endif
    