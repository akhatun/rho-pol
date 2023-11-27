#ifndef SAVEDVARFROMMASSFIT_H
#define SAVEDVARFROMMASSFIT_H

    //TTree globl pionter
    TTree *saveFitTree = NULL;
    
    //->fit parameters
    float A;
    float BoverA;
    float Nmumu;
    float mumuCost;
    float mumuExp;
    double *bgCovVec;
    double *soedingCovVec;
    float errBoverA;
    //->background fit statistics
    float chi2DF_bkg;
    int DF_bkg;
    float fitProb_bkg;
    //->soeding fit statistics
    float chi2DF_soed;
    int DF_soed;
    float fitProb_soed;
    //->number of the particles
    float nMuon;
    float nRho;
    float errNMuon;
    float errNRho;
    //->mean pt and phi of the bin
    extern float meanPt;
    extern float meanDeltaPhi;
    //->mass and width of the rho
    float mRho;
    float wRho;
    float err_mRho;
    float err_wRho;
    //parameters used for the fit
    //float lowerLimFit;
    //float upperLimFit;
    //int nBinFit;

    void createSaveFitTree(string name){
        saveFitTree = new TTree (name.c_str(),name.c_str());;

        //defining the branches of the tree
        saveFitTree->Branch("A",&A,"A/F");
        saveFitTree->Branch("BoverA",&BoverA,"BoverA/F");
        saveFitTree->Branch("Nmumu",&Nmumu,"Nmumu/F");
        saveFitTree->Branch("mumuCost",&mumuCost,"mumuCost/F");
        saveFitTree->Branch("mumuExp",&mumuExp,"mumuExp/F");
        saveFitTree->Branch("bgCovVec",&bgCovVec,"bgCovVec/D");
        saveFitTree->Branch("soedingCovVec",&soedingCovVec,"soedingCovVec/D");
        saveFitTree->Branch("errBoverA",&errBoverA,"errBoverA/F");
        saveFitTree->Branch("chi2DF_bkg",&chi2DF_bkg,"chi2DF_bkg/F");
        saveFitTree->Branch("DF_bkg",&DF_bkg,"DF_bkg/I");
        saveFitTree->Branch("fitProb_bkg",&fitProb_bkg,"fitProb_bkg/F");
        saveFitTree->Branch("chi2DF_soed",&chi2DF_soed,"chi2DF_soed/F");
        saveFitTree->Branch("DF_soed",&DF_soed,"DF_soed/I");
        saveFitTree->Branch("fitProb_soed",&fitProb_soed,"fitProb_soed/F");
        saveFitTree->Branch("nMuon",&nMuon,"nMuon/F");
        saveFitTree->Branch("nRho",&nRho,"nRho/F");
        saveFitTree->Branch("errNMuon",&errNMuon,"errNMuon/F");
        saveFitTree->Branch("errNRho",&errNRho,"errNRho/F");
        saveFitTree->Branch("meanPt",&meanPt,"meanPt/F");
        saveFitTree->Branch("meanDeltaPhi",&meanDeltaPhi,"meanDeltaPhi/F");

        saveFitTree->Branch("mRho",&mRho,"mRho/F");
        saveFitTree->Branch("wRho",&wRho,"wRho/F");
        saveFitTree->Branch("err_mRho",&err_mRho,"err_mRho/F");
        saveFitTree->Branch("err_wRho",&err_wRho,"err_wRho/F");

        //saveFitTree->Branch("lowerLimFit",&lowerLimFit,"lowerLimFit/F");
        //saveFitTree->Branch("upperLimFit",&upperLimFit,"upperLimFit/F");
        //saveFitTree->Branch("nBinFit",&nBinFit,"nBinFit/I");
    }
#endif