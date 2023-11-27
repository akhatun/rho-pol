#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TSystem.h"
#include "TDatabasePDG.h"
#include "TPad.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TFitResult.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TSystemDirectory.h"
#include "vector"
#include <string>
#include <cstdint>
#include "libraries/savedVarFromMassFit.h"
#include "libraries/treeVarDeclarationBinned.h" 
#include "TROOT.h"

//definition of constants
#ifndef PARTICLES
#define PARTICLES
TDatabasePDG *pdgdat = TDatabasePDG::Instance(); 
//TParticlePDG *partRho = pdgdat->GetParticle( 113 );
TParticlePDG *partPion = pdgdat->GetParticle( 211 );
const float mPi = partPion->Mass(); //GeV
const double mRho1 = 0.7692; //GeV
const double wRho1 = 0.1515; //GeV
const double mOmega1 = 0.7826; //GeV
const double wOmega1 = 0.0085; //GeV -- note that in the AN of the cross section they used 0.01
#endif

//FUNCTIONS FOR THE BACKGROUND
// power background from mumu
double MuMuBkgd(double x, double p1, double p2){
    return p1*pow(x,p2);
}

// fit to power background from mumu
double fitMuMuBkgd(double *x, double *par){
    return MuMuBkgd(x[0],par[0],par[1]);
}

//FUNCTIONS TO COMPOSE THE SOEDING MODEL FUNCTION
// effective width in the Breit-Wiegner in Soeding's model  
double GammaBW(double mPiPi, double mRho, double wRho, double BR_pipi){
    // assuming that mPiPi>0 and mRho^2>4mPi^2
    // is checked in the main BW function
    double fracNum = mPiPi*mPiPi-4.0*mPi*mPi;
    double fracDen = mRho*mRho-4.0*mPi*mPi;
    double frac = fracNum/fracDen;
    double frac32 = frac*sqrt(frac);
    return (BR_pipi*wRho*mRho*frac32/mPiPi);
}

// real part of Breit-Wiegner in Soeding's model
double ReBW(double mPiPi, double mRho, double wRho, double BR_pipi){
    // assuming that mPiPi>0 and mRho>mPi
    // is checked in the main BW function
    double effWidth = GammaBW(mPiPi, mRho, wRho, BR_pipi);
    double a1 = sqrt(mPiPi*mRho*effWidth);
    double a2 = mPiPi*mPiPi-mRho*mRho;
    double a3 = mRho*effWidth;
    double reNum = a1*a2;
    double reDen = a2*a2+a3*a3;
    return (reNum/reDen);
}

// imaginary part of Breit-Wiegner in Soeding's model
double ImBW(double mPiPi, double mRho, double wRho, double BR_pipi){
    // assuming that mPiPi>0 and mRho>mPi
    // this has to be checked in the main BW function
    double effWidth = GammaBW(mPiPi, mRho, wRho, BR_pipi);
    double a1 = sqrt(mPiPi*mRho*effWidth);
    double a2 = mPiPi*mPiPi-mRho*mRho;
    double a3 = mRho*effWidth;
    double imNum = a1*a3;
    double imDen = a2*a2+a3*a3;
    return (imNum/imDen);
}

// interference term in Soeding's model 
double InBW(double mPiPi, double A, double B, double mRho, double wRho, double BR_pipi){
    return (2.0*A*B*ReBW(mPiPi, mRho, wRho, BR_pipi));
}

//function to graph the interference term
double showInterference(double *x, double *par){
    return InBW(x[0],par[0],par[1], par[2], par[3],1); // is for the rho-> BR=1
}

// Soeding's model
double Soeding(double mPiPi, double A, double B, double mRho, double wRho){
    if(mPiPi<=0){
        cout<<"ERROR! mass of pion pair = 0"<<endl;
        return 0;
    }
    if(mRho<mPi){
        //cout<<"ERROR! mRho<mPi "<<endl;
    }
    double BR_pipi = 1; //for the rho
    double re = ReBW(mPiPi, mRho, wRho, BR_pipi); // real part of BW
    double im = ImBW(mPiPi, mRho, wRho, BR_pipi); // imaginary part of BW
    double in = 2.0*A*B*re;  // interference 
    double f = A*A*re*re+B*B+A*A*im*im+in;
    return f;
}

// function to fit using the Soeding's model
double fitSoedingMuMu(double *x, double *par){
    double mPiPi = x[0];
    double A = par[0];
    double BonA = par[1]; //using B/A instead of B to have the error computation yet
    double nMuMu = par[2];
    double p1 = par[3];
    double p2 = par[4];
    double mRho = par[5];
    double wRho = par[6];
    double BR_pipi = 1; //for the rho
    double signal = Soeding(mPiPi,A,BonA*A,mRho,wRho); //using B/A istead of B->is the correct way?
    double bkgd = nMuMu*MuMuBkgd(mPiPi,p1,p2);
    return (signal+bkgd);  
}

double squareExpr(double vec[], int n){
    double sumS = 0;
    for(int i=0; i<n; i++){
        sumS = sumS + vec[i]*vec[i];
        for(int k=i+1; k<n; k++){
            sumS = sumS + 2*vec[i]*vec[k];
        }
    }
    return sumS;
}

//Soding model with omega
double SodingWithOmega(double mPiPi, double A, double B, double mRho, double wRho, double C, double phiMix, double mOmega, double wOmega){
    if(mPiPi<=0){
        cout<<"ERROR! mass of pion pair = 0"<<endl;
        return 0;
    }
    double BR_rho = 1;
    double BR_omega = 0.0153;
    double reRho = ReBW(mPiPi, mRho, wRho, BR_rho); // real part of BW of the Rho
    double imRho = ImBW(mPiPi, mRho, wRho, BR_rho); // imaginary part of BW of the Rho
    double reOmega = ReBW(mPiPi, mOmega, wOmega, BR_omega); // real part of BW of the Omega
    double imOmega = ImBW(mPiPi, mOmega, wOmega, BR_omega); // imaginary part of BW of the Omega
    double sinPhi = TMath::Sin(phiMix); //sin of the mixing phase
    double cosPhi = TMath::Cos(phiMix); //cos of the mixing phase
    
    //computing the square module
    double real[4] = {A*reRho, B, C*reOmega*cosPhi, -imOmega*sinPhi};
    double realSquare = squareExpr(real,4);

    double imm[3] = {A*imRho, cosPhi*imOmega, sinPhi*reOmega};
    double immSquare = squareExpr(imm,3);

    return realSquare + immSquare;

}

//function to fit the mass distrib with soding model with the Omega
double fitSoedingWithOmega(double *x, double *par){
    double mPiPi = x[0];
    double A = par[0];
    double BonA = par[1]; //using B/A instead of B to have the error computation yet
    double nMuMu = par[2];
    double p1 = par[3];
    double p2 = par[4];
    double mRho = par[5];
    double wRho = par[6];
    double C = par[7];
    double phiMix = par[8];
    double mOmega = par[9];
    double wOmega = par[10];

    double signal = SodingWithOmega(mPiPi,A,BonA*A,mRho,wRho,C,phiMix,mOmega,wOmega); //using B/A istead of B
    double bkgd = nMuMu*MuMuBkgd(mPiPi,p1,p2);
    return (signal+bkgd);  
}

//show the interference between rho and omega
double showRhoOmegaInt(double *x, double *par){
    double mPiPi = x[0];
    double mRho = mRho1;
    double wRho = wRho1;
    double mOmega = mOmega1;
    double wOmega = wOmega1;
    double A = par[0];
    double C = par[1];
    double phiMix = par[2];
    double BR_rho = 1;
    double BR_omega = 0.0153;
    double reRho = ReBW(mPiPi, mRho, wRho, BR_rho); // real part of BW of the Rho
    double imRho = ImBW(mPiPi, mRho, wRho, BR_rho); // imaginary part of BW of the Rho
    double reOmega = ReBW(mPiPi, mOmega, wOmega, BR_omega); // real part of BW of the Omega
    double imOmega = ImBW(mPiPi, mOmega, wOmega, BR_omega); // imaginary part of BW of the Omega
    double sinPhi = TMath::Sin(phiMix); //sin of the mixing phase
    double cosPhi = TMath::Cos(phiMix); //cos of the mixing phase

    return 2*(A*reRho + C*cosPhi*reOmega + C*sinPhi*imOmega);
}

// function to fill the vector passed with the mass
void fillVector(string filePath, string treePath, bool realD, vector<float>& vec, string varName = "Mass_T"){

    //opening the file that stores the tree 
    TFile *file = new TFile(filePath.c_str());
    //calling the function that allow to access the TTree accordingly to the tipe of the data in the file
    if(realD){
        createDataTree(file,treePath.c_str());
    }else{
        createDataTree(file,treePath.c_str());
    }

    for(Long64_t ev=0; ev<globalTree->GetEntries(); ev++){
        globalTree->GetEvent(ev);
        if(varName=="Mass_T") vec.push_back(Mass_T);
        if(varName=="Pt_T") vec.push_back(Pt_T);
        if(varName=="deltaPhi") vec.push_back(deltaPhi);
    }

    file->Close();
}

//function that fills the vactor passed using the neutron class
void fillVectorWithNeutronClass(string filePath, string treePath, string dataNature, vector<float>& vec, string varName , string neutronClass){
    
    //opening the file that stores the tree 
    TFile *file = new TFile(filePath.c_str());
    //calling the function that allow to access the TTree accordingly to the tipe of the data in the file
    if(dataNature=="data") createDataTree(file,treePath.c_str());
    if(dataNature=="reco") createMCRecoTree(file,treePath.c_str());
    if(dataNature=="gen") createMCGenTree(file,treePath.c_str());

    bool fillV = false;
    bool neutron_A = false;
    bool neutron_C = false;
    
    int c_noSel = 0;
    int c_XnXn = 0;
    int c_0n0n = 0;
    int c_Xn0n = 0;
    int c_12n0n = 0;
    int c_Mn0n = 0;
    for(Long64_t ev=0; ev<globalTree->GetEntries(); ev++){
       
        fillV = false;
        neutron_A = false;
        neutron_C = false;

        globalTree->GetEvent(ev);
        
        for(int i=0; i<4; i++){
            if(ZDCAtime_T[i]<2) neutron_A = true;
            if(ZDCCtime_T[i]<2) neutron_C = true;
        }

        if(neutronClass == "noSelection"){
            fillV = true;
            c_noSel++;
        }
        else if(neutronClass == "0n0n"){
            if(neutron_A == false){
                if(neutron_C == false){
                    fillV = true;
                    c_0n0n++;
                }
            } 
        }
        else if(neutronClass == "0nXn" || neutronClass == "Xn0n"){
            if(neutron_A ^ neutron_C){ // ^ is the XOR
                fillV = true; 
                c_Xn0n++;
            } 
        }
        else if(neutronClass == "XnXn"){
            if(neutron_A && neutron_C){
                fillV = true;
                c_XnXn++;
            } 
        }
        // 1 or 2 neutron at one side and no neutron at the other
        else if(neutronClass == "12n0n"){
            // 1,2 n in ZDCA, no n in ZDCC
            if(neutron_A == true && neutron_C == false){
                if((ZNAenergy_T/2510)<2.5){
                    fillV = true;
                    c_12n0n++;
                }
            // 0 n in ZDCA, 1 or 2 n in ZDCC
            }else if(neutron_A == false && neutron_C == true){
                if((ZNCenergy_T/2510)<2.5){
                    fillV = true;
                    c_12n0n++;
                }
            }
        }
        // more than 3 neutrons at one side and no neutron(s) at the other
        else if(neutronClass == "Mn0n"){
            // more than 3 n in ZDCA, no n in ZDCC
            if(neutron_A == true && neutron_C == false){
                if((ZNAenergy_T/2510)>2.5){
                    fillV = true;
                    c_12n0n++;
                }
            // 0 n in ZDCA, more than 3 n in ZDCC
            }else if(neutron_A == false && neutron_C == true){
                if((ZNCenergy_T/2510)>2.5){
                    fillV = true;
                    c_12n0n++;
                }
            }
        }

        if(fillV){
            if(varName=="Mass_T") vec.push_back(Mass_T);
            if(varName=="Pt_T") vec.push_back(Pt_T);
            if(varName=="deltaPhi") vec.push_back(deltaPhi);
        }
        
    }

    file->Close();
}

// function to set the mean value of phi and pt inside a bin
void findBinPar(string filePath, string treePath, bool realD, float &phi, float &pt){
    //opening the file that stores the tree 
    TFile *file = new TFile(filePath.c_str());
    //calling the function that allow to access the TTree accordingly to the tipe of the data in the file
    if(realD){
        createDataTree(file,treePath.c_str());
    }else{
        createDataTree(file,treePath.c_str());
    }
    //reading only the first element since the value is the same for all events in the same bin
    globalTree->GetEvent(0);
    pt = meanPt;
    phi = meanDeltaPhi;

    file->Close();
}

// function that fill a TH1 with the content of the passed vector
// dy defalt set the x-axis title for the invariant mass
void fillHistoWithVector(TH1 *h, vector<float> v, string xTitle = "Invariant Mass [GeV/c^{2}]"){
    for(unsigned int i=0; i<v.size(); i++){
        h->Fill(v[i]);
    }
    h->Sumw2();
    h->GetXaxis()->SetTitle(xTitle.c_str());
    h->GetYaxis()->SetTitle("#count");
}

//function that performs the mass fit and saves the results in a Tree
void performFit(string fileName,string folder, TFile *imageFile,string graphicFile, TFile *saveFile, string resulFile,string treeName, string neutronClass, double &nEntr){

    //no pop up of the canvas
    //gROOT->SetBatch(kTRUE);

    //string useful to form the name of the file and folders in which are the binned data
    string filePath;
    string treePath;

    //DEFINING and FILLING VECTORS to allow the creation of fit histogram with any binning
    //real data
    vector<float> real_vecMass;
    filePath = fileName + "/real";
    treePath = folder + "/real/binnedTree";
    cout<<"fileName \t "<<fileName<<endl;
    cout<<"treePath \t "<<treePath<<endl;
    fillVectorWithNeutronClass(fileName, treePath.c_str(), "data", real_vecMass,"Mass_T",neutronClass);
    
    //MC rho data
    //reco data
    vector<float> reco_vecMass;
    filePath = fileName + "/reco";
    treePath = folder + "/reco/binnedTree";
    fillVector(fileName, treePath.c_str(), "reco", reco_vecMass);
    //gen data
    vector<float> gen_vecMass;
    filePath = fileName + "/gen";
    treePath = folder + "/gen/binnedTree";
    fillVector(fileName, treePath.c_str(), "gen", gen_vecMass);
    //Mc gamma->mu data
    vector<float> muonReco_vecMass;
    filePath = fileName + "/muonReco";
    treePath = folder + "/muonReco/binnedTree";
    fillVector(fileName, treePath.c_str(), "reco", muonReco_vecMass);


    //opening the file to save the results on the tree
    saveFile = new TFile(Form("%s/%s",neutronClass.c_str(),resulFile.c_str()),"update");

    //defining the tree in which the results will be saved
    //using a function in the header savedVarFromMassFit.h
    createSaveFitTree(treeName);    

    //declaring the variables that will be used in the fit
    double xMin;
    double xMax;
    int nBinForFit;
    //defining the limits for the random fits
    double lowLim[2] = {0.6, 0.65};
    double upperLim[2] = {0.95, 1.2};
    int binNumber[2] = {25, 35};

    //to study sistematic I do more fit instead of 1
    int nCycles = 1;
    for(int nc=0; nc<nCycles; nc++){
        //printing the progresses in the evaluation
        if(nc%(nCycles/10)==0) cout<<"evaluation at "<<(double)nc/nCycles*100<<"%"<<endl;

        //defing random ranges and number of bins
        //using the first cycle to perform the fit with optimized values
        if(nc==0){
            // old values used with phi comp with charge
            //xMin = 0.6;
            //xMax = 1.1;
            //nBinForFit = 30;

            // new values
            xMin = 0.6;
            xMax = 1.1;
            //nBinForFit = (xMax-xMin)*1000/5; //5 MeV
            nBinForFit = (xMax-xMin)*1000/50; //50 MeV
            //nBinForFit = (xMax-xMin)*1000/10; //10 MeV
            //nBinForFit = 30;

        }else{
            xMin = lowLim[0] + gRandom->Rndm()*( lowLim[1]-lowLim[0] );
            xMax = upperLim[0] + gRandom->Rndm()*( upperLim[1]-upperLim[0] );
            nBinForFit = binNumber[0] + gRandom->Rndm()*( binNumber[1]-binNumber[0] );
        }

        //saving the mean phi and pt inside the bin 
        findBinPar(fileName, treePath.c_str(), true, meanDeltaPhi,meanPt);

        //defining and filling the histograms useful for the fit
        TH1D *hMass_real = new TH1D("hMass_real","Real data mass distribution",nBinForFit,xMin,xMax);
        fillHistoWithVector(hMass_real, real_vecMass);

        TH1D *hMass_reco = new TH1D("hMass_reco","MC Reconstructed data mass distribution",nBinForFit,xMin,xMax);
        fillHistoWithVector(hMass_reco, reco_vecMass);

        TH1D *hMass_gen = new TH1D("hMass_gen","MC Generated data mass distribution",nBinForFit,xMin,xMax);
        fillHistoWithVector(hMass_gen, gen_vecMass);

        TH1D *hMass_muonReco = new TH1D("hMass_muonReco","MC muon reconstructed data mass distribution",nBinForFit,xMin,xMax);
        fillHistoWithVector(hMass_muonReco, muonReco_vecMass);

        //width of the bin of the histograms needed to obtain the number of particle form integral
        float binWidth = hMass_gen->GetXaxis()->GetBinWidth(1); 

    /*
        //INPUT CONTROL PLOTS 
        TCanvas *controlInput = new TCanvas("controlInput","controlInput",0,0,900,600);
	    controlInput->Divide(2,2);
	    controlInput->cd(1); hMass_real->DrawCopy();
	    controlInput->cd(2); hMass_reco->DrawCopy();
	    controlInput->cd(3); hMass_gen->DrawCopy();
	    controlInput->cd(4); hMass_muonReco->DrawCopy();
    */

        ///AxE pion for the fit
	    TH1D* hMass_axePion = (TH1D*) hMass_reco->Clone("hMass_axePion");
        hMass_axePion->SetTitle("AxE under pion hypotesis");
        hMass_axePion->Divide(hMass_reco,hMass_gen,1,1,"B");

        //correcting the real data
	    TH1D* hMass_wAxE = (TH1D*) hMass_reco->Clone("hMass_wAxE");	
        hMass_wAxE->Divide(hMass_real,hMass_axePion);
        hMass_wAxE->SetTitle("Real data inv. mass corrected with AxE");
        
        //compute the correction
	    TH1D* hMass_muonBkgCorr = (TH1D*)  hMass_muonReco->Clone("hMass_muonBkgCorr");	
        hMass_muonBkgCorr->Divide(hMass_muonReco,hMass_axePion);
        hMass_muonBkgCorr->SetTitle("Muon background distribution corrected");

    /*
        //CORRECTIONS CONTROL PLOTS 
        TCanvas *controlAxE = new TCanvas("controlAxE","controlAxE",0,0,800,600);
	    controlAxE->Divide(2,2);
	    controlAxE->cd(1); hMass_axePion->Draw();
	    controlAxE->cd(2); hMass_wAxE->Draw();
	    controlAxE->cd(3); hMass_muonBkgCorr->Draw();
    */

        //defining the TMatrixD that will contain the covariance matrix
        TMatrixDSym *covBgMatrix = new TMatrixDSym(2);
        TMatrixDSym *covSoedingMatrix = new TMatrixDSym(11);

        //START PERFORMING THE FITS
        //canvas that will contain the results of the fits
        string tcName = folder + + "__" + to_string(nc);
        TCanvas *cf = new TCanvas(tcName.c_str(),tcName.c_str(),400,100,900,600);

        //Fitting the BACKGROUND template
        TF1 *bkgFit = new TF1("bkgFit",fitMuMuBkgd,xMin,xMax,2);
	    bkgFit->SetParameter(0,hMass_muonBkgCorr->GetMaximum());
	    bkgFit->SetParameter(1,-1);

        //doing the fit only if there are points->otherwise set background to 0
        if(( hMass_muonBkgCorr->GetEntries() )>2){

            TFitResultPtr resFitBg =  hMass_muonBkgCorr->Fit(bkgFit,"SRM+Q");
            gStyle->SetOptFit(1);
            *covBgMatrix = resFitBg->GetCovarianceMatrix();
            //covBgMatrix->Print();
            bgCovVec = covBgMatrix->GetMatrixArray();

            //filling the structur with the statistic of the background fit
            chi2DF_bkg = bkgFit->GetChisquare()/(float)bkgFit->GetNDF();
            DF_bkg = bkgFit->GetNDF();
            fitProb_bkg = bkgFit->GetProb();

        }else{
            bkgFit->SetParameters(0,0);
        }

        //Fitting the TOTAL DISTRIB WITH THE SOEDING MODEL
        TF1 *soedFit = new TF1("soedFit",fitSoedingWithOmega,xMin,xMax,11);
        soedFit->SetParNames("A","B/A","n mumu","mumu const","mumu exponent","mRho","wRho","C","#Phi","mOmega","wOmega");
	    double par_A = 20;
        soedFit->SetParLimits(0,0,1000);
        soedFit->SetParameters(par_A,-0.5,0.1*hMass_wAxE->GetMaximum(),-10,-1.2,mRho1,wRho1,5.5,1.1,mOmega1,wOmega1); 
        soedFit->FixParameter(3,bkgFit->GetParameter(0));
        soedFit->FixParameter(4,bkgFit->GetParameter(1));
        soedFit->FixParameter(5,mRho1);
        soedFit->FixParameter(6,wRho1);
        soedFit->FixParameter(9,mOmega1);
        soedFit->FixParameter(10,wOmega1);
        double percVar = 0.3;
	    soedFit->SetParLimits(2,0,hMass_wAxE->GetMaximum());

        TFitResultPtr resFitSoeding =  hMass_wAxE->Fit(soedFit,"SRM+");

        gStyle->SetOptFit(1);

        *covSoedingMatrix = resFitSoeding->GetCovarianceMatrix();

        //for saving COV MATRIX ARRAY of SOEDING fit in the tree
        soedingCovVec = covSoedingMatrix->GetMatrixArray();

        double minIntegral = 2*mPi;
        double maxIntegral = mRho1 + 5*wRho1;

        //for saving STATISTICS of SOEDING fit in the tree
        chi2DF_soed = soedFit->GetChisquare()/(float)soedFit->GetNDF();
        DF_soed = soedFit->GetNDF();
        fitProb_soed = soedFit->GetProb();

        //for saving PARAMETERS of SOEDING fit in the tree
        A = soedFit->GetParameter(0);
        BoverA = soedFit->GetParameter(1);
        Nmumu = soedFit->GetParameter(2);
        mumuCost = soedFit->GetParameter(3);
        mumuExp = soedFit->GetParameter(4);
        errBoverA = soedFit->GetParError(1);

        //DEFINING THE FUNCTIONS TO SHOW THE VARIOUS PART OF THE FIT
        //Interference part
        TF1 *intCurve = new TF1("intCurve",showInterference,xMin,xMax,4);
        intCurve->SetLineColor(3);
        intCurve->SetLineWidth(2);
        intCurve->SetLineStyle(2); //dashed
        //parameters are A and B
        intCurve->SetParameters(soedFit->GetParameter(0),soedFit->GetParameter(1)*soedFit->GetParameter(0),mRho1,wRho1);
        
        //Breit-Wigner part
        //TF1 *breitWigner = new TF1("breitWigner",fitSoedingWithOmega,xMin,xMax,11);
        TF1 *breitWigner = new TF1("breitWigner",fitSoedingMuMu,xMin,xMax,7);
        breitWigner->SetLineColor(6);
        breitWigner->SetLineWidth(2);
        breitWigner->SetLineStyle(2); 
        //the only parameter that is non null is A
        //breitWigner->SetParameters(soedFit->GetParameter(0),0,0,0,0,mRho1,wRho1,0,0,mOmega1,wOmega1);
        breitWigner->SetParameters(soedFit->GetParameter(0),0,0,0,0,mRho1,wRho1);

        //defining sub covariance matrix: for the bw
        TMatrixDSym bwCovMat = covSoedingMatrix->GetSub(0,0,0,0);

        //obtaining the number of rho0 
        nRho = breitWigner->Integral(minIntegral,maxIntegral)/(binWidth); 
        //and the statistical error
        double errBW = breitWigner->IntegralError(minIntegral,maxIntegral,breitWigner->GetParameters(),bwCovMat.GetMatrixArray())/(binWidth); 
        errNRho = errBW;
        
        //estimating the error using the relative error on the A parameter since it is the only free
        //double errRel = soedFit->GetParError(0)/soedFit->GetParameter(0);
        //double errInt = 2*errRel*nRho; //the 2 derives from the error propagation
        //errNRho = errInt;

        //muon background part
        TF1 *muonBackground = new TF1("muonBackground",fitMuMuBkgd,xMin,xMax,2);
        muonBackground->SetLineColor(4);
        muonBackground->SetLineWidth(2);
        muonBackground->SetLineStyle(2); 
        //par[0]=amplitude=bacground const*#muons  par[1]=slope
        muonBackground->SetParameters(bkgFit->GetParameter(0)*soedFit->GetParameter(2),bkgFit->GetParameter(1));

        //breit wigner for the omega
        TF1 *bwOmega = new TF1("bwOmega",fitSoedingWithOmega,xMin,xMax,11);
        bwOmega->SetParameters(0,0,0,0,0,mRho1,wRho1,soedFit->GetParameter(7),soedFit->GetParameter(8),mOmega1,wOmega1);
        bwOmega->SetLineColor(kCyan -4);
        bwOmega->SetLineWidth(2);
        bwOmega->SetLineStyle(2); 

        //interference between the rho and the omega
        TF1 *rhoOmegaInt = new TF1("rhoOmegaInt",showRhoOmegaInt,xMin,xMax,3);
        rhoOmegaInt->SetLineColor(kOrange +8);
        rhoOmegaInt->SetLineWidth(2);
        rhoOmegaInt->SetLineStyle(2);
        rhoOmegaInt->SetParameters(soedFit->GetParameter(0),soedFit->GetParameter(7),soedFit->GetParameter(8)); 

        //SAVING THE RESULTS on the tree
        saveFitTree->Fill();

        //opening the file to save the images of the fit -> results saved in a folder corresponding to the neutron class
        imageFile = new TFile(Form("%s/%s",neutronClass.c_str(),graphicFile.c_str()),"update");

        //drawing the fits and his parts
        cf->Clear();
	    gStyle->SetOptFit(1);
	    cf->Divide(2,1);

        TLegend *legend[2];
        for(int i=1; i<=2; i++){
            if(i==2) gPad->SetLogy();
            cf->cd(i);
	        hMass_wAxE->Draw();
            soedFit->Draw("same");
            intCurve->Draw("same");
            breitWigner->Draw("same");
            muonBackground->Draw("same");
            bwOmega->Draw("same");
            rhoOmegaInt->Draw("same");

            //add legend
            legend[i-1] = new TLegend(0.7,0.77,0.98,0.94);
            legend[i-1]->SetBorderSize(0);
            legend[i-1]->AddEntry(soedFit,"fit curve","l");
            legend[i-1]->AddEntry(intCurve,"#rho - #pi#pi interference","l");
            legend[i-1]->AddEntry(breitWigner,"BW of the #rho","l");
            legend[i-1]->AddEntry(muonBackground,"muon background","l");
            legend[i-1]->AddEntry(bwOmega,"BW of the #omega","l");
            legend[i-1]->AddEntry(rhoOmegaInt,"#rho - #omega interferce","l");
            legend[i-1]->Draw();
        }
        cf->Write();
        imageFile->Write("",TObject::kOverwrite);
        imageFile->Close();

        //deleting the histograms to avoid memory leak
        hMass_real->Delete();
        hMass_reco->Delete();
        hMass_gen->Delete();
        hMass_muonReco->Delete();

    } //end of the cycle to perform more fits
    
    //saving and closing the file with the trees
    saveFile->Write();
    saveFile->Close();

} //end of performFit


//main function: cycles on the files and calls the function that
// performs the fits and saves the results in a tree
void massFitInBin(int phiBins = 15, int ptBins = 1, string neutronClass = "noSelection"){

    //calculating the running time
    TStopwatch timer; 
    timer.Start();

    //name that have the file with the bin in phi and pt selected in input
    string fileName = to_string(phiBins)+"phi_"+to_string(ptBins)+"pt";
    //path to arrive at the file 
    string filePath = "binnedData/" + fileName +".root"; //using phi computed using the charge
    
    ifstream hfile(filePath.c_str());
    if(!hfile){
        cout<<"non esiste ancora un file con "+to_string(phiBins) +" bin in phi e "+to_string(ptBins)+" bin in pt"<<endl;
        return;
    }

    //costruction of the name of the file that will store the trees and the folders with the graphs
    string resulFile = "massFitRes_" + fileName + ".root";
    string graphicFile = "images_" + fileName + ".root";

    //create a folder for the neutron class in which the results are saved
    gSystem->mkdir(neutronClass.c_str());

    //definition of the file taht will contain the folders with the graphs
    //TFile *imageFile = new TFile(graphicFile.c_str(),"recreate");
    TFile *imageFile = new TFile(Form("%s/%s",neutronClass.c_str(),graphicFile.c_str()),"recreate");
    imageFile->Close();

    //definition of the file that will contain the tree with the fit informations
    //TFile *saveFile = new TFile(resulFile.c_str(),"recreate");
    TFile *saveFile = new TFile(Form("%s/%s",neutronClass.c_str(),resulFile.c_str()),"recreate");
    saveFile->Close();
    
    //ID of the bin
    string ID;
    
    //number of entries of the histo
    double vecNEntr[phiBins];

    //loop on the pt and phi bins
    for(int j=1; j<=ptBins; j++){
        for(int i=1; i<=phiBins; i++){
            if(i>0){//i>0
                ID = "phi"+to_string(i)+"_pt"+to_string(j);
                string resultTree = "tree_" + ID; 
                performFit(filePath, ID, imageFile, graphicFile,saveFile,resulFile,resultTree, neutronClass, vecNEntr[i-1]);
            }
            
        }
    }

    timer.Stop();
    timer.Print();
} //end of the main function