// Example of Wiener-SVD unfolding
// Reference: [to be added]
// Author: Hanyu WEI    March 17, 2017

#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>

#include "TRandom3.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"

#include "WienerSVD.h" //core implementation of Wiener-SVD
#include "Util.h" // utilities

using namespace std;

int main(int argc, char** argv)
{
    // input histograms (matrices, vectors)
    string inputfile;
    // output file 
    string outputfile;
    // Options of  matrix for smoothness, 0 (unitary matrix), 1 (1st derivative matrix), 2 (2nd derivative matrix), else (3rd derivative matrix)
    Int_t C_type; 
    Float_t Norm_type;

    if( argc !=5 )
    {
        std::cout<<"Usage: "<<std::endl;
        std::cout<<"exe inputfile.root output.root matrix_C_type norm_type"<<std::endl;
        exit(1);
    }
    else
    {
        inputfile = argv[1];
        outputfile = argv[2];
        C_type = atoi(argv[3]);
        Norm_type = atof(argv[4]);
    }
    
    std::cout<<"Derivative matrix type: "<<C_type<<std::endl;
    std::cout<<"Signal normalization type (power exponent): "<<Norm_type<<std::endl;

    TFile* f = new TFile(inputfile.c_str(), "READ");
    
    // Input: signal, measure, response, covariance  
    // true signal spectrum or expected signal spectrum from model, used to construct wiener filter 
    TH1D *sig = (TH1D*)f->Get("htrue_signal"); 
    //TH1D *sig = (TH1D*)f->Get("sig"); 
    // measured spectrum
    TH1D *mes = (TH1D*)f->Get("hmeas");
    //TH1D *mes = (TH1D*)f->Get("Rmeasure");
    // response matrix in histogram format
    TH2D* res = (TH2D*)f->Get("hR");
    //TH2D* res = (TH2D*)f->Get("hresponse");
    // covariance matrix of measured spectrum in histrogram format
    TH2D* cov = (TH2D*)f->Get("hcov_tot");
    //TH2D* cov = (TH2D*)f->Get("hcov");

    Int_t n = sig->GetNbinsX();
    Double_t Nuedges[n+1];
    for(int i=0; i<n+1; i++){
        Nuedges[i] = sig->GetBinLowEdge(i+1);
    }

    Int_t m = mes->GetNbinsX();
//    Int_t mXmin = mes->GetXaxis()->GetXmin();
//    Int_t mXmax = mes->GetXaxis()->GetXmax();
 
    // construct vectors (for 1D histogram) and matrices (for 2D histogram) for input
    TVectorD signal(n);
    TVectorD measure(m);
    TMatrixD response(m, n);
    TMatrixD covariance(m, m);
   
    // convert input into mathematical formats, easy and clean to be processed. 
    // Converted defined/implemented in source files, see include/Util.h
    H2V(sig, signal);
    H2V(mes, measure);
    H2M(res, response, kTRUE);
    H2M(cov, covariance, kTRUE);

    // construct to record additinal smearing matrix and wiener filter (diagomal matrix) elements. 
    TMatrixD AddSmear(n,n);
    TVectorD WF(n);
    TMatrixD UnfoldCov(n,n);
    TH2D* smear = new TH2D("smear","Additional Smearing Matirx",n,Nuedges,n,Nuedges);
    TH1D* wiener = new TH1D("wiener","Wiener Filter Vector",n,0,n);
    TH2D* unfcov = new TH2D("unfcov","Unfolded spectrum covariance", n, Nuedges, n, Nuedges);


    // Core implementation of Wiener-SVD
    // Input as names read. AddSmear and WF to record the core information in the unfolding.
    TVectorD unfold = WienerSVD(response, signal, measure, covariance, C_type, Norm_type, AddSmear, WF, UnfoldCov);

    // output and comparison between true/expected signal spectrum with unfolded one
    TFile* file = new TFile(outputfile.c_str(), "RECREATE");
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c","",800,1000);
    c->Divide(1,2,0.0,0.0);
    c->cd(1);
    TH1D* unf = new TH1D("unf","unfolded spectrum",n,Nuedges);
    V2H(unfold, unf);
    unf->Draw();
    sig->Draw("same");
    unf->SetLineColor(kRed);
    unf->SetLineStyle(kDashed);
    TLegend* lg = new TLegend(0.4, 0.5, 0.6, 0.8,"","NDC");
    lg->AddEntry(unf,"Unfolded spectrum","lf");
    lg->AddEntry(sig,"True spectrum","lf");
    lg->SetBorderSize(0);
    lg->SetTextSize(0.08);
    lg->Draw("same");
    c->cd(2);
    TH1D* diff = new TH1D("diff","Fractional difference of unf and signal model",n, Nuedges);
    for(Int_t i=1; i<=n; i++)
    {
        Double_t s2s = 0;
        Double_t u = unf->GetBinContent(i);
        Double_t t = signal(i-1);
        if(t!=0) s2s = u/t - 1;
        else s2s = 1.; // attention to this t=0
        diff->SetBinContent(i, s2s); // in percentage 
    }
    diff->Draw("");
   
    // intrinsic bias (Ac-I)*s_bar formula
    TH1D* bias = new TH1D("bias","intrinsic bias w.r.t. model",n, Nuedges);
    TH1D* bias2 = new TH1D("bias2","intrinsic bias2 w.r.t. unfolded result",n, Nuedges);
    TMatrixD unit(n,n);
    unit.UnitMatrix();
    TVectorD intrinsicbias = (AddSmear - unit)*signal;
    TVectorD intrinsicbias2 = (AddSmear - unit)*unfold;
    for(int i=0; i<n; i++)
    {
        if(signal(i)!=0) intrinsicbias(i) = intrinsicbias(i)/signal(i);
        else intrinsicbias(i)=0.;
        if(unfold(i)!=0) intrinsicbias2(i) = intrinsicbias2(i)/unfold(i);
        else intrinsicbias2(i)=0.;
    }
    V2H(intrinsicbias, bias);
    V2H(intrinsicbias2, bias2);
    bias->Draw("same");
    bias->SetLineColor(kRed);
    TLegend* lg2 = new TLegend(0.4, 0.5, 0.6, 0.8,"","NDC");
    lg2->AddEntry(diff,"Fractional difference of unfolded and sig model","lf");
    lg2->AddEntry(bias,"Intrinsic bias","lf");
    lg2->SetBorderSize(0);
    lg2->SetTextSize(0.08);
    lg2->Draw("same");
    
    
    // diagonal uncertainty
    TH1D* fracError = new TH1D("fracError", "Fractional uncertainty", n, Nuedges);
    TH1D* absError = new TH1D("absError", "absolute uncertainty", n, Nuedges);
    for(int i=1; i<=n; i++)
    {
        fracError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1))/unfold(i-1));
        absError->SetBinContent(i, TMath::Sqrt(UnfoldCov(i-1, i-1)));
    }
    
    /// MSE
    TH1D* MSE = new TH1D("MSE", "Mean Square Error: variance+bias^2", n, Nuedges);
    TH1D* MSE2 = new TH1D("MSE2", "Mean Square Error: variance", n, Nuedges);
    for(int i=0; i<n; i++)
    {
        MSE->SetBinContent(i+1, TMath::Power(intrinsicbias2(i)*unfold(i),2)+UnfoldCov(i,i));
        MSE2->SetBinContent(i+1, UnfoldCov(i,i));
    }

    // convert matrix/vector to histogram and save
    M2H(AddSmear, smear);
    V2H(WF, wiener);
    M2H(UnfoldCov, unfcov);

    unf->Write();
    smear->Write();
    wiener->Write();
    diff->Write();
    bias->Write();
    bias2->Write();
    fracError->Write();
    absError->Write();
    unfcov->Write();
    MSE->Write();
    MSE2->Write();
    c->Write();
    file->Close();

    return 1;
}

