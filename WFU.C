// Example of SVD unfolding with Wiener Filter
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

#include "WFUnfold.h" //core implementation of Wiener Filter SVD
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

    if( argc !=4 )
    {
        std::cout<<"Usage: "<<std::endl;
        std::cout<<"exe inputfile.root output.root matrix_C_type"<<std::endl;
        exit(1);
    }
    else
    {
        inputfile = argv[1];
        outputfile = argv[2];
        C_type = atoi(argv[3]);
    }

    TFile* f = new TFile(inputfile.c_str(), "READ");
    
    // Input: signal, measure, response, covariance  
    // true signal spectrum or expected signal spectrum from model, used to construct wiener filter 
    TH1D *sig = (TH1D*)f->Get("sig"); 
    // measured spectrum
    TH1D *mes = (TH1D*)f->Get("Rmeasure");
    // response matrix in histogram format
    TH2D* res = (TH2D*)f->Get("hresponse");
    // covariance matrix of measured spectrum in histrogram format
    TH2D* cov = (TH2D*)f->Get("hcov");

    Int_t n = sig->GetNbinsX();
    Int_t sXmin = sig->GetXaxis()->GetXmin();
    Int_t sXmax = sig->GetXaxis()->GetXmax();

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
    TH2D* smear = new TH2D("smear","Additional Smearing Matirx",n,sXmin,sXmax,n,sXmin,sXmax);
    TH1D* wiener = new TH1D("wiener","Wiener Filter Vector",n,sXmin,sXmax);



    // Core implementation of Wiener-Filter SVD
    // Input as names read. AddSmear and WF to record the core information in the unfolding.
    TVectorD unfold = WFUnfold(response, signal, measure, covariance, C_type, AddSmear, WF);


    // output and comparison between true/expected signal spectrum with unfolded one
    TFile* file = new TFile(outputfile.c_str(), "RECREATE");
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c","",800,1000);
    c->Divide(1,2,0.0,0.0);
    c->cd(1);
    TH1D* unf = new TH1D("unf","unfolded spectrum",n,sXmin,sXmax);
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
    TH1D* bias = new TH1D("bias","bias",n,sXmin,sXmax);
    for(Int_t i=1; i<=n; i++)
    {
        Double_t s2s = 0;
        Double_t u = unf->GetBinContent(i);
        Double_t t = signal(i-1);
        if(t!=0) s2s = u/t - 1;
        else s2s = 1.; // attention to this t=0
        bias->SetBinContent(i, 100.*s2s); // in percentage 
    }
    bias->Draw("");
    bias->GetYaxis()->SetRangeUser(-20,20); //adjustable

    // convert matrix/vector to histogram and save
    M2H(AddSmear, smear);
    V2H(WF, wiener);

    unf->Write();
    smear->Write();
    wiener->Write();
    bias->Write();
    c->Write();
    file->Close();

    return 1;
}

