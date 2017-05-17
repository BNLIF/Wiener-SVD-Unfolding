#include "Util.h"

#include <iostream>
#include "stdlib.h"

#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorD.h"

using namespace std;

TMatrixD Matrix(Int_t row, Int_t column)
{
    TMatrixD MMM(row, column);
    for(Int_t i=0; i<row; i++)
    {
        for(Int_t j=0; j<column; j++)
        {
            Double_t eee;
            cout<<"Row:"<<i+1<<" Column:"<<j+1<<endl;
            cin>>eee;
            MMM(i, j) = eee;
        }
    }
    MMM.Print();
    return MMM;
}

void MatrixMatirx(TMatrixD M1, TMatrixD M2)
{
    cout<<"Matrix 1 ========="<<endl;
    M1.Print();
    cout<<"Matrix 2 ========="<<endl;
    M2.Print();
    TMatrixD M = M1*M2;
    cout<<"Matrix 1*2 ========="<<endl;
    M.Print();
}

TVectorD Vector(Int_t row)
{
    TVectorD V(row);
    for(Int_t i=0; i<row; i++)
    {
        Double_t eee;
        cout<<"Row "<<i+1<<endl;
        cin>>eee;
        V(i) = eee;
    }
    V.Print();
    return V;
}


void MatrixVector(TMatrixD M, TVectorD V)
{
    cout<<"Matrix 1 ========="<<endl;
    M.Print();
    cout<<"Vector 1 ========="<<endl;
    V.Print();
    cout<<"Matrix * Vector =========="<<endl;
    TVectorD MM = M*V;
    MM.Print();
}


void SVD(TMatrixD M)
{
    Int_t nrow = M.GetNrows();
    Int_t ncol = M.GetNcols();

    TDecompSVD svd(M);
    TMatrixD U = svd.GetU();
    TMatrixD V (TMatrixD::kTransposed, svd.GetV());
    TVectorD S = svd.GetSig();

    cout<<"U ========="<<endl;
    U.Print();
    TMatrixD uu(U, TMatrixD::kMultTranspose, U);
    uu.Print();
    cout<<"V ========="<<endl;
    V.Print();
    TMatrixD vv(V, TMatrixD::kMultTranspose, V);
    vv.Print();
    cout<<"D ========="<<endl;
    TMatrixD D(nrow, ncol);
    for(Int_t i=0; i<nrow; i++)
    {
        for(Int_t j=0; j<ncol; j++)
        {
            D(i, j) = 0;
            if(i == j)
            {
                D(i, j) = S(i);
            }
        }
    }
    D.Print();

    cout<<"Orignal ========="<<endl;
    M.Print();
    cout<<"Reassembly ========"<<endl;
    TMatrixD MM = U*D*V;
    MM.Print();
    MM.Draw("colz");
}


void H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn)
{
    // Fill 2D histogram into matrix
    // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        for(Int_t j=0; j<histo->GetNbinsY(); j++)
        {
            if(rowcolumn) mat(i, j) = histo->GetBinContent(i+1, j+1);
            else mat(j, i) = histo->GetBinContent(i+1, j+1);
        }
    }
}

void H2V(const TH1D* histo, TVectorD& vec)
{
    // Fill 1D histogram into matrix
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        vec(i) = histo->GetBinContent(i+1);
    }
}

void M2H(const TMatrixD mat, TH2D* histo)
{
    // Fill matrix to histogram
    for(Int_t i=0; i<mat.GetNrows(); i++)
    {
        for(Int_t j=0; j<mat.GetNcols(); j++)
        {
            histo->SetBinContent(i+1, j+1, mat(i, j));
        }
    }
}


void V2H(const TVectorD vec, TH1D* histo)
{
    // Fill vector to histogram,
    for(Int_t i=0; i<vec.GetNrows(); i++)
    {
        histo->SetBinContent(i+1, vec(i));
    }
}
