// Core implementation of Wiener filter SVD
// Author: Hanyu WEI,   March 17, 2017

#ifndef WFUnfold_H
#define WFUnfold_H
#include "TMatrixD.h"
#include "TVectorD.h"

// construction of additional matrix for smoothness/improvement
TMatrixD Matrix_C(Int_t n, Int_t type);

// wiener filter unfolding
// first five parameters are inputs as names read
// AddSmear is additional smearing matrix after unfolding
// WF is the elements of Wiener Filter
TVectorD WFUnfold(TMatrixD Response, TVectorD Signal, TVectorD Measure, TMatrixD Covariance, Int_t C_type, TMatrixD& AddSmear, TVectorD& WF);

#endif


