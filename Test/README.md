Test files. All in histogram format.  
Directly used for WFU.C in main directory.  


Test.root: square response matrix.  
Test\_ns.root: non-square response matrix.


1. sig - true signal spectrum, S.
2. measure - expected measured spectrum, M, equal to true signal spectrum multiplied with response matrix [M = R*S].
3. hcov - covariance matrix of measured spectrum. Same method of converting a matrix to histogram as hresponse. (Statistical only - Poisson)
4. Rmeasure - one measured spectrum with fluctuation.
5. hresponse - response matrix, R.   
    5.1 In histogram format, x/y axis refer to matrix row/column, e.g. hresponse->GetBinContent(i+1, j+1) = R(i, j), i,j belong to [0, n-1].  
    5.2 Column normalized.  
    5.3 Rotate the histrogram 90 degree to the right, it is consistent with the matrix R in writting format.   
