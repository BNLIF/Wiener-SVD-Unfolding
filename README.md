Implementaion of Wiener-SVD unfolding    
Reference: 	JINST, 12, P10002 (2017) \[arXiv:1705.03568](https://arxiv.org/abs/1705.03568)  

Pre-requiste: ROOT, some ROOT classes (e.g. TMatrix) will be used.  
Core algorithms implementation in src/  
Header files in include/  
Use Makefile to build executables.  

INSTRUCTION

Input: 
1. true/expected signal spectrum (models); 
2. measured spectrum; 
3. covariance matrix of the uncertainty of the measured spectrum; 
4. response matrix from signal (corresponding to matrix column, TH2D Y-axis) to measurement (corresponding to matrix row, TH2D X-axis), M = R*S; 
5. Choice of addional matrix for smoothness, e.g. 2nd derivative matrix.  

Output:
1. Unfolded spectrum
2. Additional smearing matrix
3. Wiener filter
4. Covariance matrix of the unfolded spectrum


An example of application can be found in Example.C  
Usage:  
make clean;  
make;  
Example input.root output.root 2 0 
(2 means: 2nd derivative matrix; 0 means the number of measured events in each bin normalized by m(i)^0)

* 2: 2nd-order derivative matrix
* input.root:
** TH1D htrue_signal; a model of true signal
** TH1D hmeas; measured spectrum
** TH2D hcov_tot; covariance matrix of measured spectrum uncertainty
** TH2D hR; resposne matrix (X-axis: measured, Y-axis: true)
* output.root
** TH1D unf; unfolded spectrum
** TH2D smear; additional smearing matrix, unfolded = smear*r^-1*m
** TH1D wiener; Wiener regularization
** TH1D bias; intrinsic bias, (smear - I)*s/s, s is sig in input.root
** TH1D fracError; fractional uncertainy (corresponding to diagonal elements of unfcov) of the unfolded spectrum
** TH2D unfcov; covariance matrix of the unfolded spectrum uncertainty

Notes:
1. All input in mathematical format (matrix, vector) for Wiener-SVD. However, tools are provided to convert histogram/matrix to matrix/histogram.  
2. Need a choice of addtional matrix in unfolding for smoothness, e.g. 2nd derivative matrix (recommended). An infinitesimal value is possibly needed to add on the diagonal elements and makes the matrix inversible.
3. Current version can handle non-square matrix (only number of rows greater or equal to number of columns is inversible/solvable).
4. Be careful about histogram x/y-axis meaning after coversion. Two options: hist(i+1, j+1) = matrix(i, j) or hist(i+1, j+1)=matrix(j, i). No units considered in the mutual coversion, to be customized by users in the definiation of histograms.
5. The response matrix should include truth event efficiency, otherwise the users should pay attention by themselves. Two ways to include the efficiency:
  5.1. normalize each column j (i.e. truth bin j) of the response matrix and then apply efficiency from truth bin j to each element of this column bin.
  5.2. Construct response matrix R_ij = (number of reco events in bin i from true events in bin j)/(number of true events in bin j)
Other constant factors e.g. flux, target mass, etc. count on the users' discretion. Remember this equation M = R\*S must hold unless you have other purposes or procedures.
6. Output includes the addtional smearing matrix and wiener filter elements after unfolding. These are all info we need to know for this method. 

More details:
Please read the comments in the source files correspondingly.

