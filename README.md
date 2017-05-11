Implementaion of Wiener-SVD unfolding    
Reference: [arXiv:1705.03568](https://arxiv.org/abs/1705.03568)]  
Version 2.0  

Pre-requiste of ROOT installed, some ROOT classes (e.g. TMatrix) will be used.  
Core algorithms implementation in src/  
Header files in include/  
Use Makefile to build excutables.  

An example of application in Example.C  
Input file for Example.C in Test/  
Usage:  
make clean;  
make;  
Example Test/Test.root output.root 2  


INSTRUCTION

Input: 
1. true/expected signal spectrum (in reality predicted by models); 
2. measured spectrum; 
3. covariance matrix of measured spectrum; 
4. response matrix from signal (corresponding to matrix column) to measurement (corresponding to matrix row), M = R*S, column normalized (pdf);
5. Choice of addional matrix for smoothness, e.g. 2nd derivative matrix.  

Important notes:
1. All input in mathematical format (matrix, vector) for Wiener-SVD. However, tools are provided to convert histogram/matrix to matrix/histogram. 
2. Need a choice of addtional matrix in unfolding for smoothness, e.g. 2nd derivative matrix (recommended). An infinitesimal value is possibly needed to add on the diagonal elements and makes the matrix inversible.
3. Current version can handle non-square matrix (only number of rows greater or equal to number of columns is inversible/solvable).
4. Be careful about histogram x/y-axis meaning after coversion. Two options: hist(i+1, j+1) = matrix(i, j) or hist(i+1, j+1)=matrix(j, i). No units considered in the mutual coversion, to be customized by users in the definiation of histograms.
5. Output includes the addtional smearing matrix and wiener filter elements after unfolding (key information defined in reference [to be added]). Easy to perform an error propogation with it as well as other applications. 

More details:
Please read the comments in the source files correspondingly.

