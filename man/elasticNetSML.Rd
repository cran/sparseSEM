\name{elasticNetSML}
\alias{elasticNetSML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The Elastic Net penalty for SEM
}
\description{
	Fit the elastic Net penalized SEM model with input data X, Y: Y = BY + fX + e.
	
	For users new to this package, elasticNetSML provides the simplified entry point: 
	Missing matrix can be all 0 (none or uknown), B matrix can be all 0 (unknow connections in the network), 
	thus only Y, X needed.
	
	Underlying the function, the function obtained the optimal hyperparameter (alpha, lambda) from 
	k-fold cross validation (CV) with fixed k= 5. and the grid search of candiate alpha and number of 
	steps for candate lambda are also fixed.  Specifically, for each alpha from 0.95 to 0.05 at
	a step of 0.05, the function perform 5 fold CV for lambda_max to lambda_min in 20 step 
	to determine the optimal alpha and lambda for the data. 
	
	Once users get the idea of what the package and the function does, the package also
	provide functions with more flexibility to allow users to tune the parameters.  Please
	see function elasticNetSMLcv() and elasticNetSMLpoint().
	
	Note that the function is computationally expensive. Generally, the following steps are followed: \cr
        Step 1. SEM-ridge regression (L2 penalty) with k-fold CV: this step find the optimal ridge hyperparameter rho \cr
        Step 2. fit SEM reidge regression model (L2 penalty) with rho from Step 1, obtain the initial status (non-sparse)
                of network structure (B_ridge); \cr
        Step 3.SEM-elastic net regression with k-fold CV: this step finds the optimal hyperparameter (alpha, lambda) \cr
        Step 4. fit SEM-elastic net model with (alpha, lambda) from Step 3. \cr
        Step 5. result calculation for PD, FDR, provide the output.
        
  
  For large scale network inference, a standalone C/C++ software with openMPI for 
  parallel computation is also available upon request.
  
  If you are interested in implementating the algorithm to integrit with Spark, please reach me out.
	}
\usage{
elasticNetSML(Y, X, Missing, B, Verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	The observed node response data with dimension of M by N. Y is normalized inside the function.
}
  \item{X}{
	The network node attribute matrix with dimension of M by N, with M being the number of nodes,
      and N being the number of samples. Theoretically, X can be L by N matrix, with L being the total
      node attributes. However, in current implementation, each node only allows one and only one attribute.
      If you have more than one attributes for some nodes,  please consider selecting the top one by either
      correlation or principal component methods.
      X is normalized inside the function.
}
  \item{Missing}{
	M by N matrix corresponding to elements of Y. 0 denotes no missing, while 1 denotes missing.
      If a node j in sample i has a missing label (Missing[j,i] = 1), the node response Y[j,i] is set to 0.
}
  \item{B}{
	 For a network with M nodes, B is the M by M adjacency matrix.
      If data is simulated/with known true network topology (i.e., known adjacency matrix), the Power
      of detection (PD) and False Discovery Rate (FDR) is computed in the output parameter 'statistics'.
  
      If the true network topology is unknown, B is optional, and the PD/FDR in output parameter
      'statistics' should be ignored.
}
  \item{Verbose}{
	describe the information output from 0 - 10, larger number means more output
}
}
\details{
	the function perform CV and parameter inference, calculate power and FDR
}
\value{

	\item{Bout}{
	the computed weights for the network topology. B[i,j] = 0 means there is no edge between node i and j;
        B[i,j]!=0 denotes an (undirected) edge between note i and j.
}
	\item{fout}{
	f is 1 by M array keeping the weight for X (in SEM: Y = BY + FX + e). Theoretically, F can be M by L matrix,
        with M being the number of nodes, and L being the total node attributes. However, in current implementation,
        each node only allows one and only one attribute.
        If you have more than one attributes for some nodes, please consider selecting the top one by either
        correlation or principal component methods.
	
	}
	\item{stat}{
	statistics is 1x6 array keeping record of:
                1. correct positive
                2. total positive
                3. false positive
                4. positive detected
                5. Power of detection (PD) = correct positive/total positive
                6. False Discovery Rate (FDR) = false positive/positive detected
	
	}
	\item{simTime}{computational time}

}
\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068. \cr
2. Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation Chapter 7. University of Miami(1186). \cr
}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}
\note{
		Difference in three functions:\cr
		1) elasticNetSML: Default alpha = 0.95: -0.05: 0.05; default 20 lambdas \cr
		2) elasticNetSMLcv: user supplied alphas (one or more), lambdas; compute the optimal parameters and network parameters \cr
		3) elasticNetSMLpoint: user supplied one alpha and one lambda, compute the network parameters 
		
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
	library(sparseSEM)
	data(B);
	data(Y);
	data(X);
	data(Missing);
	#Example not run for CRAN's CPU time limit requirement.
	\donttest{
	  OUT <- elasticNetSML(Y, X, Missing, B, Verbose = 1); 
  }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Elastic Net}
\keyword{sparseSEM}
