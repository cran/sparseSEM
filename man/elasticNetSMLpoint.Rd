\name{elasticNetSMLpoint}
\alias{elasticNetSMLpoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The Elastic Net penalty for SEM
}
\description{
	For user supplied one alpha in range of (0,1) and one lambda, 
	the function perform selection path from lambda_max to lambda 
	to determine the optimal network topology. 
	
	In the case of the grid provided in elasticNetSMLcv() function may not 
	be granular enough and user would like to explore/twist (alpha, lambda)
	a little bit, this function provides the solution.
	
	Note: While in elasticNetSML() and elasticNetSMLcv(), k-fold CV was conducted 
	to determine the optimal (alpha, lambda), elasticNetSMLpoint() does not.
	
	}
\usage{
	elasticNetSMLpoint(Y, X, Missing, B, alpha_factor, lambda_factor, Verbose)
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
\item{alpha_factor}{
	alpha_factor: in range of (0, 1); must be scalar
}
  \item{lambda_factor}{
	penalty lambda_factor: in range of (0, 1); must be scalar
}
  \item{Verbose}{
	describe the information output from 0 - 10, larger number means more output
}
}
\details{
	the function perform selection path from lambda_max to lambda, calculate power and FDR
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
2. Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).\cr \cr	
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
\dontrun{	OUT <- elasticNetSMLpoint(Y, X, Missing, B,
		alpha_factor = 0.5, lambda_factor = 0.1, Verbose = 1);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Elastic Net}
\keyword{sparseSEM}
