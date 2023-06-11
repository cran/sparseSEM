\name{elasticNetSEMcv}
\alias{elasticNetSEMcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The Elastic Net penalty for SEM with user supplied alphas and lambdas
}
\description{
Function elasticNetSEMcv allows users to set their own grid search through combination of a set of user provided
	alpha an lambda.
 }
\usage{
elasticNetSEMcv(Y, X, Missing, B, alpha_factors,lambda_factors,kFold, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	The observed node response data with dimension of M (nodes) by N (samples). Y is normalized inside the function.
}
  \item{X}{
	The network node attribute matrix with dimension of M by N. Theoretically, X can be L by N matrix, with L being the total
      node attributes. In current implementation, each node only allows one and only one attribute. \cr
      If you have more than one attributes for some nodes,  please consider selecting the top one by either
      correlation or principal component methods.  \cr
      If for some nodes there is no attribute available, fill in the rows with all zeros.  See the yeast data `yeast.rda` for example. \cr
      X is normalized inside the function.
}
  \item{Missing}{
	Optional M by N matrix corresponding to elements of Y. 0 denotes not missing, and 1 denotes missing.
      If a node i in sample j has the label missing (Missing[i,j] = 1), then Y[i,j] is set to 0.
}
  \item{B}{
	 Optional input. For a network with M nodes, B is the M by M adjacency matrix.
      If data is simulated/with known true network topology (i.e., known adjacency matrix), the Power
      of detection (PD) and False Discovery Rate (FDR) is computed in the output parameter 'statistics'.
  
      If the true network topology is unknown, B is optional, and the PD/FDR in output parameter
      'statistics' should be ignored.
}
  \item{alpha_factors}{
	The set of candidate alpha values.  Default is seq(start = 0.95, to = 0.05, step = -0.05)
}
  \item{lambda_factors}{
	The set of candidate lambda values. Default is 10^seq(start =1, to = 0.001, step = -0.2)
}
  \item{kFold}{
  k-fold cross validation, default k=5
}
  \item{verbose}{
	describe the information output from -1 - 10, larger number means more output
}

}
\details{
	the function perform CV and parameter inference, calculate power and FDR
}
\value{
	\item{cv}{
	dataframe stores the minimum Mean Square Error (MSE) for each alpha and the corresponding lambda from the selection path [lambda_max, ...., lambda_min]. \cr
	col1: alpha \cr
	col2: lambda (With the given alpha, this is the lambda having minimum MSE) \cr
	col3: MSE \cr
	col4: STE \cr
	
	The final (alpha, lambda) is set at the (alpha, lambda) that is within 1ste of the min(MSE) with higher level of penalty on the likehood function. 
	
  }

	\itemize{
    \item{fit}{the model fit with optimal (alpha,lambda) from cv}
    \itemize{ 
	    \item{Bout}{
  	      the computed weights for the network topology. B[i,j] = 0 means there is no edge between node i and j;
          B[i,j]!=0 denotes an (undirected) edge between note i and j with B[i,j] being the weight of the edge.
          }
      \item{fout}{
      	f is 1 by M array keeping the weight for X (in SEM: Y = BY + FX + e). Theoretically, F can be M by L matrix,
              with M being the number of nodes, and L being the total node attributes. However, in current implementation,
              each node only allows one and only one attribute.
              If you have more than one attributes for some nodes, please consider selecting the top one by either
              correlation or principal component methods.
      	
      	}
    	\item{stat}{
    	statistics is 1x6 array keeping record of: \cr
                    1. correct positive \cr
                    2. total positive   \cr
                    3. false positive   \cr
                    4. positive detected  \cr
                    5. Power of detection (PD) = correct positive/total positive  \cr
                    6. False Discovery Rate (FDR) = false positive/positive detected
    	
    	}
	  \item{simTime}{computational time}
    \item{call}{the call that produced this object}	
    }
  }
}
\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068. \cr
2. Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).\cr \cr	 \cr
}
\author{Anhui Huang; Dept of Electrical and Computer Engineering, Univ of Miami, Coral Gables, FL}
\note{
		Difference in three functions:\cr
		1) elasticNetSML: Default alpha = 0.95: -0.05: 0.05; default 20 lambdas \cr
		2) elasticNetSEMcv: user supplied alphas (one or more), lambdas; compute the optimal parameters and network parameters \cr
		3) elasticNetSMLpoint: user supplied one alpha and one lambda, compute the network parameters 
		
		User is responsible to set the random seed to guarantee repeatable results.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
	library(sparseSEM)
	data(B);
	data(Y);
	data(X);
	data(Missing);
	\dontrun{OUT <- elasticNetSEMcv(Y, X, Missing, B, alpha_factors = c(0.75, 0.5, 0.25),
	lambda_factors=c(0.1, 0.01, 0.001), kFold = 5, verbose  = 1);
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Elastic_Net}
\keyword{sparseSEM}
