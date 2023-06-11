\name{enSEM_stability_selection}
\alias{enSEM_stability_selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Stability Selection for the Elastic Net penalized SEM
}
\description{
	Fit the elastic-net penalized structureal Equation Models (SEM) with input data (X, Y): Y = BY + fX + e.
	Perform Stability Selection (STS) on the input dataset. This function implements STS described in
	Meinshausen N. and Buhlmann P (2010) and Shah R. and Samworth R (2013). 
	
	Underlying the function, the program obtains the performs n rounds of boostraping each with half of the original 
	sample size, and run the selection path of hyperparameter (alpha, lambda). The following stability selection scores are 
	calculated: \cr
	1. E(v): the upper bound of the expected number of falsely selected variables  \cr
	2. pre-comparison error rate = E(v)/p where p is the total number of model parameters (in SEM, p = M*M -M) \cr
	3. E(v)_ShaR the expected number of falsely selected variables described in Shah R. and Samworth R (2013) \cr
	4. FDR: False discovery rate = E(v)/nSelected \cr
	5. FDR_ShaR: FDR described in Shah R. and Samworth R (2013) \cr
	
  The final output is based on Scores described in described in Shah R. and Samworth R (2013), and original scores
  described in Meinshausen N. and Buhlmann P (2010) are provided for reference.
  
	}
\usage{
enSEM_stability_selection(Y,X, Missing,B,
                           alpha_factors, 
                           lambda_factors, 
                           kFold,
                           nBootstrap,
                           verbose)
                                      
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
  k-fold cross validation, default k=3.  Note STS result is not based on CV. However, fitting l1/l2 regularized SEM will
  run the first step described in elasticNetSEM() function: 
  Step 1. SEM-ridge regression (L2 penalty) with k-fold CV: this step find the optimal ridge hyperparameter rho to provide an initial values for l1/l2 regularized SEM. \cr
}
  \item{nBootstrap}{
  bootstrapping parameter. default nBootstrap = 100.
}
  \item{verbose}{
	describe the information output from -1 - 10, larger number means more output
}
}
\details{
	the function perform STS
}
\value{

	\item{STS}{
	The stable effects are those effects selected by STS, i.e., the non-zero values in matrix B.
	
}
	\item{statistics}{
	the final STS scores with components of: \cr
                1. threshold: denoted as pi in  Meinshausen N. and Buhlmann P (2010) \cr
                2. pre-comparison error rate   \cr
                3. E(v)   \cr
                4. E(v)_ShahR \cr
                5. nSTS: final number of stable effects with pi that leads to minimum FDR  \cr
                6. FDR \cr
                7. FDR_ShahR \cr
                
	
	}
	\item{STS data}{Bootstrapping details.}
  \item{call}{the call that produced this object}	
}
\references{
[1]: Meinshausen, N. and Buhlmann, P., 2010. Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), pp.417-473.

[2] Shah, R.D. and Samworth, R.J., 2013. Variable selection with error control: another look at stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 75(1), pp.55-80.
}
\author{Anhui Huang}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
	library(sparseSEM)
	data(B);
	data(Y);
	data(X);
	data(Missing);
	#Example
	\donttest{
output = enSEM_stability_selection(Y,X, Missing,B,
                                     alpha_factors = seq(1,0.05, -0.05), 
                                     lambda_factors =10^seq(-0.2,-4,-0.2), 
                                     kFold = 3,
                                     nBootstrap = 100,
                                     verbose = -1)
  }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Elastic_Net}
\keyword{sparseSEM}
\keyword{stability_selection}
