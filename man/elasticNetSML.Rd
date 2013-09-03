\name{elasticNetSML}
\alias{elasticNetSML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The Elastic Net penalty for SEM
}
\description{
	For each alpha from 0.95 to 0.05 at a step of 0.05, the function perform 5 fold CV for 
	lambda_max to lambda_min in 20 step to determine the optimal alpha and lambda for the data. }
\usage{
elasticNetSML(Y, X, Missing, B, Verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	gene expression M by N matrix
}
  \item{X}{
	cis_eQTL M by N matrix
}
  \item{Missing}{
	missing data in Y
}
  \item{B}{
	true network topology if available
}
  \item{Verbose}{
	describe the information output from 0 - 10, larger number means more output
}
}
\details{
	the function perform CV and parameter inference, calculate power and FDR
}
\value{

	\item{Bout}{the matrix B from sparseSEM}
	\item{fout}{f: the weight for matrix X}
	\item{stat}{compute the power and FDR statistics if the ture topology is provided}
	\item{simTime}{computational time}

}
\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068. \cr
2. Huang A., Bazerque J. A., Giannakis G. B., DUroso G., Myers C. L., Cai X., Elastic Net algorithm for inferring gene regulatory networks based on structural equation models, to be submitted. \cr
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
	#Example is commented to satisfy CRAN's CPU time limit requirement.
	#OUT <- elasticNetSML(Y, X, Missing, B, Verbose = 1); 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Elastic Net}
\keyword{sparseSEM}
