\name{lassoSML}
\alias{lassoSML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The Lasso penalty for SML
}
\description{
	Upon lambda_max to lambda_min in 20 step, the function compute 5 fold CV 
	to determine the optimal lambda for the data. 
}
\usage{
	lassoSML(Y, X, Missing, B, Verbose = 5)
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

	\item{Bout}{the matrix B from SEM}
	\item{fout}{f: the weight for matrix X}
	\item{stat}{compute the power and FDR statistics if the ture topology is provided}
	\item{simTime}{computational time}

}
\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068. \cr
2. Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).\cr \cr	
}
\author{Anhui Huang}


\examples{
	library(sparseSEM)
	data(B);
	data(Y);
	data(X);
	data(Missing);
	\dontrun{OUT <- lassoSML(Y, X, Missing, B, Verbose = 0); 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Lasso}
\keyword{sparseSEM}
