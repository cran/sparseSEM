\name{X}
\alias{X}
\docType{data}
\title{
	Genotype matrix
}
\description{
	X is the M by N matrix corresponding to M cis-eQTL for M genes
}
\usage{data(X)}
\format{
  The format is:
 int [1:30, 1:200] 2 1 3 1 2 3 2 1 2 2 ...
}
\details{
	current implementation only consider 1 eQTL per gene
}

\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068.
2. Huang A., Bazerque J. A., Giannakis G. B., DUroso G., Myers C. L., Cai X., Elastic Net algorithm for inferring gene regulatory networks based on structural equation models, to be submitted.
}
\examples{
data(X)

}
\keyword{datasets}
