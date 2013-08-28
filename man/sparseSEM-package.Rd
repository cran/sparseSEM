\name{sparseSEM-package}
\alias{sparseSEM-package}
\alias{sparseSEM}
\docType{package}
\title{
	sparseSEM: Sparse-aware Maximum Likelihood for Structural Equation Models in Inferring Gene Regulatory Networks
}
\description{
	state-of-the-art sparse aware maximum likelihood function for structural equation models.
	Two penalty functions including the Lasso and Elastic net available.
	Version 2.3: updated after MPI command line C/C++ program is finalized.
				examples were in comment to reduce CPU time to <5s.
				
}
\details{
\tabular{ll}{
Package: \tab sparseSEM\cr
Type: \tab Package\cr
Version: \tab 2.3\cr
Date: \tab 2014-08-28\cr
License: \tab GPL\cr
}

}
\author{
Anhui Huang

Maintainer: Anhui Huang <a.huang1@umiami.edu>

}
\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068.
2. Huang A., Bazerque J. A., Giannakis G. B., DUroso G., Myers C. L., Cai X., Elastic Net algorithm for inferring gene regulatory networks based on structural equation models, to be submitted.
}
\keyword{sparse aware maximum likelihood, structural equation models, Lasso, Elastic Net}

\examples{
library(sparseSEM)
}
