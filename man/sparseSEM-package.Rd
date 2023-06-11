\name{sparseSEM-package}
\alias{sparseSEM-package}
\alias{sparseSEM}
\docType{package}
\title{
	sparseSEM: Sparse-Aware Maximum Likelihood for Structural Equation Models in Inferring Gene Regulatory Networks
}
\description{
	state-of-the-art sparse aware maximum likelihood function for structural equation models.
	Two penalty functions including the Lasso and Elastic-net available. \cr
	For users new to this package, function elasticNetSML() provides the simplified entry point: 
	Missing matrix can be all 0 (none or uknown), so as B matrix (unknow connections in the network), 
	thus only Y and X are mandatory. \cr
	
		Once users get the idea of what the package and the function does, the package also
	provides other functions with more flexibility to allow fine tuning the parameters: \cr
	Function elasticNetSMLcv() user provides alphas (one or more) and lambdas; the function then computes the optimal parameters and network parameters; \cr
	Function elasticNetSMLpoint(). user provides one alpha and one lambda, and the function computes the network parameters.
				
	Note:
	 Regularization path: the lasso-strong rule is applied, thus the path is run through lambda_max to lambda_opt,
        where lambda_max is the lambda that keeps only 1 non-zero edge, and lambda_opt is the optimal lambda from CV. \cr
        
				
}
\details{
\tabular{ll}{
Package: \tab sparseSEM\cr
Type: \tab Package\cr
Version: \tab 3.8-2\cr
Date: \tab 2023-04-17\cr
License: \tab GPL\cr
}

}
\author{
Anhui Huang

Maintainer: Anhui Huang <anhuihuang@gmail.com>

}
\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068. \cr
2. Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).\cr \cr	 \cr
}
\keyword{structural equation models}

\examples{
library(sparseSEM)
}
