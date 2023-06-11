\name{X}
\alias{X}
\docType{data}
\title{
	Genotype matrix
}
\description{
	X is the M by N matrix corresponding to M network nodes from N samples.
}
\usage{data(X)}
\format{
  The format is:
 int [1:30, 1:200] 2 1 3 1 2 3 2 1 2 2 ...
}
\details{
	current implementation only consider 1 inpepedent attribute per node. If users 
	have more than one attributes for some nodes,  please consider selecting the 
	top one by either correlation or principal component methods.
}

\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068. \cr
2. Huang, A. (2014). "Sparse model learning for inferring genotype and phenotype associations." Ph.D Dissertation. University of Miami(1186).\cr \cr	
}
\examples{
data(X)

}
\keyword{datasets}
