\name{B}
\alias{B}
\docType{data}
\title{
	True network edges
}
\description{
	B is the M by M matrix defining network topology
}
\usage{data(B)}
\format{
	The format is: M by M, where M is the number of vertices (genes)
	num [1:30, 1:30] 0 0 0 0 0 ...
}
\details{
	If B is not available (real data): the stat output that describes the true accuracy and FDR should be ignored.
}

\references{
1. Cai, X., Bazerque, J.A., and Giannakis, G.B. (2013). Inference of Gene Regulatory Networks with Sparse Structural Equation Models Exploiting Genetic Perturbations. PLoS Comput Biol 9, e1003068.
2. Huang A., Bazerque J. A., Giannakis G. B., DUroso G., Myers C. L., Cai X., Elastic Net algorithm for inferring gene regulatory networks based on structural equation models, to be submitted.
}
\examples{
data(B)
}
\keyword{datasets}
