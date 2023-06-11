\name{yeast}
\alias{yeast cisQTL GRN}
\docType{data}
\title{
	Yeast cis-QTL Gene Regulatory Network Dataset
}
\description{
	The dataset (Y,X) are two matrices each with 3380 rows and 112 columns as described in the Vignettee "Elastic Net Enabled Sparse-Aware Maximum Likelihood for
Structural Equation Models in Inferring Gene Regulatory Networks". 

  The Yeast expression trait data set was obatained from Brem R.B. and Kruglyak L (2005), and has been screened through:\cr
  1. screen and keep the ORF names in Kellis's ORF list (Kellis M et al, 2003); \cr
  2. screen out ORF with more than 5\% of missing expression date; \cr
  3. perform eQTL mapping by Wilcoxon test, adjust by qvalue; keep the top one cisQTL; \cr
  4. fill in zeros for ORF without cisQTL in matrix X \cr
}
\usage{data(yeast)}
\format{
  The format is matrix.
 
}
\details{
	Yeast cis-QTL Gene Regulatory Network Dataset 
}
\references{
1. Brem RB, Kruglyak L: The landscape of genetic complexity across 5,700 gene expression traits in yeast. Proceedings of the National Academy of Sciences of the United States of America 2005, 102:1572-1577.
 \cr
2. Kellis M, Patterson N, Endrizzi M, Birren B, Lander ES: Sequencing and comparison of yeast species to identify genes and regulatory elements. Nature 2003, 423:241-254\cr \cr	
}
\examples{
data(yeast)
}
\keyword{datasets}
