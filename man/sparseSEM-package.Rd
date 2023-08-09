\name{sparseSEM-package}
\alias{sparseSEM-package}
\alias{sparseSEM}
\docType{package}
\title{
	sparseSEM: Elastic Net Penalized Maximum Likelihood for Structural Equation Models with Network GPT Framework
}
\description{
	State-of-the-art Elastic Net Penalized Maximum Likelihood for Structural Equation Models implemented with a Network Generative Pre-training Transformer (Network GPT).
	Two penalty functions including Lasso and Elastic-net are available.
	
	For users new to this package, function elasticNetSEM() provides the simplified entry point: 
	Missing matrix can be all 0 (none or uknown), so as B matrix (unknow connections in the network), 
	thus only Y and X are mandatory. Then model will fit SEM: Y = BY + fX + e. See the reference for model details. 
	
  The package also provides other functions with more flexibility to allow fine tuning the parameters: \cr
	- elasticNetSEMcv(): user provides alphas (one or more) and lambdas; the function then computes the optimal parameters and network parameters; \cr
	- elasticNetSEMpoint(): user provides one pair of (alpha, lambda), and the function computes the network parameters. \cr
  - enSEM_stability_selection(): stability selection via bootstrapping. \cr
	

# Network Generative Pre-training Transformer (GPT) Framework: \cr 
In all functions the "Network GPT" framework is deployed behind the scene. Specifically, a pre-trained network was built in the following steps:
	 
	 - Step 1. Pre-train the model with ridge (L2 penalty) SEM with k-fold CV: this step find the optimal ridge hyperparameter rho;
        
   - Step 2. Generate a complete graph by fitting the SEM ridge regression model with rho from Step 1, obtain the initial status of a (non-sparse, fully connected)
                complete network structure (B_hat). 
                
	Note that the term "Transformer" does not carry the same meaning as the "transformer architecture" commonly used in Natural Language Processing (NLP). In Network GPT, the term means the creation and generation of the complete graph. \cr
	
# Regularization path: 
A lasso-strong rule is developed for SEM and applied to the selection path. In each step from lambda_max to lambda_min,
        where lambda_max is the lambda that keeps only 1 non-zero edge, and lambda_min is the smallest lambda set arbitrarily (eg., from CV, or 0.001 * lambda_max), the elements in B are pre-set to 0 if they meet with the discarding rule. Those elements will not be computed again in the block coordinate ascent algorithm, resulting in reducing computational costs. See the Vignettes and references for more details. \cr
        
				
}
\details{
\tabular{ll}{
Package: \tab sparseSEM\cr
Type: \tab Package\cr
Version: \tab 4.0\cr
Date: \tab 2023-07-28\cr
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
\concept{structural equation models}
\concept{Network GPT}

\examples{
library(sparseSEM)
}
