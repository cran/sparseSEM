Elastic Net Penalized Maximum Likelihood for Structural Equation Models
with Netowrk GPT Framework
================

- [References](#references)

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/sparseSEM)](https://cran.r-project.org/package=sparseSEM)[![](https://cranlogs.r-pkg.org/badges/sparseSEM)](https://CRAN.R-project.org/package=sparseSEM)

We provide extremely efficient procedures for fitting the lasso and
elastic net regularized Structural Equation Models (SEM). The model
output can be used for inferring network structure (topology) and
estimating causal effects. Key features include sparse variable
selection and effect estimation via **l1** and **l2** penalized maximum
likelihood estimator (MLE) implemented with BLAS/Lapack routines. The
implementation enables extremely efficient computation. Details can be
found in Huang A. ([2014](#ref-dissertation)).

To achieve high performance accuracy, the software implements a Network
Generative Pre-traning Transformer (GPT) framework:

- Perform a `Network GPT` that generates a complete (fully connected)
  graph from **l2** penalized SEM (i.e., ridge SEM); and
- Use the complete graph as the initial state and fit the `elastic net`
  (**l1** and **l2**) penalized SEM.

Note that the term `Transformer` does not carry the same meaning as the
`transformer architecture` commonly used in Natural Language Processing
(NLP). In Network GPT, the term refers to the creation and generation of
the complete graph.

Version 4.0:

- Enhanced documentation with a new vignette
  `Network Inferrence via sparseSEM` to enable quick setup and running
  of the package;
- Added a new `yeast GRN` real dataset that was used to generate the
  graph in the vignettes;
- Added the dataset preprocessing description in the vignette; and
- further streamline function input and output from both C/C++ and R
  functions

Version 3.8:

- simplified user interface with central functions and simple parameters
  setup;
- stability selection function with both serial and parallel
  bootstrapping;
- streamlined function output.

Version 3 is a major release that updates BLAS/Lapack routines according
to R-API change.

## References

<div id="refs" class="references">

<div id="ref-package">

<div id="ref-dissertation">

<p>
Huang Anhui. (2014) <br> Sparse Model Learning for Inferring Genotype
and Phenotype Associations. <br> Ph.D Dissertation, University of Miami,
Coral Gables, FL, USA.
</p>

</div>

</div>
