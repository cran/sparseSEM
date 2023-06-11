Elastic Net Penalized Maximum Likelihood for Structural Equation Models
================

- [References](#references)

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/sparseSEM)](https://cran.r-project.org/package=sparseSEM)[![](https://cranlogs.r-pkg.org/badges/sparseSEM)](https://CRAN.R-project.org/package=sparseSEM)

We provide extremely efficient procedures for fitting the lasso and
elastic net regularized Structural Equation Models (SEM). The model
output can be used for inferring network structure (topology) and
estimating causal effects. Key features include sparse variable
selection and effect estimation via **l1** and **l2** penalized maximum
likelihood estimator (MLE) implemented with BLAS/Lapack routines. The
implementation enables extremely efficient computation. Details may be
found in Huang A. ([2014](#ref-dissertation)).

Version 3.8 is a major release with several new features, including:

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
Huang A. (2014) <br> Sparse Model Learning for Inferring Genotype and
Phenotype Associations. <br> Ph.D Dissertation, University of Miami,
Coral Gables, FL, USA.
</p>

</div>

</div>
