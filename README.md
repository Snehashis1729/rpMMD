# rpMMD

`rpMMD` is an R package in development for random-projection maximum mean discrepancy (MMD) methods for multivariate two-sample problems. The intended input will be two numeric matrices, with observations in rows and variables in columns. The package will generate random projection directions, project the samples to one dimension, compute MMD statistics, combine information across projections, and use permutation calibration for testing.

## Installation

The development version will be available from GitHub with:

```r
install.packages("remotes")
remotes::install_github("Snehashis1729/rpMMD")
```

## Current development status

The package skeleton, DESCRIPTION file, project structure, and initial documentation are currently being set up. For the remainder of the semester, the project will add random projection utilities, MMD calculations, the main testing workflow, permutation-based p-values, unit tests, documentation, a vignette, and computational improvements using Rcpp.

## Reference

Shi, X., Du, J., and Xie, T. (2025). *Multivariate two-sample tests via random projection*. Statistics and Computing, 35, Article 98.