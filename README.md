# rpMMD

## Overview

rpMMD is an R package for random-projection maximum mean
discrepancy methods for multivariate two-sample testing.
The package generates random projection directions,
projects multivariate observations onto these directions,
computes Gaussian-kernel squared MMD estimates, constructs
PCvM-type and PKS-type combined statistics, and calculates
permutation-based p-values. Computational acceleration
is provided through an Rcpp implementation.

## Installation

Install the development version directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("Snehashis1729/rpMMD")
```

A working C++ compiler is required for source installation.

To install the package with its vignette:

```r
remotes::install_github(
  "Snehashis1729/rpMMD",
  build_vignettes = TRUE
)
```

## Example

```r
library(rpMMD)

set.seed(123)

X <- matrix(rnorm(150), ncol = 3)
Y <- matrix(rnorm(150, mean = 0.5), ncol = 3)

fit <- rpmmd_test(
  X, Y,
  statistic = "cvm",
  nproj = 20,
  nperm = 199,
  engine = "cpp"
)

print(fit)
summary(fit)
plot(fit)
```

## Available functions

| Function | Description |
|----------|-------------|
| `rpmmd_directions()` | Generate random unit directions |
| `rpmmd_project()` | Project multivariate observations |
| `rpmmd_mmd()` | Compute squared Gaussian-kernel MMD |
| `rpmmd_test()` | Perform random-projection testing |
| `print()` | Display test results |
| `summary()` | Summarize test results |
| `plot()` | Visualize the permutation distribution |

## Vignette

```r
vignette("introduction", package = "rpMMD")
```

## Current development status

- [x] R package structure
- [x] Random projection directions
- [x] Multivariate sample projection
- [x] Gaussian-kernel squared MMD estimation
- [x] PCvM-type and PKS-type statistics
- [x] Permutation-based p-values
- [x] Rcpp computational engine
- [x] Permutation computation optimization
- [x] S3 constructor and validator
- [x] S3 print, summary and plot methods
- [x] Initial unit tests
- [x] Initial vignette source
- [ ] Systematic statistical simulation studies
- [ ] Computational benchmarking
- [ ] Expanded testing and numerical validation
- [ ] Final package validation

## Remaining work

The remainder of the project will focus on empirical Type I
error and power studies, computational benchmarking of the
R and C++ engines, further optimization, verification against
the motivating paper, expanded unit testing, and final
documentation. The package is currently an intermediate
development release.

## Reference

Shi, X., Du, J., and Xie, T. (2025).
*Multivariate two-sample tests via random projection*.
Statistics and Computing, 35, Article 98.

https://doi.org/10.1007/s11222-025-10627-5
