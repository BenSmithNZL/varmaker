---
---
---

# About

varmaker simulates observations from the vector autoregressive framework and also returns the theoretical properties of those observations.

# Installation

This package has not been released to CRAN yet, but you can install the development version of varmaker with:

``` r
install.packages("devtools")
devtools::install_github("BenSmithNZL/varmaker")
```

# Example

## VAR

To simulate observations from a VAR(2) with two series, you can run:

``` r
library(varmaker)

cm_1 <- list(c(2, 0),
             matrix(c(0.5, 0.1,
                      0.4, 0.5),
                    nrow = 2,
                    ncol = 2,
                    byrow = TRUE),
             matrix(c(0, 0,
                      0.25, 0),
                    nrow = 2,
                    ncol = 2,
                    byrow = TRUE))

Sigma_a_1 <- matrix(c(0.09, 0,
                      0, 0.04),
                    nrow = 2,
                    ncol = 2,
                    byrow = TRUE)

data_1 <- create_var(cm_1, Sigma_a_1, n = 1000)
```

![](https://drive.google.com/uc?id=1dKfcRPyWczK-8sBBlExNAFlL1Cu5TuH5)

The object `data_1` contains the simulated observations themselves, along with theoretical properties of the process such as the mean, autocovariance, autocorrelation, and the Granger-causalities of the series.

## VMA

To simulate observations from a VMA(2) with two series, you can run:

``` r
library(varmaker)

cm_1 <- list(c(2, 0),
             matrix(c(0.5, 0.1,
                      0.4, 0.5),
                    nrow = 2,
                    ncol = 2,
                    byrow = TRUE),
             matrix(c(0, 0,
                      0.25, 0),
                    nrow = 2,
                    ncol = 2,
                    byrow = TRUE))

Sigma_a_1 <- matrix(c(0.09, 0,
                      0, 0.04),
                    nrow = 2,
                    ncol = 2,
                    byrow = TRUE)

data_2 <- create_vma(cm_1, Sigma_a_1, n = 1000)
```

![](https://drive.google.com/uc?id=1Q03Jt7NLV9UqPfDfHJQt0Mn9cu1VqE8U)
