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

A VAR(2) with two series takes the form: $$ z_t = \phi_0 + \phi_1 z_{t_1}+ \phi_2 z_{t_2} + a_t$$with $a_t \sim N_2(0, \Sigma_a)$. As an example, we can set $\phi_0 = \begin{bmatrix} 2 \\ 0 \\ \end{bmatrix}, \phi_1 = \begin{bmatrix} 0.5 & 0.1 \\ 0.4 & 0.5 \\ \end{bmatrix}, \phi_2 = \begin{bmatrix} 0 & 0 \\ 0.25 & 0 \\ \end{bmatrix}$ and $\Sigma_a = \begin{bmatrix} 0.09 & 0 \\ 0 & 0.04 \\ \end{bmatrix}$. To simulate observations from this model, you can run:

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

A VMA(2) with two series takes the form: $$ z_t = \theta_0 + \theta_1 a_{t-1} + \theta_2 a_{t-2} + a_t$$with $a_t \sim N_2(0, \Sigma_a)$. As an example, we can set $\theta_0 = \begin{bmatrix} 2 \\ 0 \\ \end{bmatrix}, \theta_1 = \begin{bmatrix} 0.5 & 0.1 \\ 0.4 & 0.5 \\ \end{bmatrix}, \theta_2 = \begin{bmatrix} 0 & 0 \\ 0.25 & 0 \\ \end{bmatrix}$ and $\Sigma_a = \begin{bmatrix} 0.09 & 0 \\ 0 & 0.04 \\ \end{bmatrix}$. To simulate observations from this model, you can run:

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
