#' Create a multivariate time series from a VAR
#'
#' @param coefficient_matrices A list of the coefficient matrices. The first element of the list must be a K-length vector, all following elements must be K-by-K matrices.
#' @param Sigma_a A K-by-K covariance matrix of the innovations.
#' @param n An integer of the number of periods.
#' @param seed Seed for the random generation of the innovations. Seed is equal to 56143868 by default.
#'
#' @return \code{create_var} returns a [list()] object.
#'
#' @export
#'
#' @examples
#' create_var(list(c(1, 0.5, 0.25),
#'                 matrix(c(0.2, 0, 0.1,
#'                         -0.3, 0, 0,
#'                          0, 0, -0.1),
#'                          nrow = 3,
#'                          ncol = 3,
#'                          byrow = TRUE),
#'                 matrix(c(0, 0.3, 0,
#'                          0, 0, 0,
#'                          0, 0, 0.1),
#'                          nrow = 3,
#'                          ncol = 3,
#'                          byrow = TRUE)),
#'            matrix(c(1.25, 0, 0,
#'                     0, 1, 0,
#'                     0, 0, 0.5),
#'                     nrow = 3,
#'                     ncol = 3,
#'                     byrow = TRUE),
#'            1000)
create_var <- function(coefficient_matrices, Sigma_a, n, seed = 56143868){

  K = length(coefficient_matrices[[1]])
  H = length(coefficient_matrices) - 1
  set.seed(seed)

  ## Check the dimensions of the coefficient matrices

  # Check that phi_0 is not a scalar.
  if (length(coefficient_matrices[[1]]) == 1) {
    stop("phi_0 is not a vector.")
  }

  # Check if phi_0 is a column vector
  if (is.vector(coefficient_matrices[[1]]) == FALSE) {
    stop("phi_0 is not a vector.")
  }

  # Check if they're all a K by K matrix
  for (h in 2:(H+1)) {
    if (isFALSE(all(dim(coefficient_matrices[[h]]) == K))) {
      stop(paste0("phi_", h-1, " is not a ", K, " by ", K, ' matrix.',  sep = ""))
    }
  }

  # Check if sigma is a K by K matrix
  if (isFALSE(all(dim(Sigma_a) == K))) {
    stop("Sigma is not a K by K matrix.")
  }

  ## Calculate the VAR
  burn_ins = min(round(n/2), 500)
  a = MASS::mvrnorm(n = n + burn_ins, mu = rep(0, K), Sigma_a)
  z = matrix(0, nrow = n + burn_ins, ncol = K)

  for (t in (H+1):(n + burn_ins)) {

    z[t, ] <- t(coefficient_matrices[[1]])

    for (h in 1:H) {

      z[t, ] <- z[t, ] + z[t-h, ] %*% coefficient_matrices[[h+1]]

    }

    z[t, ] <- z[t, ]+  a[t, ]

  }

  a <- utils::tail(a, n)
  rownames(a) <- 1:n
  z <- utils::tail(z, n)
  rownames(z) <- 1:n

  ## Creating the coefficient matrix beta
  beta = t(coefficient_matrices[[1]])

  for (h in 1:H) {
    beta = rbind(beta, coefficient_matrices[[h+1]])
  }

  ## Calculating the expected value of the DGP.
  mu = diag(x = 1, nrow = K, ncol = K)
  for (h in 1:H) {
    mu = mu - coefficient_matrices[[h+1]]
  }
  mu = solve(mu)
  mu = mu %*% coefficient_matrices[[1]]

  ## Return output
  return(list(a = a,
              beta = beta,
              burn_ins = burn_ins,
              H = H,
              K = K,
              mu = mu,
              n = n,
              z = z))
}
