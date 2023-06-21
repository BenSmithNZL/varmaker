#' Create a multivariate time series from a VAR
#'
#' @param coefficient_matrices A list of the coefficient matrices
#' @param Sigma_a Covariance matrix of the innovations
#' @param n Number of periods
#' @param seed Seed for the random generation of the innovations
#'
#' @return A list
#' @export
#'
#' @examples
#' create_var(list(matrix(c(1, 0.5, 0.25),
#'                          nrow = 3,
#'                          ncol = 1,
#'                          byrow = TRUE),
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
#'                 diag(x = 1, nrow = 3, ncol = 3),
#'                 1000)
create_var <- function(coefficient_matrices, Sigma_a, n, seed = 56143868){

  K = dim(coefficient_matrices[[1]])[1]
  H = length(coefficient_matrices) - 1
  set.seed(seed)

  ## Check the dimensions of the coefficient matrices

  # Check if phi_0 is a column vector
  if (dim(coefficient_matrices[[1]])[2] != 1) {
    stop("phi_0 is not a column vector.")
  }

  # Check if they're all a K by K matrix
  for (h in 2:(H+1)) {
    if (isFALSE(all(dim(coefficient_matrices[[h]]) == K))) {
      stop(paste0("phi_", h-1, " is not a ", K, " by ", K, ' matrix.',  sep = ""))
    }
  }

  # Check if sigma is a K by K matrix
  if (dim(Sigma_a)[1] != K) {
    stop("Sigma is not a K by K matrix.")
  }
  if (dim(Sigma_a)[2] != K) {
    stop("Sigma is not a K by K matrix.")
  }

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

  a <- tail(a, n)
  z <- tail(z, n)

  return(list(K = K,
              H = H,
              a = a,
              z = z))
}
