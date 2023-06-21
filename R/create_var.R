#' Create a multivariate time series from a VAR
#'
#' @param x Some sort of input
#'
#' @return A character string
#' @export
#'
#' @examples
#' create_var("Hello")
create_var <-function(coefficient_matrices, Sigma_a, n, seed = 56143868){

  K = dim(coefficient_matrices[[1]])[1]
  H = length(coefficient_matrices) - 1

  ## Check the dimensions of the coefficient matrices

  # Check if phi_0 is a column vector
  if (dim(coefficient_matrices[[1]])[2] != 1) {
    stop("phi_0 is not a column vector.")
  }

  # Check if they're all a K by K matrix
  for (h in 2:(H+1)) {
    if (isFALSE(dim(coefficient_matrices[[1]]) == K)) {
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

  return(list(K = K,
              H = H,
              a = a,
              z = z))
}
