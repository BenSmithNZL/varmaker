#' Simulate observations from a VAR.
#'
#' Creates observations from the canonical VAR.
#'
#' @param coefficient_matrices A list of the coefficient matrices. The first element of the list must be a K-length vector, all following elements must be K-by-K matrices.
#' @param Sigma_a A K-by-K covariance matrix of the innovations.
#' @param n An integer of the number of periods.
#' @param seed Seed for the random generation of the innovations. Seed is equal to 56143868 by default.
#'
#' @return \code{create_var} returns a [list()] object.
#' * \code{a} is a matrix of the innovations.
#' * \code{autocovariance} is a list of the autocovariances.
#' * \code{beta} is the matrix of coefficients.
#' * \code{companion_matrix} is the coefficients in companion matrix form.
#' * \code{H} is an integer for the number of lags in the VAR.
#' * \code{K} is an integer for the number of time series.
#' * \code{mu} is the theoretical mean of the process.
#' * \code{n} is the number of periods.
#' * \code{z} is the process itself.
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
#' @references Lütkepohl, H. (2005) \emph{New Introduction to Multiple Time Series Analysis.} Springer.
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
  for (h in 1:H) {
    if (isFALSE(all(dim(coefficient_matrices[[h+1]]) == K))) {
      stop(paste0("phi_", h, " is not a ", K, " by ", K, ' matrix.',  sep = ""))
    }
  }

  # Check if sigma is a K by K matrix
  if (isFALSE(all(dim(Sigma_a) == K))) {
    stop("Sigma is not a K by K matrix.")
  }

  # Check Sigma_a is symmetric
  if (isSymmetric(Sigma_a) == FALSE) {
    stop("Sigma_a is not a symmetric matrix")
  }

  ## Calculate the VAR
  burn_ins = min(round(n/2), 500)
  a = MASS::mvrnorm(n = n + burn_ins, mu = rep(0, K), Sigma_a)
  z = matrix(0, nrow = n + burn_ins, ncol = K)

  for (t in (H+1):(n + burn_ins)) {

    z[t, ] <- t(coefficient_matrices[[1]])

    for (h in 1:H) {

      z[t, ] <- z[t, ] + coefficient_matrices[[h+1]] %*% z[t-h, ]

    }

    z[t, ] <- z[t, ]+  a[t, ]

  }

  a <- utils::tail(a, n)
  z <- utils::tail(z, n)
  rownames(a) <- 1:n
  rownames(z) <- 1:n

  ## Creating the coefficient matrix beta
  beta = t(coefficient_matrices[[1]])

  if (H > 0) {

    for (h in 1:H) {

      beta = rbind(beta, coefficient_matrices[[h+1]])

    }

  }

  ## Calculating the expected value of the DGP.
  if (H == 0) {

    mu = as.matrix(coefficient_matrices[[1]])

  } else {

    mu = diag(x = 1, nrow = K, ncol = K)

    for (h in 1:H) {

      mu = mu - coefficient_matrices[[h+1]]

    }

    mu = solve(mu)
    mu = mu %*% coefficient_matrices[[1]]

  }

  ## Calculating companion-form coefficient matrix.
  if (H == 0) {

    companion_matrix = NULL

  } else if (H == 1){

    companion_matrix = coefficient_matrices[[2]]

  } else {

    companion_matrix = matrix(data = 0,
                              nrow = K*H,
                              ncol = K*H)

    for (h in 1:H) {

      companion_matrix[1:K, (K*(h-1)+1):(K*h)] <- coefficient_matrices[[h+1]]

    }

    companion_matrix[cbind(c((K+1):(K*H)), c(1:(K*(H-1))))] <- 1

  }

  ## Calculating autocovariances
  Sigma_A = matrix(data = 0,
                   nrow = K*2,
                   ncol = K*2)

  Sigma_A[1:K, 1:K] <- Sigma_a

  Gamma_Z_0 = solve(diag(x = 1, nrow = (K*H)^2, ncol = (K*H)^2) - kronecker(companion_matrix, companion_matrix)) %*% c(Sigma_A)
  Gamma_Z_0 = matrix(Gamma_Z_0, nrow = K*H, ncol = K*H, byrow = TRUE)

  autocovariance = list(Gamma_Z_0[1:K, 1:K],
                        Gamma_Z_0[1:K, (K+1):(K*H)])



  ## Return output
  return(list(a = a,
              autocovariance = autocovariance,
              beta = beta,
              burn_ins = burn_ins,
              companion_matrix = companion_matrix,
              H = H,
              K = K,
              mu = mu,
              n = n,
              z = z))

}

