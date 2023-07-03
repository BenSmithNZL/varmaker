#' Simulate observations from a VAR.
#'
#' Creates observations from the canonical VAR. \deqn{z_t = \phi_0 + \phi_1 z_{t-1} + \ldots + \phi_h z_{t-h} + a_t}
#' where \eqn{a_t \sim N(0, \Sigma_a)}
#'
#' @param coefficient_matrices A list of the coefficient matrices. The first element of the list corresponds to
#' \eqn{\phi_0}, and must be a K-length vector. All following elements must be K-by-K matrices, with the second element corresponding to
#' \eqn{\phi_1}, the third to \eqn{\phi_2}, and so on. Only \eqn{\phi_0} is required; all further coefficient matrices are optional.
#' @param Sigma_a A K-by-K covariance matrix of the innovations. Must be symmetrical.
#' @param n An integer of the number of periods.
#' @param seed Seed for the random generation of the innovations. Seed is equal to 56143868 by default.
#'
#' @return \code{create_var} returns a [list()] object with the following elements:
#' * \code{a} is a n-by-K matrix of the innovations.
#' * \code{autocovariance} is a list of the autocovariances.
#' * \code{beta} is the matrix of coefficients.
#' * \code{burn_ins} is the number of burn-in observations used.
#' * \code{companion_matrix} is the coefficients in companion matrix form.
#' * \code{eigens} is list containing the eigenvalues and the moduli of the companion matrix.
#' * \code{H} is an integer for the number of lags in the VAR.
#' * \code{K} is an integer for the number of time series.
#' * \code{mu} is K-vector of the theoretical mean of the process.
#' * \code{n} is the number of periods.
#' * \code{stable} logical. If \code{TRUE} then the process is stable.
#' * \code{z} is the process itself.
#'
#' @export
#'
#' @details
#' There are a certain number of burn-ins that are used.
#'
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
  max_autocovariance = 12

  #### Checks ####
  # Check if seed is an integer
  if (is.numeric(seed) == FALSE | length(seed) != 1) {
    stop("seed must be an integer.")
  } else {
    set.seed(seed)
  }

  # Check if n is an integer
  if (is.numeric(n) == FALSE | length(n) != 1) {
    stop("n must be an integer.")
  }

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

  #### Companion form coefficients ####
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

  #### Eigenvalues and stability ####
  eigenvalues <- eigen(companion_matrix)$values
  moduli <- sqrt(Re(eigenvalues)^2 + Im(eigenvalues)^2)
  if (any(sqrt(Re(eigenvalues)^2 + Im(eigenvalues)^2) >= 1)) {
    warning("Process is unstable.")
    stable <- FALSE
  } else {
    stable <- TRUE
  }

  #### Simulate observations ####
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

  #### Beta coefficient matrix ####
  beta = t(coefficient_matrices[[1]])

  if (H > 0) {

    for (h in 1:H) {

      beta = rbind(beta, coefficient_matrices[[h+1]])

    }

  }

  #### Expected value ####
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



  #### Autocovariance ####

  if (H == 0) {

    autocovariance = NULL

  } else {

    Sigma_A = matrix(data = 0,
                     nrow = K*H,
                     ncol = K*H)

    Sigma_A[1:K, 1:K] <- Sigma_a

    Gamma_Z_0 = solve(diag(x = 1, nrow = (K*H)^2, ncol = (K*H)^2) - kronecker(companion_matrix, companion_matrix)) %*% c(Sigma_A)
    Gamma_Z_0 = matrix(Gamma_Z_0, nrow = K*H, ncol = K*H, byrow = TRUE)

    autocovariance = list()

    for (h in 1:H) {

      autocovariance[[h]] <- Gamma_Z_0[1:K, (K*(h-1)+1):(K*h)]

    }

    # Recursively finding autocovariances
    for (g in (H+1):(max_autocovariance+1)) {

      autocovariance[[g]] <- matrix(data = 0, nrow = K, ncol = K)

      for (h in 1:H) {

        autocovariance[[g]] <- autocovariance[[g]] + coefficient_matrices[[h+1]] %*% autocovariance[[g-h]]

      }

    }

    names(autocovariance) <- 0:max_autocovariance

  }

  #### Return output ####
  return(list(a = a,
              autocovariance = autocovariance,
              beta = beta,
              burn_ins = burn_ins,
              companion_matrix = companion_matrix,
              eigens = list(values = eigenvalues,
                            moduli = moduli),
              H = H,
              K = K,
              mu = mu,
              n = n,
              stable = stable,
              z = z))

}
