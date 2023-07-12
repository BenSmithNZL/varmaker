#' Simulate observations from a VMA(Q).
#'
#' Simulate observations from the canonical VMA of order Q and dimension K.
#'
#' @param theta A list of the coefficient matrices.
#' @param Sigma_a A K-by-K covariance matrix of the innovations.
#' @param n An integer of the number of periods.
#' @param seed Seed for the random generation of the innovations. Set to 56143868 by default.
#' @param burn_ins The number of burn in observations to use. Set to 500 by default.
#'
#' @return \code{create_var} returns a [list()] object with the following elements:
#' * \code{a} Matrix of the innovations.
#' * \code{autocorrelations} List of the theoretical autocorrelations.
#' * \code{autocovariance} List of the theoretical autocovariances.
#' * \code{beta} Matrix of coefficients.
#' * \code{burn_ins} Integer of the number of burn-in observations used.
#' * \code{companion_matrix} Matrix of coefficients in companion form.
#' * \code{K} Integer for the dimension of the VAR.
#' * \code{mu} Vector of the theoretical mean of the process.
#' * \code{n} Integer for the number of observations.
#' * \code{Q} Integer for the order of the VAR.
#' * \code{theta} List containing the coefficients used.
#' * \code{Sigma_a} Matrix of the covariances of the innovations.
#' * \code{z} Matrix of the simulated observations.
#'
#' @export
#'
#' @details
#' This function simulates observations from the canonical VAM model of order Q and dimension K. The model is given by:
#' \deqn{z_t = \theta_0 + \theta_1 a_{t-1} + \ldots + \theta_Q a_{t-Q} + a_t}
#' where \eqn{z_t} is a K-length vector, \eqn{\theta_0} is a K-length intercept vector, \eqn{\theta_h} is a K-by-K coefficient matrix, \eqn{a_t} is a K-length vector of innovations, and \eqn{a_t \sim N(0, \Sigma_a)}.
#'
#' The first element of \code{theta} corresponds to
#' \eqn{\theta_0} and must be a K-length vector. All following elements must be K-by-K matrices, with the second element corresponding to
#' \eqn{\theta_1}, the third to \eqn{\theta_2}, and so on. Only \eqn{\theta_0} is required; all further coefficient matrices are optional.
#' \code{Sigma_a} must be symmetrical or else an error will be thrown.
#'
#' A certain number of burn-ins observations are used when simulation. These observations are discarded once the observations are generated.
#'
#'
#' @examples
#' create_vma(list(c(1, 0.5, 0.25),
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
#' @references Lütkepohl, Q. (2005) \emph{New Introduction to Multiple Time Series Analysis.} Springer.
create_vma <- function(theta, Sigma_a, n, seed = 56143868, burn_ins = 500){

  K = length(theta[[1]])
  Q = length(theta) - 1

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

  # Check if burn_ins is an integer
  if (is.numeric(n) == FALSE | length(n) != 1) {

    stop("burn_ins must be an integer.")

  }

  # Check that theta_0 is not a scalar.
  if (length(theta[[1]]) == 1) {

    stop("theta_0 is not a vector.")

  }

  # Check if theta_0 is a column vector
  if (is.vector(theta[[1]]) == FALSE) {

    stop("theta_0 is not a vector.")

  }

  # Check if they're all a K by K matrix
  if (Q > 1) {

    for (q in 1:Q) {

      if (isFALSE(all(dim(theta[[q+1]]) == K))) {

        stop(paste0("theta_", q, " is not a ", K, " by ", K, ' matrix.',  sep = ""))

      }

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

  # Warn if there are non-zero off-diagonal elements of Sigma_a
  if (any(Sigma_a[row(Sigma_a)!=col(Sigma_a)] == 0) == FALSE) {

    warning("Sigma_a contains non-zero off-diagonal elements.")

  }

  #### Companion form coefficients ####
  if (Q == 0) {

    companion_matrix = NULL

  } else if (Q == 1){

    companion_matrix = theta[[2]]

  } else {

    companion_matrix = matrix(data = 0,
                              nrow = K*Q,
                              ncol = K*Q)
    for (q in 1:Q) {

      companion_matrix[1:K, (K*(q-1)+1):(K*q)] <- theta[[q+1]]

    }

    companion_matrix[cbind(c((K+1):(K*Q)), c(1:(K*(Q-1))))] <- 1

  }

  #### Simulate observations ####
  a = MASS::mvrnorm(n = n + burn_ins, mu = rep(0, K), Sigma_a)
  z = matrix(0, nrow = n + burn_ins, ncol = K)

  if (Q == 0) {

    for (k in 1:K) {

      z[, k] <- theta[[1]][k] + a[, k]

    }

  } else {

    for (t in (Q+1):(n + burn_ins)) {

      z[t, ] <- t(theta[[1]])

      for (q in 1:Q) {

        z[t, ] <- z[t, ] + theta[[q+1]] %*% a[t-q, ]

      }

      z[t, ] <- z[t, ] + a[t, ]

    }

  }

  a <- utils::tail(a, n)
  z <- utils::tail(z, n)
  rownames(a) <- 1:n
  rownames(z) <- 1:n

  #### Beta coefficient matrix ####
  beta = t(theta[[1]])

  if (Q > 0) {

    for (q in 1:Q) {

      beta = rbind(beta, theta[[q+1]])

    }

  }

  #### Autocovariance ####
  if (Q == 0) {

    autocovariance <- NULL

  } else {

    autocovariance <- list()
    theta_modified <- theta
    theta_modified[[1]] <- diag(1, K, K)

    for (g in 0:Q) {

      autocovariance[[g+1]] <- matrix(0, K, K)

      for (j in g:Q) {

        autocovariance[[g+1]] = autocovariance[[g+1]] + (theta_modified[[j+1]] %*% Sigma_a %*% t(theta_modified[[j-g+1]]))

      }

    }

    names(autocovariance) <- paste0('Gamma_', 0:(length(autocovariance) - 1))

  }

  #### Autocorrelation ####
  autocorrelation <- list()

  D_inverse <- solve(sqrt(diag(diag(autocovariance[[1]]), K, K)))

  for (g in 1:length(autocovariance)) {

    autocorrelation[[g]] <- D_inverse %*% autocovariance[[g]] %*% D_inverse

  }

  names(autocorrelation) <- paste0('R_', 0:(length(autocorrelation) - 1))


  #### Output ####
  return(list(a = a,
              autocorrelation = autocorrelation,
              autocovariance = autocovariance,
              beta = beta,
              burn_ins = burn_ins,
              companion_matrix = companion_matrix,
              K = K,
              mu = as.matrix(theta[[1]]),
              n = n,
              Q = Q,
              theta = theta,
              Sigma_a = Sigma_a,
              z = z))

}
