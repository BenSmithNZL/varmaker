#' Simulate observations from a VAR(H).
#'
#' Simulate observations from the canonical VAR of order H and dimension K.
#'
#' @param phi A list of the coefficient matrices.
#' @param Sigma_a A K-by-K covariance matrix of the innovations..
#' @param n An integer of the number of periods.
#' @param seed Seed for the random generation of the innovations. Set to 56143868 by default.
#' @param burn_ins The number of burn in observations to use. Set to 500 by default.
#'
#' @return \code{create_var} returns a [list()] object with the following elements:
#' * \code{a} Matrix of the innovations.
#' * \code{autocovariance} List of the theoretical autocovariances.
#' * \code{beta} Matrix of coefficients.
#' * \code{burn_ins} Integer of the number of burn-in observations used.
#' * \code{companion_matrix} Matrix of coefficients in companion form.
#' * \code{eigens} List containing the eigenvalues and the moduli of the companion matrix.
#' * \code{H} Integer for the order of the VAR.
#' * \code{K} Integer for the dimension of the VAR.
#' * \code{mu} Vector of the theoretical mean of the process.
#' * \code{n} Integer for the number of observations.
#' * \code{stable} Logical. If \code{TRUE} then the process is stable.
#' * \code{z} Matrix of the simulated observations.
#'
#' @export
#'
#' @details
#' This function simulates observations from the canonical VAR model of order H and dimension K. The model is given by:
#' \deqn{z_t = \phi_0 + \phi_1 z_{t-1} + \ldots + \phi_H z_{t-H} + a_t}
#' where \eqn{z_t} is a K-length vector, \eqn{\phi_0} is a K-length intercept vector, \eqn{\phi_h} is a K-by-K coefficient matrix, \eqn{a_t} is a K-length vector of innovations, and \eqn{a_t \sim N(0, \Sigma_a)}.
#'
#' The first element of \code{phi} corresponds to
#' \eqn{\phi_0} and must be a K-length vector. All following elements must be K-by-K matrices, with the second element corresponding to
#' \eqn{\phi_1}, the third to \eqn{\phi_2}, and so on. Only \eqn{\phi_0} is required; all further coefficient matrices are optional.
#' \code{Sigma_a} must be symmetrical or else an error will be thrown.
#'
#' A certain number of burn-ins observations are used when simulation. These observations are discarded once the observations are generated.
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
create_var <- function(phi, Sigma_a, n, seed = 56143868, burn_ins = 500){

  K = length(phi[[1]])
  H = length(phi) - 1
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

  # Check if burn_ins is an integer
  if (is.numeric(n) == FALSE | length(n) != 1) {

    stop("burn_ins must be an integer.")

  }

  # Check that phi_0 is not a scalar.
  if (length(phi[[1]]) == 1) {

    stop("phi_0 is not a vector.")

  }

  # Check if phi_0 is a column vector
  if (is.vector(phi[[1]]) == FALSE) {

    stop("phi_0 is not a vector.")

  }

  # Check if they're all a K by K matrix
  if (H > 1) {

    for (h in 1:H) {

      if (isFALSE(all(dim(phi[[h+1]]) == K))) {

        stop(paste0("phi_", h, " is not a ", K, " by ", K, ' matrix.',  sep = ""))

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
  if (H == 0) {

    companion_matrix = NULL

  } else if (H == 1){

    companion_matrix = phi[[2]]

  } else {

    companion_matrix = matrix(data = 0,
                              nrow = K*H,
                              ncol = K*H)
    for (h in 1:H) {

      companion_matrix[1:K, (K*(h-1)+1):(K*h)] <- phi[[h+1]]

    }

    companion_matrix[cbind(c((K+1):(K*H)), c(1:(K*(H-1))))] <- 1

  }

  #### Eigenvalues and stability ####
  if (H == 0) {

    eigens = NULL
    stable = TRUE

  } else {

    eigens <- list(eigen(companion_matrix)$values)
    eigens[[2]] <- sqrt(Re(eigens[[1]])^2 + Im(eigens[[1]])^2)
    names(eigens) <- c('values', 'moduli')

    if (any(sqrt(Re(eigens$values)^2 + Im(eigens$values)^2) >= 1)) {

      stable <- FALSE
      warning("Process is unstable.")

    } else {

      stable <- TRUE

    }

  }

  #### Simulate observations ####
  a = MASS::mvrnorm(n = n + burn_ins, mu = rep(0, K), Sigma_a)
  z = matrix(0, nrow = n + burn_ins, ncol = K)

  if (H == 0) {

    for (k in 1:K) {

      z[, k] <- phi[[1]][k] + a[, k]

    }

  } else {

    for (t in (H+1):(n + burn_ins)) {

      z[t, ] <- t(phi[[1]])

      for (h in 1:H) {

        z[t, ] <- z[t, ] + phi[[h+1]] %*% z[t-h, ]

      }

      z[t, ] <- z[t, ] + a[t, ]

    }

  }

  a <- utils::tail(a, n)
  z <- utils::tail(z, n)
  rownames(a) <- 1:n
  rownames(z) <- 1:n

  #### Beta coefficient matrix ####
  beta = t(phi[[1]])

  if (H > 0) {

    for (h in 1:H) {

      beta = rbind(beta, phi[[h+1]])

    }

  }

  #### Expected value ####
  if (H == 0) {

    mu = as.matrix(phi[[1]])

  } else {

    mu = diag(x = 1, nrow = K, ncol = K)

    for (h in 1:H) {

      mu = mu - phi[[h+1]]

    }

    mu = solve(mu)
    mu = mu %*% phi[[1]]

  }

  #### Autocovariance ####

  if (H == 0) {

    autocovariance = NULL

  } else {

    Sigma_A = matrix(data = 0,
                     nrow = K*H,
                     ncol = K*H)

    Sigma_A[1:K, 1:K] <- Sigma_a
    Gamma_Z_0 <- solve(diag(x = 1, nrow = (K*H)^2, ncol = (K*H)^2) - kronecker(companion_matrix, companion_matrix)) %*% c(Sigma_A)
    Gamma_Z_0 <- matrix(Gamma_Z_0, nrow = K*H, ncol = K*H, byrow = TRUE)
    autocovariance <- list()

    for (h in 1:H) {

      autocovariance[[h]] <- Gamma_Z_0[1:K, (K*(h-1)+1):(K*h)]

    }

    # Recursively finding autocovariances
    for (g in (H+1):(max_autocovariance+1)) {

      autocovariance[[g]] <- matrix(data = 0, nrow = K, ncol = K)

      for (h in 1:H) {

        autocovariance[[g]] <- autocovariance[[g]] + phi[[h+1]] %*% autocovariance[[g-h]]

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
              eigens = eigens,
              H = H,
              K = K,
              mu = mu,
              n = n,
              phi = phi,
              Sigma_a = Sigma_a,
              stable = stable,
              z = z))

}
