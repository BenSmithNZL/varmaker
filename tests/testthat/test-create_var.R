#### Constants ####
cm_1 <- list(c(0, 0),
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


cm_2 <- list(c(0, 0, 0), #2.1.14
             matrix(c(0.5, 0, 0,
                      0.1, 0.1, 0.3,
                      0, 0.2, 0.3),
                    nrow = 3,
                    ncol = 3,
                    byrow = TRUE))

Sigma_a_2 <- matrix(c(2.25, 0, 0,
                      0, 1, 0.5,
                      0, 0.5, 0.74),
                    nrow = 3,
                    ncol = 3,
                    byrow = TRUE)

autocovariance_1_0 <- matrix(c(0.131, 0.066,
                               0.066, 0.181),
                             nrow = 2,
                             ncol = 2,
                             byrow = TRUE)

autocovariance_1_1 <- matrix(c(0.072, 0.051,
                               0.104, 0.143),
                             nrow = 2,
                             ncol = 2,
                             byrow = TRUE)

autocovariance_1_2 <- matrix(c(0.046, 0.040,
                               0.113, 0.108),
                             nrow = 2,
                             ncol = 2,
                             byrow = TRUE)

autocovariance_1_3 <- matrix(c(0.035, 0.031,
                               0.093, 0.083),
                             nrow = 2,
                             ncol = 2,
                             byrow = TRUE)

autocovariance_2_0 <- matrix(c(3.000, 0.161, 0.019,
                               0.161, 1.172, 0.674,
                               0.019, 0.674, 0.954),
                             nrow = 3,
                             ncol = 3,
                             byrow = TRUE)

autocovariance_2_1 <- matrix(c(1.500, 0.080, 0.009,
                               0.322, 0.335, 0.355,
                               0.038, 0.437, 0.421),
                             nrow = 3,
                             ncol = 3,
                             byrow = TRUE)

autocovariance_2_2 <- matrix(c(0.750, 0.040, 0.005,
                               0.194, 0.173, 0.163,
                               0.076, 0.198, 0.197),
                             nrow = 3,
                             ncol = 3,
                             byrow = TRUE)

autocorrelation_1_0 <- matrix(c(1, 0.086, 0.011,
                                0.086, 1, 0.637,
                                0.011, 0.637, 1),
                              nrow = 3,
                              ncol = 3,
                              byrow = TRUE)

autocorrelation_1_1 <- matrix(c(0.500, 0.043, 0.006, #Autocorrelations are different from Lutkepohl's book due to rounding errors
                                0.172, 0.286, 0.336,
                                0.022, 0.413, 0.441),
                              nrow = 3,
                              ncol = 3,
                              byrow = TRUE)

autocorrelation_1_2 <- matrix(c(0.250, 0.021, 0.003,
                                0.103, 0.147, 0.154,
                                0.045, 0.187, 0.207),
                              nrow = 3,
                              ncol = 3,
                              byrow = TRUE)


#### Check parameters ####
test_that("Check that first element is a column vector", {
    expect_error(
        create_var(
          list(2,
             matrix(c(0.2, 0, 0.1,
                      -0.3, 0, 0,
                      0, 0, -0.1),
                    nrow = 3,
                    ncol = 3,
                    byrow = T),
             matrix(c(0, 0.3, 0,
                      0, 0, 0,
                      0, 0, 0.1),
                    nrow = 3,
                    ncol = 3,
                    byrow = T)),
        diag(x = 1, nrow = 3, ncol = 3),
        1000),
      c("phi_0 is not a vector."))
})

test_that("Check that all the other coefficient matricies are K by K.", {
  expect_error(create_var(list(c(1, 0.5, 0.25),
                               matrix(c(0.2, 0,
                                        -0.3, 0),
                                      nrow = 2,
                                      ncol = 2,
                                      byrow = T),
                               matrix(c(0, 0.3, 0,
                                        0, 0, 0,
                                        0, 0, 0.1),
                                      nrow = 3,
                                      ncol = 3,
                                      byrow = T)),
                          diag(x = 1, nrow = 3, ncol = 3),
                          1000),
               c("phi_1 is not a 3 by 3 matrix."))
  expect_error(create_var(list(c(1, 0.5, 0.25),
                               matrix(c(0.2, 0, 0.1,
                                        -0.3, 0, 0,
                                        0, 0, -0.1),
                                      nrow = 3,
                                      ncol = 3,
                                      byrow = T),
                               matrix(c(0, 0.3,
                                        0, 0,
                                        0, 0),
                                      nrow = 3,
                                      ncol = 2,
                                      byrow = T)),
                          diag(x = 1, nrow = 3, ncol = 3),
                          1000),
               c("phi_2 is not a 3 by 3 matrix."))
})

test_that("Check Sigma_a is symmetric", {
  expect_error(create_var(list(c(1, 0.5, 0.25),
                               matrix(c(0.2, 0, 0.1,
                                        -0.3, 0, 0,
                                        0, 0, -0.1),
                                      nrow = 3,
                                      ncol = 3,
                                      byrow = TRUE),
                               matrix(c(0, 0.3, 0,
                                        0, 0, 0,
                                        0, 0, 0.1),
                                      nrow = 3,
                                      ncol = 3,
                                      byrow = TRUE)),
                          matrix(c(1.25, 1, 0,
                                   0, 1, 0,
                                   0, 0, 0.5),
                                 nrow = 3,
                                 ncol = 3,
                                 byrow = TRUE),
                          1000),
               c("Sigma_a is not a symmetric matrix"))
})

test_that("Check n is an integer", {
  expect_error(create_var(cm_1, Sigma_a_1, '1001'),
               c("n must be an integer."))
  expect_error(create_var(cm_1, Sigma_a_1, c(100, 285)),
               c("n must be an integer."))
})

#### Check answers ####
test_that("Check autocovariances", {
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance$Gamma_0, 3),
               autocovariance_1_0)
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance$Gamma_1, 3),
               autocovariance_1_1)
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance$Gamma_2, 3),
               autocovariance_1_2)
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance$Gamma_3, 3),
               autocovariance_1_3)
  expect_equal(round(create_var(cm_2, Sigma_a_2, 1001)$autocovariance$Gamma_0, 3),
               autocovariance_2_0)
  expect_equal(round(create_var(cm_2, Sigma_a_2, 1001)$autocovariance$Gamma_1, 3),
               autocovariance_2_1)
  expect_equal(round(create_var(cm_2, Sigma_a_2, 1001)$autocovariance$Gamma_2, 3),
               autocovariance_2_2)
})

test_that("Check autocorrelations", {
  expect_equal(round(create_var(cm_2, Sigma_a_2, 1001)$autocorrelation$R_0, 3),
               autocorrelation_1_0)
  expect_equal(round(create_var(cm_2, Sigma_a_2, 1001)$autocorrelation$R_1, 3),
               autocorrelation_1_1)
  expect_equal(round(create_var(cm_2, Sigma_a_2, 1001)$autocorrelation$R_2, 3),
               autocorrelation_1_2)
})

test_that("Check reproducability of results", {
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$z[100, ], 7),
               c(0.2554283, -0.3114040))
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001, seed = 100)$z[100, ], 7),
               c(-0.4943526, 0.8541217))
})


