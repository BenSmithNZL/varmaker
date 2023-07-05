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
                          matrix(c(1.25, 0, 0,
                                   0, 1, 0,
                                   0, 0, 0.5),
                                 nrow = 3,
                                 ncol = 3,
                                 byrow = TRUE),
                          '1000'),
               c("n must be an integer."))
})

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

autocovariance_0 <- matrix(c(0.131, 0.066,
                             0.066, 0.181),
                     nrow = 2,
                     ncol = 2,
                     byrow = TRUE)

autocovariance_1 <- matrix(c(0.072, 0.051,
                             0.104, 0.143),
                           nrow = 2,
                           ncol = 2,
                           byrow = TRUE)

autocovariance_2 <- matrix(c(0.046, 0.040,
                             0.113, 0.108),
                           nrow = 2,
                           ncol = 2,
                           byrow = TRUE)

autocovariance_3 <- matrix(c(0.035, 0.031,
                             0.093, 0.083),
                           nrow = 2,
                           ncol = 2,
                           byrow = TRUE)

test_that("Check autocovariances", {
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance[['0']], 3),
               autocovariance_0)
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance[['1']], 3),
               autocovariance_1)
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance[['2']], 3),
               autocovariance_2)
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$autocovariance[['3']], 3),
               autocovariance_3)
})

test_that("Check reproducability of results", {
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001)$z[100, ], 7),
               c(0.2554283, -0.3114040))
  expect_equal(round(create_var(cm_1, Sigma_a_1, 1001, seed = 100)$z[100, ], 7),
               c(-0.4943526, 0.8541217))
})


