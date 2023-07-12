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

#### Check parameters ####
test_that("Check that first element is a column vector", {
  expect_error(
    create_vma(
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
    c("theta_0 is not a vector."))
})

test_that("Check that all the other coefficient matricies are K by K.", {
  expect_error(create_vma(list(c(1, 0.5, 0.25),
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
               c("theta_1 is not a 3 by 3 matrix."))
  expect_error(create_vma(list(c(1, 0.5, 0.25),
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
               c("theta_2 is not a 3 by 3 matrix."))
})

test_that("Check Sigma_a is symmetric", {
  expect_error(create_vma(list(c(1, 0.5, 0.25),
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
  expect_error(create_vma(cm_1, Sigma_a_1, '1001'),
               c("n must be an integer."))
  expect_error(create_vma(cm_1, Sigma_a_1, c(100, 285)),
               c("n must be an integer."))
})
