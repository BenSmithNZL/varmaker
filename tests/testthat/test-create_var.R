test_that("Check that first element is a column vector", {
    expect_error(
      create_var(
        list(matrix(c(1, 0.5, 0.25, 2),
                    nrow = 2,
                    ncol = 2,
                    byrow = T),
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
      c("phi_0 is not a column vector."))
})

test_that("Check that all the other coefficient matricies are K by K.", {
  expect_error(create_var(list(matrix(c(1, 0.5, 0.25),
                                      nrow = 3,
                                      ncol = 1,
                                      byrow = T),
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
})

