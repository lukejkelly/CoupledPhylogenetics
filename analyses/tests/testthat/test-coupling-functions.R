library("testthat")

test_that("get tau", {
    # Synthetic data and arbitrary coupling time
    n <- 10 + rgeom(1, 1e-2)
    x <- rnorm(n + 2)
    y <- rnorm(n + 1)
    t <- floor(n * rbeta(1, 6, 3))
    x[ind_x(t, n + 1)] <- y[ind_y(t - 1, n)]
    expect_equal(get_tau(x, y), t)
    expect_true(is.na(get_tau(rnorm(n + 2), rnorm(n + 1))))
})
