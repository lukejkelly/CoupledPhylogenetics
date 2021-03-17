library("testthat")

test_that("get tau when offset is 1", {
    n <- 10 + rgeom(1, 1e-2)
    lag_offset <- 1
    x <- rnorm(n + 1)
    y <- rnorm(n + 1 - lag_offset)
    t <- floor(length(x) * rbeta(1, 6, 3))
    x[ind0(t, n)] <- y[ind0(t - lag_offset, n - lag_offset)]
    expect_equal(get_tau(x, y, lag_offset), t)
    expect_error(get_tau(rnorm(length(x)), rnorm(length(y)), lag_offset))
})

test_that("get tau when offset is 3", {
    n <- 10 + rgeom(1, 1e-2)
    lag_offset <- 3
    x <- rnorm(n + 1)
    y <- rnorm(n + 1 - lag_offset)
    t <- floor(length(x) * rbeta(1, 6, 3))
    x[ind0(t, n)] <- y[ind0(t - lag_offset, n - lag_offset)]
    expect_equal(get_tau(x, y, lag_offset), t)
    expect_error(get_tau(rnorm(length(x)), rnorm(length(y)), lag_offset))
})
