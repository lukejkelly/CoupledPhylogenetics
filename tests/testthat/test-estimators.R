library("testthat")

test_that("index x", {
    k <- rgeom(1, 1e-1)
    m <- k + rgeom(1, 1e-1)
    inds_obs <- ind_x(k, m)
    inds_exp <- (k + 1):(m + 1)
    expect_equal(inds_obs, inds_exp)
})

test_that("index y", {
    k <- rgeom(1, 1e-1)
    m <- k + rgeom(1, 1e-1)
    inds_obs <- ind_y(k, m)
    inds_exp <- (k + 1):(m + 1)
    expect_equal(inds_obs, inds_exp)
})

test_that("index z", {
    k <- rgeom(1, 1e-1)
    m <- k + rgeom(1, 1e-1)
    inds_obs <- ind_z(k, m)
    inds_exp <- (k + 1):(m + 1)
    expect_equal(inds_obs, inds_exp)
})

test_that("monte carlo estimator", {
    n <- rgeom(1, 1e-2)
    x <- rnorm(n + 2)
    expect_equal(monte_carlo_estimator(x, 0, n + 1), mean(x))
    m <- floor(n * runif(1))
    k <- floor(m * runif(1))
    expect_equal(monte_carlo_estimator(x, k, m), mean(x[(k + 1):(m + 1)]))
})

test_that("unbiased estimator", {
    for (i in seq_len(20)) {
        n <- 10 + rgeom(1, 1e-2)
        x <- rnorm(n + 2)
        y <- rnorm(n + 1)
        # boundary checks
        k0 <- floor(n / 2) - 2
        obs0a <- unbiased_estimator(x, y, k0, n, k0 + 1)
        expect_equal(obs0a$mc, obs0a$ue)
        expect_equal(obs0a$bc, 0)
        obs0b <- unbiased_estimator(x, y, k0, n, k0 + 2)
        expect_false(obs0b$mc == obs0b$ue)
        expect_false(obs0b$bc == 0)
        # t < k < m
        m1 <- floor(n * rbeta(1, 6, 3))
        k1 <- floor(m1 * rbeta(1, 6, 3))
        t1 <- floor(k1 * rbeta(1, 6, 3))
        expect_equal(unbiased_estimator(x, y, k1, m1, t1)$bc, 0)
        # k < t < m
        m2 <- floor(n * rbeta(1, 6, 3))
        t2 <- floor(m2 * rbeta(1, 6, 3))
        k2 <- floor(t2 * rbeta(1, 6, 3))
        if (t2 - 1 < k2 + 1) {
            exp2 <- 0
        } else {
            w2 <- seq.int(1, t2 - 1 - k2) / (m2 - k2 + 1)
            d2 <- x[ind_x(k2 + 1, t2 - 1)] - y[ind_y(k2, t2 - 2)]
            exp2 <- w2 %*% d2
        }
        expect_equal(unbiased_estimator(x, y, k2, m2, t2)$bc, exp2)
        # k < m < t
        t3 <- floor(n * rbeta(1, 6, 3))
        m3 <- floor(t3 * rbeta(1, 6, 3))
        k3 <- floor(m3 * rbeta(1, 6, 3))
        i3 <- seq.int(1, t3 - 1 - k3) / (m3 - k3 + 1)
        if (t3 - 1 < m3 + 1) {
            w3 <- i3
        } else {
            w3 <- pmin(1, i3)
        }
        d3 <- x[ind_x(k3 + 1, t3 - 1)] - y[ind_y(k3, t3 - 2)]
        exp3 <- w3 %*% d3
        expect_equal(unbiased_estimator(x, y, k3, m3, t3)$bc, exp3)
    }
})
