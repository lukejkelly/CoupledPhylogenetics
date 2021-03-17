library("testthat")

test_that("indexing from 0", {
    k <- rgeom(1, 1e-1)
    m <- k + rgeom(1, 1e-1)
    inds_obs1 <- ind0(k, m)
    inds_obs2 <- ind0(k:m)
    inds_exp <- (k + 1):(m + 1)
    expect_equal(inds_obs1, inds_exp)
    expect_equal(inds_obs2, inds_exp)
})


test_that("monte carlo estimator", {
    n <- rgeom(1, 1e-2)
    x <- rnorm(n + 1)
    expect_equal(monte_carlo_estimator(x, 0, n), mean(x))
    m <- floor(n * runif(1))
    k <- floor(m * runif(1))
    expect_equal(monte_carlo_estimator(x, k, m), mean(x[(k + 1):(m + 1)]))
})

test_that("unbiased estimator", {
    warning("implement additional checks for lag_offset > 1")
    for (i in c(rep(1, 100), rep(seq_len(10), 10))) {
        lag_offset <- i
        n <- 20 + lag_offset + rgeom(1, 1e-2)
        x <- rnorm(n + 1)
        y <- rnorm(n + 1 - i)
        # boundary checks
        k0 <- floor(n / 2) - 2
        m0 <- n
        t0a <- k0 + lag_offset
        d0a <- ind0(t0a) - length(x)
        if (d0a > 0) {
            x <- c(x, rnorm(d0a))
            y <- c(y, rnorm(d0a))
        }
        obs0a <- unbiased_estimator(x, y, k0, m0, t0a, lag_offset)
        expect_equal(obs0a$mc, obs0a$ue)
        expect_equal(obs0a$bc, 0)
        if (lag_offset == 1) {
            exp0a <- unbiased_estimator_(x, y, k0, m0, k0 + 1)
            expect_equal(obs0a$mc, exp0a$mc)
            expect_equal(obs0a$bc, exp0a$bc)
            expect_equal(obs0a$ue, exp0a$ue)
        }
        t0b <- k0 + lag_offset + 1
        d0b <- ind0(t0a) > length(x)
        if (d0b - 0) {
            x <- c(x, rnorm(d0b))
            y <- c(y, rnorm(d0b))
        }
        obs0b <- unbiased_estimator(x, y, k0, m0, t0b, lag_offset)
        expect_false(obs0b$mc == obs0b$ue)
        expect_false(obs0b$bc == 0)
        if (lag_offset == 1) {
            exp0b <- unbiased_estimator_(x, y, k0, m0, k0 + 2)
            expect_equal(obs0b$mc, exp0b$mc)
            expect_equal(obs0b$bc, exp0b$bc)
            expect_equal(obs0b$ue, exp0b$ue)
        }
        # t < k < m
        m1 <- floor(n * rbeta(1, 6, 3))
        k1 <- floor(m1 * rbeta(1, 6, 3))
        t1 <- floor(k1 * rbeta(1, 6, 3))
        b1 <- unbiased_estimator(x, y, k1, m1, t1, lag_offset)$bc
        expect_equal(b1, 0)
        if (lag_offset == 1) {
            exp1 <- unbiased_estimator_(x, y, k1, m1, t1)$bc
            expect_equal(b1, exp1)
        }
        # k < t < m
        m2 <- floor(n * rbeta(1, 6, 3))
        t2 <- floor(m2 * rbeta(1, 6, 3))
        k2 <- floor(t2 * rbeta(1, 6, 3))
        b2 <- unbiased_estimator(x, y, k2, m2, t2, lag_offset)$bc
        if (t2 - lag_offset <= k2) {
            expect_equal(b2, 0)
        } else {
            expect_false(b2 == 0)
        }
        if (lag_offset == 1) {
            if (t2 - 1 < k2 + 1) {
                exp2 <- 0
            } else {
                w2 <- seq.int(1, t2 - 1 - k2) / (m2 - k2 + 1)
                d2 <- x[ind0(k2 + 1, t2 - 1)] - y[ind0(k2, t2 - 2)]
                exp2 <- sum(w2 * d2)
            }
            expect_equal(b2, exp2)
            expect_equal(b2, unbiased_estimator_(x, y, k2, m2, t2)$bc)
        }
        # k < m < t
        t3 <- floor(n * rbeta(1, 6, 3))
        m3 <- floor(t3 * rbeta(1, 6, 3))
        k3 <- floor(m3 * rbeta(1, 6, 3))
        b3 <- unbiased_estimator(x, y, k3, m3, t3, lag_offset)$bc
        if (t3 > lag_offset + k3) {
            expect_false(b3 == 0)
        } else {
            expect_true(b3 == 0)
        }
        if (lag_offset == 1) {
            i3 <- seq.int(1, t3 - 1 - k3) / (m3 - k3 + 1)
            if (t3 - 1 < m3 + 1) {
                w3 <- i3
            } else {
                w3 <- pmin(1, i3)
            }
            d3 <- x[ind0(k3 + 1, t3 - 1)] - y[ind0(k3, t3 - 2)]
            exp3 <- sum(w3 * d3)
            expect_equal(b3, exp3)
            expect_equal(b3, unbiased_estimator_(x, y, k3, m3, t3)$bc)
        }
    }
})
