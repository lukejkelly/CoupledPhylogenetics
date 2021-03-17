unbiased_estimator_ <- function(x, y, k, m, t) {
    # monte carlo estimator + bias correction
    mc <- monte_carlo_estimator(x, k, m)
    if (k + 1 > t - 1) {
        bc <- 0
    } else {
        w_inds <- seq.int(k + 1, t - 1)
        x_inds <- ind0(k + 1, t - 1)
        y_inds <- ind0(k, t - 2)
        w <- pmin(1, (w_inds - k) / (m - k + 1))
        d <- x[x_inds] - y[y_inds]
        bc <- sum(w * d)
    }
    ue <- mc + bc
    return(list(mc = mc, bc = bc, ue = ue))
}
