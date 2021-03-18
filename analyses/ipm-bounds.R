tv_bound_estimator <- function(tau, lag, t) {
    d <- pmax(0, ceiling((tau - lag - t) / lag))
    return(d)
}

w1_bound_estimator <- function(x, y, t, tau, lag) {
    j_max <- ceiling((tau - lag - t) / lag)
    if (j_max > 0) {
        j_inds <- seq_len(j_max)
        x_inds <- ind0(t + j_inds * lag)
        y_inds <- ind0(t + (j_inds - 1) * lag)
        w1 <- sum(abs(x[x_inds] - y[y_inds]))
    } else {
        w1 <- 0
    }
    return(w1)
}
