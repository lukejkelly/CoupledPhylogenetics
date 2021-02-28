tv_bound_estimator <- function(tau, lag, t) {
    d <- pmax(0, ceiling((tau - lag - t) / lag))
    return(d)
}

w1_bound_estimators <- function(tau, lag, iters, x, y) {
    ds <- tv_bound_estimator(tau, lag, iters)
    ws <- double(length(iters))
    for (i in seq_along(ws)) {
        x_inds <- iters[i] + 1 + seq_len(ds[i])
        y_inds <- x_inds - 1
        ws[i] <- sum(abs(x[x_inds] - y[y_inds]))
    }
    return(ws)
}
