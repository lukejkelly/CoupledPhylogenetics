samples_k_to_m <- function(x, k, m) {
    # x indexed from 0 so x[1] = x_0, we want x_k, ..., x_m
    inds <- seq.int(k, m) + 1
    x_km <- x[inds]
    return(x_km)
}

monte_carlo_estimator <- function(x, k, m) {
    # Assume x includes initial state so x[1] = x_0
    x_km <- samples_k_to_m(x, k, m)
    mc <- mean(x_km)
    return(mc)
}

unbiased_estimator <- function(x, y, k, m, t) {
    # Assume x and y both include initial state so x[1] = x_0, y[1] = y_0 and
    # proceed at lag/subsample l, so x[i] = x_{(i - 1) * l}, y is likewise
    mc <- monte_carlo_estimator(x, k, m)
    n_inds <- max((t - 1) - (k + 1) + 1, 0)
    b_inds <- seq.int(k + 1, t - 1, length.out = n_inds)
    w <- pmin(1, (b_inds - k) / (m - k + 1))
    d <- x[(b_inds + 1)] - y[(b_inds + 1) - 1]
    bc <- w %*% d
    ue <- mc + bc
    return(c(mc = mc, bc = bc, ue = ue))
}
