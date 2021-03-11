# functions to index sample output
# TODO: add initial lag/subsample once TraitLab is updated
# currently, lag l = subsample interval s so
#   x[1] = x_0, x[2] = x_s = x_l, x[3] = y_2s, ..., x[n+2] = x_(n+1)s
#               y[2] = y_0,       y[2] = y_s, ...,  y[n+1] = y_ns,
# where we have vertically aligned the coupled states
# this will shortly change to allow l = ts so
#   x[1] = x_0, ..., x[t+1] = x_ts = x_l, ..., x_[t+n+1] = x_(n+t)s
#                      y[1] = y_0,        ...,   y_[n+1] = y_ns
# all else equal, the ground truth chain z has the same indexing as y

ind_x <- function(k, m) {
    x_inds <- seq.int(k, m) + 1
    return(x_inds)
}
ind_y <- function(k, m) {
    y_inds <- seq.int(k, m) + 1
    return(y_inds)
}
ind_z <- ind_y

monte_carlo_estimator <- function(x, k, m) {
    # (x_k + ... + x_m) / (m - k + 1)
    mc <- mean(x[ind_x(k, m)])
    return(mc)
}

unbiased_estimator <- function(x, y, k, m, t) {
    # monte carlo term mc(x_k:m) + bias correction bc(x_(k+1):(t-1), y_k:(t-2))
    mc <- monte_carlo_estimator(x, k, m)
    if (k + 1 > t - 1) {
        bc <- 0
    } else {
        w_inds <- seq.int(k + 1, t - 1)
        x_inds <- ind_x(k + 1, t - 1)
        y_inds <- ind_y(k, t - 2)
        w <- pmin(1, (w_inds - k) / (m - k + 1))
        d <- x[x_inds] - y[y_inds]
        bc <- sum(w * d)
    }
    ue <- mc + bc
    return(list(mc = mc, bc = bc, ue = ue))
}
