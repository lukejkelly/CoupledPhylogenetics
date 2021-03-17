# functions to index sample output
# previously, lag l = subsample interval s so
#   x[1] = x_0, x[2] = x_s = x_l, x[3] = y_2s, ..., x[n+1] = x_ns
#               y[1] = y_0,       y[2] = y_s, ...,    y[n] = y_(n-1)s,
# where we have vertically aligned the coupled states
# this has now changed to allow l = ks so
#   x[1] = x_0, ..., x[k+1] = x_ks = x_l, ...,   x_[n+1] = x_ns
#                      y[1] = y_0,        ..., y_[n+1-k] = y_(n-k)s
# all else equal, the ground truth chain z has the same indexing as y

ind0 <- function(k, m = NULL) {
    # outputs are indexed from 0, k can be a sequence if m is NULL
    if (is.null(m)) {
        inds <- k + 1
    } else {
        inds <- seq.int(k, m) + 1
    }
    return(inds)
}

monte_carlo_estimator <- function(x, k, m) {
    # (x_k + ... + x_m) / (m - k + 1)
    mc <- mean(x[ind0(k, m)])
    return(mc)
}

unbiased_estimator <- function(x, y, k, m, tau, lag) {
    # monte carlo estimator + bias correction
    mc <- monte_carlo_estimator(x, k, m)
    bc <- bias_correction(x, y, k, m, tau, lag)
    ue <- mc + bc
    return(list(mc = mc, bc = bc, ue = ue))
}

bias_correction <- function(x, y, k, m, tau, lag) {
    bc <- 0
    for (t in seq.int(k, m)) {
        j_max <- ceiling((tau - lag - t) / lag)
        if (j_max > 0) {
            j_inds <- seq_len(j_max)
            x_inds <- ind0(t + j_inds * lag)
            y_inds <- ind0(t + (j_inds - 1) * lag)
            bc_t <- sum(x[x_inds] - y[y_inds])
        } else {
            bc_t <- 0
        }
        bc <- bc + bc_t
    }
    bc <- bc / (m - k + 1)
    return(bc)
}
