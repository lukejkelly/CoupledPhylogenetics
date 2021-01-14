unbiased_estimator <- function(x, y, k, m, t) {
    # Assume x and y both include initial state so x[1] = x_0, y[1] = y_0 and
    # proceed at lag/subsample l, so x[i] = x_{(i - 1) * l}, y is likewise, and
    # x[i] coupled to y[i - 1]
    # if (t > m) {
    #     mc <- NA
    #     bc <- NA
    #     ue <- NA
    # } else {
        x_inds <- seq.int(k, m)
        b_inds <- seq.int(k + 1, t - 1,
                          length.out = max((t - 1) - (k + 1) + 1, 0))
        mc <- mean(x[x_inds + 1])

        w <- pmin(1, (b_inds - k) / (m - k + 1))
        d <- x[(b_inds + 1)] - y[(b_inds + 1) - 1]
        bc <- w %*% d

        ue <- mc + bc
    # }
    return(c(mc = mc, bc = bc, ue = ue))
}
