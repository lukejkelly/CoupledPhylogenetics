if (Sys.info()["nodename"] == "kelly.local") {
    source("../../RSPR/R/rspR.r")
} else {
    source("/home/users/kelly/rspr/R/rspR.r")
}

compute_tree_distances <- function(out_dir, grid_a) {
    # Rooted SPR distance between pairs of trees drawn from coupled kernel
    for (i in seq_len(nrow(grid_a))) {
        svMisc::progress(i, nrow(grid_a))
        grid_a_i <- grid_a[i, ]

        x <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_x")
        y <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_y")

        lag_offset <- grid_a_i$lag / grid_a_i$sample_interval

        # d <- purrr::map2_dbl(x[-seq_len(lag_offset)], y, ape::dist.topo,
        #                      "score")
        d <- purrr::map2_dbl(x[-seq_len(lag_offset)], y, rspr)

        output_file_name <- make_file_name_(out_dir, grid_a_i, grid_a_i$lag,
                                            grid_a_i$c, "", "dist")
        readr::write_lines(d, output_file_name)
    }
    message("tree distances computed")
}

compute_tree_jump_distances <- function(out_dir, grid_a) {
    # Rooted SPR distance between pairs of trees in a chain
    grid_s <- tibble(grid_a, mean = NA_real_, q90 = NA_real_)
    for (i in seq_len(nrow(grid_a))) {
        svMisc::progress(i, nrow(grid_a))
        grid_s_i <- grid_s[i, ]

        lag_offset <- grid_s_i$lag / grid_s_i$sample_interval
        x <- get_trees_(out_dir, grid_s_i, grid_s_i$lag, grid_s_i$c, "_x")
        x0 <- x[-seq_len(lag_offset)]
        p <- sample.int(min(200, 2 * floor(length(x0) / 2)))
        x1 <- x0[p[seq_len(length(p) / 2)]]
        x2 <- x0[p[-seq_len(length(p) / 2)]]
        d <- purrr::map2_dbl(x1, x2, rspr)

        grid_s$mean[i] <- mean(d)
        grid_s$q90[i] <- quantile(d, 0.9)
    }
    message("tree distances computed")
    return(grid_s)
}
