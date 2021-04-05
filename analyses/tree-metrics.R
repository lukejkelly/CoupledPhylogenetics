compute_tree_distances <- function(out_dir, grid_a) {

    for (i in seq_len(nrow(grid_a))) {
        svMisc::progress(i, nrow(grid_a))
        grid_a_i <- grid_a[i, ]

        x <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_x")
        y <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_y")

        lag_offset <- grid_a_i$lag / grid_a_i$sample_interval

        d <- purrr::map2_dbl(x[-seq_len(lag_offset)], y, ape::dist.topo,
                             "score")

        output_file_name <- make_file_name_(out_dir, grid_a_i, grid_a_i$lag,
                                            grid_a_i$c, "", "dist")
        readr::write_lines(d, output_file_name)
    }
    message("tree distances computed")
}
