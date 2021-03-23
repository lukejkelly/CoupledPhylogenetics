get_rwty_trees <- function(out_dir, grid_a) {
    grid_e <- tibble(grid_a, rwty = list(NA))
    for (i in seq_len(nrow(grid_e))) {
        svMisc::progress(i, nrow(grid_e))

        g_i <- grid_e[i, ]
        trees <- get_trees_(out_dir, g_i, g_i$lag, g_i$c, "_x") %>%
            `[`(seq_len(ind0(g_i$run_length / g_i$sample_interval))) %>%
            map(unroot)
        class(trees) <- "multiPhylo"

        output <- list(trees = trees, ptable = NULL, gens.per.tree = 1)
        class(output) <- "rwty.chain"

        grid_e$rwty[[i]] <- output
    }
    message("rwty trees read")
    return(grid_e)
}
