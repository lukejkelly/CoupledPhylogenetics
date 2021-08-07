library("dplyr")

get_rwty_output <- function(out_dir, grid_a) {
    # the structure of rwty::load.trees makes it difficult to load trees
    # without a parameter file so we do it manually instead
    # rwty operates on unrooted trees
    grid_e <- tibble(grid_a, rwty = list(NA))
    for (i in seq_len(nrow(grid_e))) {
        svMisc::progress(i, nrow(grid_e))

        grid_e_i <- grid_e[i, ]
        trees <- get_rwty_trees(out_dir, grid_e_i)
        pars <- get_rwty_pars(out_dir, grid_e_i)

        output <- list(trees = trees, ptable = pars, gens.per.tree = 1)
        class(output) <- "rwty.chain"
        grid_e$rwty[[i]] <- output
    }
    message("rwty trees read")
    return(grid_e)
}

get_rwty_trees <- function(out_dir, grid_i) {
    trees <- get_trees_(out_dir, grid_i, grid_i$lag, grid_i$c, "_x") %>%
        `[`(ind0(0, grid_i$run_length / grid_i$sample_interval)) %>%
        map(unroot)
    class(trees) <- "multiPhylo"
    return(trees)
}

get_rwty_pars <- function(out_dir, grid_i) {
    pars <- get_pars_(out_dir, grid_i, grid_i$lag, grid_i$c, "_x") %>%
        select(-c(rho, log_likelihood)) %>%
        janitor::remove_constant() %>%
        slice(ind0(0, grid_i$run_length / grid_i$sample_interval)) %>%
        as.data.frame()
    return(pars)
}

win_size <- function(y, n_win = 5) {
    n_steps <- y$run_length / y$sample_interval
    ws <- floor(n_steps / n_win)
    return(ws)
}
x_breaks <- function(y, n_win = 5) {
    xb <- seq.int(1, y$run_length / y$sample_interval, win_size(y, n_win)) %>%
        as.character()
    return(xb)
}
