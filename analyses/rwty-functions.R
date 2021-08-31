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
        `[`(ind0(1, grid_i$run_length / grid_i$sample_interval)) %>%
        map(unroot)
    class(trees) <- "multiPhylo"
    return(trees)
}

get_rwty_pars <- function(out_dir, grid_i) {
    pars <- get_pars_(out_dir, grid_i, grid_i$lag, grid_i$c, "_x") %>%
        select(-c(rho, log_likelihood)) %>%
        janitor::remove_constant() %>%
        slice(ind0(1, grid_i$run_length / grid_i$sample_interval)) %>%
        as.data.frame()
    return(pars)
}

win_size <- function(y, n_win = 10) {
    n_steps <- y$run_length / y$sample_interval
    ws <- floor(n_steps / n_win)
    return(ws)
}
x_end <- function(y, n_win = 10) {
    xb <- seq_len(n_win) * win_size(y, n_win)
    return(xb)
}

my_asdsf_plot <- function(grid_e, fig_rwty, chains, burnin, n_win,
                          min.freq) {
    list_L <- grid_e %>% select(L) %>% unique() %>% pull(L) %>% as.list()
    fig_rwty_stats <- map(
        fig_rwty,
        ~ .$data %>% select(Generation, ASDSF, upper.95, lower.95)
    )
    my_fig_rwty_data <- map2(
        list_L,
        fig_rwty_stats,
        ~ tibble(L = .x, .y, xend = x_end(grid_e[1, ], n_win))
    ) %>%
        bind_rows()
    my_fig_rwty <- my_fig_rwty_data %>%
        ggplot(aes(x = Generation, xend = xend)) +
        geom_segment(aes(y = ASDSF, yend = ASDSF, colour = "ASDSF")) +
        # geom_segment(aes(y = lower.95, yend = lower.95, colour = "95%")) +
        # geom_segment(aes(y = upper.95, yend = upper.95, colour = "95%")) +
        # geom_pointrange(aes(ymin = lower.95, ymax = upper.95,
        #                     colour = "asdsf"),
        #                 show.legend = TRUE) +
        geom_hline(aes(yintercept = 0.01), alpha = 0.25, linetype = "dashed") +
        ylim(0, NA) +
        labs(title = "ASDSF",
             subtitle = sprintf("minimum %.02e iterations, replications = %d",
                                grid_e$run_length[1], n_distinct(grid_e$c)),
             x = sprintf("iteration / %d", grid_e$sample_interval[1]),
             y = "ASDSF",
             linetype = NULL, colour = NULL
         )
    if (length(list_L) > 1) {
        for (scales in c("free", "fixed")) {
            fig_p <- my_fig_rwty +
                facet_wrap(~ L, ncol = length(list_L), scales = scales,
                           labeller = "label_both")
            ggsave(sprintf(fig_template, sprintf("asdsf_axes-%s", scales)),
                   fig_p,
                   width = 3 * length(list_L) + 2,
                   height = 3)
        }
    } else {
        ggsave(sprintf(fig_template, "asdsf"), my_fig_rwty, width = 5,
               height = 3)
    }
}
