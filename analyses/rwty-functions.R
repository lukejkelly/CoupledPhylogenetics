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

my_asdsf_plot <- function(grid_e, fig_rwty) {

    list_L <- grid_e %>% select(L) %>% unique() %>% pull(L)
    fig_rwty_stats <- map(
        fig_rwty,
        ~ .$data %>% select(Generation, ASDSF)
    )
    my_fig_rwty_data <- tibble(L = list_L, fig_rwty_stats) %>%
        unnest(fig_rwty_stats)

    my_fig_rwty <- my_fig_rwty_data %>%
        ggplot(aes(x = Generation, y = ASDSF)) +
        geom_line(aes(linetype = "ASDSF"), alpha = 0.5) +
        geom_point(aes(shape = "ASDSF")) +
        geom_hline(aes(yintercept = 0.01), alpha = 0.25, linetype = "dashed") +
        labs(
            title = "ASDSF",
            subtitle = sprintf(
                "minimum %.02e iterations, replications = %d",
                grid_e$run_length[1],
                n_distinct(grid_e$c)
            ),
            x = sprintf("iteration / %d", grid_e$sample_interval[1]),
            y = "ASDSF",
            linetype = NULL,
            shape = NULL
         ) +
         xlim(0, NA) +
         scale_y_log10() +
         theme_light()
    if (length(list_L) > 1) {
        for (scales in c("free", "fixed")) {
            fig_p <- my_fig_rwty +
                facet_wrap(
                    ~ L,
                    ncol = length(list_L),
                    scales = scales,
                    labeller = "label_both"
                )
            ggsave(
                sprintf(fig_template, sprintf("asdsf_axes-%s", scales)),
                fig_p,
                width = 3 * length(list_L) + 2,
                height = 3
            )
        }
    } else {
        ggsave(
            sprintf(fig_template, "asdsf"),
            my_fig_rwty,
            width = 5,
            height = 3
        )
    }
}
