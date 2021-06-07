# average distances between components up to coupling
get_par_data <- function(out_dir, grid_a) {
    grid_p <- tibble(grid_a, d_log_prior = NA, d_log_llkd = NA,
                     d_root_time = NA, d_mu = NA, d_spr = NA)
    for (i in seq_len(nrow(grid_p))) {
        svMisc::progress(i, nrow(grid_p))
        grid_p_i <- grid_p[i, ]

        tau_i <- grid_p_i$tau
        lag_offset <- grid_p_i$lag / grid_p_i$sample_interval

        # sampled x and y
        x <- get_pars_(out_dir, grid_p_i, grid_p_i$lag, grid_p_i$c, "_x") %>%
            select(log_prior, integrated_llkd, root_time, mu) %>%
            slice(ind0(lag_offset, tau_i))
        y <- get_pars_(out_dir, grid_p_i, grid_p_i$lag, grid_p_i$c, "_y")  %>%
            select(log_prior, integrated_llkd, root_time, mu) %>%
            slice(ind0(0, tau_i - lag_offset))
        # pre-computed tree distances
        z <- read.table(make_file_name_(out_dir, grid_p_i, grid_p_i$lag,
                                        grid_p_i$c, "", "dist")) %>%
            slice(ind0(0, tau_i - lag_offset))

        grid_p$d_log_prior[i] <- mean(x$log_prior == y$log_prior)
        grid_p$d_log_llkd[i] <- mean(x$integrated_llkd == y$integrated_llkd)
        grid_p$d_root_time[i] <- mean(x$root_time == y$root_time)
        grid_p$d_mu[i] <- mean(x$mu == y$mu)
        grid_p$d_spr[i] <- mean(z$V1 == 0)
    }
    message("parameter distances computed")
    return(grid_p)
}

# reshape data and plot rates
make_par_plot <- function(fig_template, grid_p) {
    warning("this function was for a specific example")
    fig_par_data <- grid_p %>%
        pivot_longer(starts_with("d_"), "param", names_prefix = "d_",
                     values_to = "coupling_rate")

    fig_par <- fig_par_data %>%
        ggplot(aes(x = coupling_rate, colour = param)) +
        stat_ecdf(pad = FALSE, alpha = 0.75, size = 1.25) +
        labs(title = "ECDF of proportion of individual couplings before overall coupling",
             subtitle = sprintf("minimum %.02e iterations, sampling interval $.02e, replications = %d",
                                grid_a$run_length[1],
                                grid_a$sample_interval[1],
                                n_distinct(grid_a$c)),
             x = "proportion of samples before coupling that individual parameters met",
             y = "Fhat",
             colour = "parameter") +
        facet_wrap(~ L + lambda, ncol = 2) +
         ggsave(sprintf(fig_template,
                        sprintf("par-coupling-ecdf_axes-%s", "fixed")),
                width = 3 * 2 + 2,
                height = 3 * 2)
}
