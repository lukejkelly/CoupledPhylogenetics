## not overlapping lags
grid_a$lag[grid_a$L == 200] <- 5e5

grid_a <- grid_a %>%
    mutate(facet_label = sprintf("L: %i, lag: %.4g", L, lag))

## mu varies rather than lag
make_tau_ecdf <- function(grid_a) {
    fig_tau <- grid_a %>%
        ggplot(aes(
            x = tau - lag / sample_interval,
            linetype = factor(
                sprintf("%.4g", mu),
                levels = sprintf("%.4g", unique(sort(mu)))
            )
        )) +
        stat_ecdf(pad = FALSE, alpha = 0.75) +
        labs(
            title = "ECDF of coupling time tau",
            subtitle = sprintf(
                "minimum %.02e iterations, replications = %d",
                grid_a$run_length[1],
                n_distinct(grid_a$c)
            ),
            x = sprintf("(tau - lag) / %d", grid_a$sample_interval[1]),
            y = "ECDF",
            linetype = "mu"
        ) +
        theme_light()
    if (n_distinct(grid_a$L) > 1) {
        for (scales in c("free_x", "fixed")) {
            fig_p <- fig_tau +
                facet_wrap(
                    ~ facet_label,
                    ncol = n_distinct(grid_a$L),
                    scales = scales #,
                    # labeller = "label_both"
                )
            ggsave(
                sprintf(fig_template, sprintf("tau-ecdf_axes-%s", scales)),
                fig_p,
                width = 3 * n_distinct(grid_a$L) + 2,
                height = 3
            )
        }
    } else {
        ggsave(
            sprintf(fig_template, "tau-ecdf"),
            fig_tau,
            width = 5,
            height = 3
        )
    }
}

make_tau_eccdf <- function(grid_a) {
    fig_tau <- grid_a %>%
        ggplot(aes(
            x = tau - lag / sample_interval,
            linetype = factor(
                sprintf("%.4g", mu),
                levels = sprintf("%.4g", unique(sort(mu)))
            )
        )) +
        geom_step(aes(y = 1 - ..y..), stat = "ecdf", alpha = 0.75) +
        labs(
            title = "EECDF of coupling time tau",
            subtitle = sprintf(
                "minimum %.02e iterations, replications = %d",
                grid_a$run_length[1],
                n_distinct(grid_a$c)
            ),
            x = sprintf("(tau - lag) / %d", grid_a$sample_interval[1]),
            y = "1 - ECDF",
            linetype = "mu"
        ) +
        scale_y_log10() +
        theme_light()
    if (n_distinct(grid_a$L) > 1) {
        for (scales in c("free_x", "fixed")) {
            fig_p <- fig_tau +
                facet_wrap(
                    ~ facet_label,
                    ncol = n_distinct(grid_a$L),
                    scales = scales #,
                    # labeller = "label_both"
                )
            ggsave(
                sprintf(fig_template, sprintf("tau-eccdf_axes-%s", scales)),
                fig_p,
                width = 3 * n_distinct(grid_a$L) + 2,
                height = 3
            )
        }
    } else {
        ggsave(
            sprintf(fig_template, "tau-eccdf"),
            fig_tau,
            width = 5,
            height = 3
        )
    }
}

## tv up to last meeting time
tv_data <- grid_a %>%
    nest(n_data = c(c, tau)) %>%
    mutate(i_max = map_dbl(n_data, ~max(.$tau))) %>%
    mutate(iter = map(i_max, ~round(seq.int(0, ., length.out = 10 * n_c)))) %>%
    select(-i_max) %>%
    unnest(cols = iter) %>%
    unnest(cols = n_data) %>%
    mutate(tv = NA_real_)
tv_data$tv <- tv_bound_estimator(
    tv_data$tau,
    tv_data$lag / tv_data$sample_interval,
    tv_data$iter
)
fig_tv_data <- tv_data %>%
    select(-tau) %>%
    nest(s = c(c, tv)) %>%
    mutate(tv = map_dbl(s, ~mean(.$tv)))

fig_tv <- fig_tv_data %>%
    ggplot(aes(
        x = iter,
        y = tv,
        linetype = factor(
            sprintf("%.4g", mu),
            levels = sprintf("%.4g", unique(sort(mu)))
            # sprintf("%.4g", lag),
            # levels = sprintf("%.4g", unique(sort(lag)))
        )
    )) +
    geom_line(alpha = 0.75) +
    labs(
        title = "TV upper bound",
        subtitle = sprintf(
            "minimum %.02e iterations, replications = %d",
            grid_a$run_length[1],
            n_distinct(grid_a$c)
        ),
        x = sprintf("iteration / %d", grid_a$sample_interval[1]),
        y = "TV upper bound",
        linetype = "mu"
    ) +
    theme_light()
if (n_distinct(grid_a$L) > 1) {
    for (scales in c("free", "fixed")) {
        fig_p <- fig_tv +
            facet_wrap(
                ~ facet_label,
                ncol = n_L,
                scales = scales
                # labeller = "label_both"
            )
        ggsave(
            sprintf(fig_template, sprintf("tv_axes-%s", scales)),
            fig_p,
            width = 3 * n_distinct(grid_a$L) + 2,
            height = 3
        )
    }
} else {
    ggsave(sprintf(fig_template, "tv"), fig_tv, width = 5, height = 3)
}
