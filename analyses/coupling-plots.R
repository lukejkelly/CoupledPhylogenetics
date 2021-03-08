library("ape")
library("tidyverse")
library("fs")

source("estimators.R")
source("coupling-functions.R")
source("ipm-bounds.R")

# Set target, e.g. target <- "20210113"
target_dir <- file.path("..", "CoupledPhylogeneticFiles", target)

# Make grids of config and run settings
config_file <- file.path(target_dir, "config.R")
grids <- make_grid(config_file)
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
grid_b <- grids$grid_b
grid_c <- grids$grid_c

grid_d <- grid_a
grid_d$cl <- rep(list(c("5", "6"), c("10", "9")), each = 2)
grid_d$tr <- grid_a %>%
    select(L, root_time, lambda, mu, beta) %>%
    pmap(function(...) read_nexus_file(target_dir, ...))

# Target figure and output templates and directories
fig_dir <- file.path(target_dir, "figs")
dir_create(fig_dir)
fig_template <- sprintf("%s/%%s.pdf", fig_dir)

out_dir <- file.path(target_dir, "output")

################################################################################
# Coupling times
fig_tau_data <- get_coupling_times(out_dir, grid_a, grid_c)

fig_tau <- fig_tau_data %>%
    ggplot(aes(x = tau, fill = as.factor(lambda), colour = NULL)) +
    stat_ecdf(aes(colour = as.factor(lambda)), size = 1.5, pad = FALSE,
              alpha = 0.75) +
    guides(colour = guide_legend(title = "lambda")) +
    labs(title = "ECDF of coupling time tau as number of leaves L and birth rate lambda vary",
         subtitle = sprintf(
             "%.02e iterations, coupling lag = %.02e, replications = %d",
             grid_a$run_length[1], grid_a$sample_interval[1], length(grid_c)),
         x = "tau",
         y = "Fhat")
for (scales in c("free", "fixed")) {
    fig_tau +
    facet_wrap(~ L, ncol = 1, scales = scales, labeller = "label_both") +
    ggsave(sprintf(fig_template, sprintf("tau-ecdf_axes-%s", scales)),
           width = 8, height = 10)
}

# Integral probablity metrics
iters <- seq.int(0, grid_a$run_length[1] / grid_a$sample_interval[1])
lag <- 1
tv_data <- expand_grid(fig_tau_data, iter = iters, tv = NA_integer_)
tv_data$tv <- tv_bound_estimator(tv_data$tau, lag, tv_data$iter)

fig_tv_data <- tv_data %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval,
             iter) %>%
    summarise(tv = mean(tv), .groups = "drop")

fig_tv <- fig_tv_data %>%
    ggplot(aes(x = iter, y = tv, colour = as.factor(lambda))) +
    geom_line(size = 1.5, alpha = 0.75) +
    guides(colour = guide_legend(title = "lambda")) +
    labs(title = "TV upper bound as number of leaves L and birth rate lambda vary",
         subtitle = sprintf(
             "%.02e iterations, coupling lag = %.02e, replications = %d",
             grid_a$run_length[1], grid_a$sample_interval[1], length(grid_c)),
         x = "iteration",
         y = "d_tv")
 for (scales in c("free", "fixed")) {
     fig_tv +
     facet_wrap(~ L, ncol = 1, scales = scales, labeller = "label_both") +
     ggsave(sprintf(fig_template, sprintf("tv_axes-%s", scales)),
            width = 8, height = 10)
 }
fig_tv +
    ylim(0, 2) +
    facet_wrap(~ L, ncol = 1) +
    ggsave(sprintf(fig_template, "tv_axes-clipped"), width = 8, height = 10)

make_w1_figure(out_dir, fig_tau_data, grid_d, 1, iters)

################################################################################
# Marginal histograms
make_marginal_hist(out_dir, grid_a, grid_b, grid_c, "integrated_llkd", "llkd")
make_marginal_hist(out_dir, grid_a, grid_b, grid_c, "root_time", "root")

################################################################################
# Estimators

make_estimator_figs(out_dir, grid_a, grid_b, grid_c, NULL, "root_time", "root")
make_estimator_figs(out_dir, grid_a, grid_b, grid_c, grid_d, "clade support",
                    "clade")
make_estimator_figs(out_dir, grid_a, grid_b, grid_c, grid_d, "topology support",
                    "topology")
