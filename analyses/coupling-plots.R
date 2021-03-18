library("ape")
library("tidyverse")
library("fs")

source("estimators.R")
source("coupling-functions.R")
source("ipm-bounds.R")
source("tree_metrics.R")

# Set target, e.g. target <- "20210311"
target_dir <- file.path("..", target)

# Make grids of config and run settings
config_file <- file.path(target_dir, "config.R")
grids <- make_grid(config_file)
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
grid_b <- grids$grid_b

grid_d <- grid_a %>% nest(s = c(lag, c)) %>% select(-s)
# grid_d$cl <- rep(list(c("5", "6"), c("10", "9")), each = 2)
grid_d$cl <- list(c("5", "6"))
grid_d$tr <- grid_d %>%
    select(L, root_time, lambda, mu, beta) %>%
    pmap(function(...) read_nexus_file(target_dir, ...))

# Target figure and output templates and directories
fig_dir <- file.path(target_dir, "figs")
dir_create(fig_dir)
fig_template <- sprintf("%s/%%s.pdf", fig_dir)

out_dir <- file.path(target_dir, "output")

# Some useful constants
n_L <- n_distinct(grid_a$L)
n_lambda <- n_distinct(grid_a$lambda)
n_c <- n_distinct(grid_a$c)

rl_a <- grid_a$run_length[1]
si_a <- grid_a$sample_interval[1]

################################################################################
# Compute tree distances and write to file
compute_tree_distances(out_dir, grid_a)

################################################################################
# Coupling times
grid_a$tau <- get_coupling_times(out_dir, grid_a)

fig_tau <- grid_a %>%
    ggplot(aes(x = tau, colour = as.factor(lag))) +
    stat_ecdf(size = 1.5, pad = FALSE, alpha = 0.75) +
    labs(title = "ECDF of coupling time tau",
         subtitle = sprintf("%.02e iterations, replications = %d", rl_a, n_c),
         x = sprintf("tau / %d", si_a),
         y = "Fhat",
         colour = "lag")
for (scales in c("free", "fixed")) {
    fig_tau +
    facet_wrap(~ L + lambda, ncol = n_lambda, scales = scales,
               labeller = "label_both") +
    ggsave(sprintf(fig_template, sprintf("tau-ecdf_axes-%s", scales)),
           width = 3 * n_lambda + 2, height = 3 * n_L)
}

# Integral probablity metrics
iters <- seq.int(0, rl_a / si_a)
tv_data <- expand_grid(grid_a, iter = iters, tv = NA_real_)
tv_data$tv <- tv_bound_estimator(tv_data$tau,
                                 tv_data$lag / tv_data$sample_interval,
                                 tv_data$iter)

fig_tv_data <- tv_data %>%
    select(-tau) %>%
    nest(s = c(c, tv)) %>%
    mutate(tv = map_dbl(s, ~mean(.$tv)))

fig_tv <- fig_tv_data %>%
    ggplot(aes(x = iter, y = tv, colour = as.factor(lag))) +
    geom_line(size = 1.5, alpha = 0.75) +
    labs(title = "TV upper bound",
         subtitle = sprintf("%.02e iterations, replications = %d", rl_a, n_c),
         x = sprintf("iteration / %d", si_a),
         y = "d_tv")
 for (scales in c("free", "fixed")) {
     fig_tv +
     facet_wrap(~ L + lambda, ncol = n_lambda, scales = scales,
                labeller = "label_both") +
     ggsave(sprintf(fig_template, sprintf("tv_axes-%s", scales)),
            width = 3 * n_lambda + 2, height = 3 * n_L)
 }
fig_tv +
    ylim(0, 2) +
    facet_wrap(~ L, ncol = n_lambda, labeller = "label_both") +
    ggsave(sprintf(fig_template, "tv_axes-clipped"),
           width = 3 * n_lambda + 2, height = 3 * n_L)

make_w1_figure(out_dir, grid_a, grid_d, iters)

################################################################################
# Marginal histograms
make_marginal_hist(out_dir, grid_a, grid_b, "integrated_llkd", "llkd")
make_marginal_hist(out_dir, grid_a, grid_b, "root_time", "root")

################################################################################
# Estimators

make_estimator_figs(out_dir, grid_a, grid_b, NULL, "root_time", "root")
make_estimator_figs(out_dir, grid_a, grid_b, grid_d, "clade support",
                    "clade")
make_estimator_figs(out_dir, grid_a, grid_b, grid_d, "topology support",
                    "topology")

trace_estimator(out_dir, grid_a, grid_b, NULL, "root_time", "root")
trace_estimator(out_dir, grid_a, grid_b, grid_d, "clade support",
                    "clade")
trace_estimator(out_dir, grid_a, grid_b, grid_d, "topology support",
                    "topology")
