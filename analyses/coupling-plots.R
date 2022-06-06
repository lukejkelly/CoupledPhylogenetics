library("ape")
library("ggplot2")
library("dplyr")
library("purrr")
library("fs")
library("rwty") # modified to implement makeplot.asdsf.mb

# Run from current experiment folder at same level in hierarchy as analyses/
source_analysis_file <- function(s) source(file.path("..", "analyses", s))
source_analysis_file("estimators.R")
source_analysis_file("coupling-functions.R")
source_analysis_file("ipm-bounds.R")
# source_analysis_file("tree-metrics.R")
source_analysis_file("rwty-functions.R")

# Make grids of config and run settings
grids <- make_grid("config.R")
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
grid_b <- grids$grid_b

# grid_d <- grid_b # grid_a %>% nest(s = c(lag, c)) %>% select(-s)
# # grid_d$cl <- rep(list(c("5", "6", "7", "8"), c("5", "6", "9", "10")), each = 1)
# # grid_d$tr <- grid_d %>%
# #     select(L, root_time, lambda, mu, beta) %>%
# #     pmap(function(...) read_nexus_file(target_dir, ...))
# grid_d$cl <- NA
# grid_d$tr <- NA

# Target figure and output templates and directories
out_dir <- "output"
fig_dir <- "figs"
dir_create(fig_dir)
fig_template <- sprintf("%s/%%s.pdf", fig_dir)

# Some useful constants
n_L <- n_distinct(grid_a$L)
n_lambda <- n_distinct(grid_a$lambda)
n_c <- n_distinct(grid_a$c)

rl_a <- grid_a$run_length[1]
si_a <- grid_a$sample_interval[1]

################################################################################
# # Compute tree distances and write to file
# compute_tree_distances(out_dir, grid_a)

################################################################################
# Coupling times
grid_a$tau <- get_coupling_times(out_dir, grid_a)
make_tau_ecdf(grid_a)
make_tau_eccdf(grid_a)

################################################################################
# Integral probablity metrics
# iters <- seq.int(0, max(grid_a$tau, rl_a / si_a))
# iters <- seq.int(0, max(grid_a$tau - grid_a$lag / grid_a$sample_interval) + 1)
iters <- seq.int(0, max(grid_a$lag / grid_a$sample_interval))
make_tv_figure(out_dir, grid_a, iters)

################################################################################
# Marginal histograms
make_marginal_hist(out_dir, grid_a, grid_b, "log_prior", "prior")
make_marginal_hist(out_dir, grid_a, grid_b, "integrated_llkd", "llkd")
make_marginal_hist(out_dir, grid_a, grid_b, "root_time", "root")
make_marginal_hist(out_dir, grid_a, grid_b, "ncat", "ncat")
make_marginal_hist(out_dir, grid_a, grid_b, "mu", "mu")
make_marginal_hist(out_dir, grid_a, grid_b, "beta", "beta")

################################################################################
# # Estimators
# make_estimator_figs(out_dir, grid_a, grid_b, NULL, "root_time", "root")
# # make_estimator_figs(out_dir, grid_a, grid_b, grid_d, "clade support",
# #                     "clade")
# make_estimator_figs(out_dir, grid_a, grid_b, grid_d, "topology support",
#                     "topology")
#
# trace_estimator(out_dir, grid_a, grid_b, NULL, "root_time", "root")
# # trace_estimator(out_dir, grid_a, grid_b, grid_d, "clade support",
# #                     "clade")
# trace_estimator(out_dir, grid_a, grid_b, grid_d, "topology support",
#                 "topology")

################################################################################
# ASDSF

# read trees and parameters
if (all(grid_a$run_length > 0)) {
    grid_e <- get_rwty_output(out_dir, grid_a)
} else {
    grid_e <- grid_a %>%
        filter(lag == max(lag)) %>%
        mutate(run_length = lag) %>%
        get_rwty_output(out_dir, .)
}
fig_rwty_data <- grid_e %>%
    select(-c(c, tau, lag)) %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval)

# calculate ASDSF for each
burnin <- 0
window.lookback <- 0.75
window.number <- 10
min.freq <- 0.1

fig_rwty <- fig_rwty_data %>%
    group_map(
        ~makeplot.asdsf.mb(
            .x$rwty, burnin, window.lookback, window.number, min.freq
        )$asdsf.plot +
        labs(
            title = sprintf("window lookback = %g", window.lookback),
            x = sprintf("iteration / %d", .y$sample_interval),
            y = "ASDSF"
        )
    )

# package plot
fig_rwty %>%
    gridExtra::grid.arrange(
        grobs = .,
        ncol = n_distinct(fig_rwty_data$L),
        nrow = 1
    ) %>%
    ggsave(
        sprintf(fig_template, "asdsf-rwty"),
        plot = .,
        width = 3 * n_distinct(fig_rwty_data$L) + 2,
        height = 3
    )

# modified plot
my_asdsf_plot(grid_e, fig_rwty)
