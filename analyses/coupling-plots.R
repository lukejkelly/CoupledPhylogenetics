library("ape")
library("ggplot2")
library("dplyr")
library("purrr")
library("fs")
library("rwty")

source("estimators.R")
source("coupling-functions.R")
source("ipm-bounds.R")
source("tree-metrics.R")
source("rwty-functions.R")

# Set target, e.g. target <- "20210428"
target_dir <- file.path("..", target)

# Make grids of config and run settings
config_file <- file.path(target_dir, "config.R")
grids <- make_grid(config_file)
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
grid_b <- grids$grid_b

grid_d <- grid_a %>% nest(s = c(lag, c)) %>% select(-s)
grid_d$cl <- rep(list(c("5", "6", "7", "8"), c("5", "6", "9", "10")), each = 1)
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
make_tau_ecdf(grid_a)

################################################################################
# Integral probablity metrics
iters <- seq.int(0, max(grid_a$tau, rl_a / si_a))
make_tv_figure(out_dir, grid_a, iters)
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

################################################################################
# Are We There Yet?
grid_e <- get_rwty_output(out_dir, grid_a)

fig_awty_data <- grid_e %>%
    select(-c(c, tau)) %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval, lag)
fig_awty <- fig_awty_data %>%
    group_map(~makeplot.asdsf(.x$rwty, 1, win_size(.y), 0.1)$asdsf.plot +
               labs(title = sprintf("L = %d, lambda = %g, lag = %g", .y$L,
                                    .y$lambda, .y$lag),
                    subtitle = sprintf("window size = %d", win_size(.y)),
                    x = sprintf("iteration / %d", .y$sample_interval)))

gridExtra::grid.arrange(grobs = fig_awty,
                        ncol = n_distinct(fig_awty_data$lag),
                        nrow = prod(n_distinct(fig_awty_data$L),
                                    n_distinct(fig_awty_data$lambda))) %>%
    ggsave(sprintf(fig_template, "asdsf"),
           plot = .,
           width = 3 * n_distinct(fig_awty_data$lag) + 2,
           height = 3 * prod(n_distinct(fig_awty_data$L),
                             n_distinct(fig_awty_data$lambda)))
