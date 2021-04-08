library("ape")
library("tidyverse")
library("fs")
library("rwty")

source("estimators.R")
source("coupling-functions.R")
source("ipm-bounds.R")
source("tree-metrics.R")
source("rwty-functions.R")

# Set target, e.g. target <- "20210402"
target_dir <- file.path("..", target)

# Make grids of config and run settings
config_file <- file.path(target_dir, "config.R")
grids <- make_grid(config_file)
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
grid_b <- grids$grid_b

grid_d <- grid_a %>% nest(s = c(lag, c)) %>% select(-s)
grid_d$cl <- rep(list(c("5", "6"), c("10", "9")), each = 2)
# grid_d$cl <- list(c("5", "6"))
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
grid_e <- get_rwty_trees(out_dir, grid_a)

fig_awty_data <- grid_e %>%
    select(-c) %>%
    group_by(L, root_time, lambda, mu, beta)
fig_awty_keys <- group_keys(fig_awty_data)
fig_awty <- fig_awty_data %>%
    group_map(~makeplot.asdsf(.x$rwty)$asdsf.plot +
               labs(title = sprintf("L = %d : lambda = %f", .y$L, .y$lambda)))

for (i in seq_len(nrow(fig_awty_keys))) {
    fig_awty[[i]] <- fig_awty[[i]]$asdsf.plot +
        labs(title = sprintf("L = %d : lambda = %f", fig_awty_keys$L[i],
                             fig_awty_keys$lambda[i]))
}
gridExtra::grid.arrange(grobs = fig_awty, nrow = n_distinct(fig_awty_keys$L),
                        ncol = n_distinct(fig_awty_keys$lambda)) %>%
    ggsave(sprintf(fig_template, "asdsf"),
           plot = .,
           width = 3 * n_distinct(fig_awty_keys$lambda) + 2,
           height = 3 * n_distinct(fig_awty_keys$L))
