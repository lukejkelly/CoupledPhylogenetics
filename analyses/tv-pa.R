# Plotting total variation bound against acceptance proportions
# We may include this in coupling-plots at some stage

library("ape")
library("tidyverse")
library("fs")

source("estimators.R")
source("coupling-functions.R")
source("ipm-bounds.R")
source("tree-metrics.R")

# Set target, e.g. target <- "20210402"
target_dir <- file.path("..", target)

# Make grids of config and run settings
config_file <- file.path(target_dir, "config.R")
grids <- make_grid(config_file)
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
grid_b <- grids$grid_b

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

# Integral probablity metrics
iters <- seq.int(0, max(tv_data$tau, rl_a / si_a))
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
    ylim(0, 5) +
    facet_wrap(~ L + lambda, ncol = n_lambda, labeller = "label_both") +
    ggsave(sprintf(fig_template, "tv_axes-clipped"),
           width = 3 * n_lambda + 2, height = 3 * n_L)

# getting acceptance rates
lol <- get_pa_xy(out_dir, pa_data[1, ], pa_data$lag[1], pa_data$c[1],
           active_moves)

active_moves <- 1:7
pa_data <- grid_a %>% select(-tau)
pa_data$pa_xy <- map(
    seq_len(nrow(pa_data)),
    ~get_pa_xy(out_dir, pa_data[., ], pa_data$lag[.], pa_data$c[.],
               active_moves)
)
fig_pa_data <- pa_data %>%
    unnest(pa_xy) %>%
    pivot_longer(starts_with("m_"), "move", names_prefix = "m_",
                 values_to = "pa") %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval, lag,
             t, move) %>%
    summarise(pa = mean(pa), .groups = "drop")
fig_pa_data$move <- move_names(as.numeric(fig_pa_data$move))

fig_pa <- fig_pa_data %>%
    ggplot(aes(x = t, y = pa, colour = move)) +
    geom_line(alpha = 0.75) +
    labs(title = "coupled acceptance rates by move type",
         subtitle = sprintf("n = %d coupled chains", n_distinct(grid_a$c)),
         x = sprintf("iter / %d", grid_a$sample_interval[1]),
         y = "acceptance rate",
         colour = "move") +
     facet_wrap(~ L + lambda + lag, ncol = n_distinct(grid_a$lag),
                scales = "free", labeller = "label_both") +
     ggsave(sprintf(fig_template, "pa-trace"), width = 3 * 3 + 2,
            height = 3 * n_distinct(grid_a$L) * n_distinct(grid_a$lambda))
