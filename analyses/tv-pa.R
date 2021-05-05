# Plotting total variation bound against acceptance proportions
# We may include this in coupling-plots at some stage
# Requires samples at subsample rate 1

library("ape")
library("tidyr")
library("dplyr")
library("purrr")
library("fs")
library("ggplot2")

source("estimators.R")
source("coupling-functions.R")
source("ipm-bounds.R")
source("tree-metrics.R")

# Set target, e.g. target <- "20210428"
target_dir <- file.path("..", target)

# Make grids of config and run settings
config_file <- file.path(target_dir, "config.R")
grids <- make_grid(config_file)
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
if (any(grid_a$sample_interval != 1)) {
    stop("Requires chains with every sample")
}
# grid_b <- grids$grid_b

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
iters <- seq.int(0, max(grid_a$tau, rl_a / si_a, na.rm = TRUE))
make_tv_figure(out_dir, grid_a, iters)

# getting acceptance rates
active_moves <- c(1:8, 12)
pa_data <- grid_a # %>%
    # filter(!is.na(tau)) %>%
    # select(-tau)
pa_data$pa_xy <- map(
    seq_len(nrow(pa_data)),
    ~get_pa_xy(out_dir, pa_data[., ], pa_data$lag[.], pa_data$c[.],
               active_moves)
)
pa_data <- pa_data %>% select(-tau)

# plot average acceptance rates
fig_rate_data <- pa_data %>%
    unnest(pa_xy) %>%
    pivot_longer(starts_with("m_"), "move", names_prefix = "m_",
                 values_to = "pa") %>%
    filter(!is.na(pa)) %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval, lag,
             move) %>%
    summarise(pa = mean(pa), .groups = "drop")
fig_rate_data$move <- move_names(as.numeric(fig_rate_data$move))

fig_rate <- fig_rate_data %>%
    ggplot(aes(x = move, y = pa, fill = move)) +
    geom_col() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = sprintf("average acceptance rate by move across n = %d coupled chains",
                         n_distinct(grid_a$c)),
         x = sprintf("iter / %d", grid_a$sample_interval[1]),
         y = "acceptance rate",
         colour = "move") +
     facet_wrap(~ L + lambda + lag, ncol = n_distinct(grid_a$lag),
                scales = "free", labeller = "label_both") +
     ggsave(sprintf(fig_template, "pa-hist"),
            width = 3 * n_distinct(grid_a$lag) + 2,
            height = 3 * n_distinct(grid_a$L) * n_distinct(grid_a$lambda))

# same plot when lambda dummy indexes the experiments
lambda_map <- c("vary mu = N, clade = N", "vary mu = Y, clade = N",
                "vary mu = N, clade = Y", "vary mu = Y, clade = Y")
fig_rate_data$lambda2 <- lambda_map[fig_rate_data$lambda]

fig_rate <- fig_rate_data %>%
    ggplot(aes(x = move, y = pa, fill = move)) +
    geom_col() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = sprintf("average acceptance rate by move across n = %d coupled chains",
                         n_distinct(grid_a$c)),
         x = sprintf("iter / %d", grid_a$sample_interval[1]),
         y = "acceptance rate",
         colour = "move") +
     facet_wrap(~ lambda2, ncol = length(lambda_map) / 2, scales = "fixed",
                labeller = "label_both") +
     ggsave(sprintf(fig_template, "pa-hist"),
            width = 3 * length(lambda_map) / 2 + 2,
            height = 3 * length(lambda_map) / 2)

# plot final move (only works for subsample = 1 and chains stop at coupling)
fig_final_data <- pa_data %>%
    unnest(pa_xy) %>%
    pivot_longer(starts_with("m_"), "move", names_prefix = "m_",
                 values_to = "pa") %>%
    filter(!is.na(pa)) %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval, lag,
             c) %>%
    slice_tail() %>%
    ungroup()
fig_final_data$move <- move_names(as.numeric(fig_final_data$move))

fig_final <- fig_final_data %>%
    ggplot(aes(x = move, fill = move)) +
    geom_bar() +
    # geom_bar(aes(y = (..count..) / sum(..count..))) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = "last move before coupling",
         subtitle = sprintf("across n = %d coupled chains",
                            n_distinct(grid_a$c)),
         y = "count",
         colour = "move") +
    facet_wrap(~ L + lambda + lag, ncol = n_distinct(grid_a$lag),
               scales = "free", labeller = "label_both") +
    ggsave(sprintf(fig_template, "final-hist"),
           width = 3 * n_distinct(grid_a$lag) + 2,
           height = 3 * n_distinct(grid_a$L) * n_distinct(grid_a$lambda))

# getting distances
pa_data$dist <- map(
    seq_len(nrow(pa_data)),
    ~read.table(make_file_name_(out_dir, pa_data[., ], pa_data$lag[.],
                                pa_data$c[.], "", "dist")) %>%
        pull(V1) %>%
        diff()
)

fig_dist_data <- pa_data %>%
    unnest(c(pa_xy, dist)) %>%
    pivot_longer(starts_with("m_"), "move", names_prefix = "m_",
                 values_to = "pa") %>%
    filter(!is.na(pa)) %>%
    filter(dist != 0) %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval, lag,
             c) %>%
    slice_tail(prop = 0.5) %>%
    ungroup()
fig_dist_data$move <- move_names(as.numeric(fig_dist_data$move))
fig_dist_data$pa <- ifelse(fig_dist_data$pa == 1, "coupled", "uncoupled")

fig_dist <- fig_dist_data %>%
    ggplot(aes(x = move, y = dist, colour = pa)) +
    geom_boxplot(alpha = 0.5) +
    labs(title = sprintf("change in score distance by move across n = %d coupled chains",
                         n_distinct(grid_a$c)),
         subtitle = "first half of chain ignored",
         x = "move",
         y = "change in score distance",
         colour = element_blank()) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
     facet_wrap(~ L + lambda + lag, ncol = n_distinct(grid_a$lag),
                scales = "free", labeller = "label_both") +
     ggsave(sprintf(fig_template, "dist-box"), width = 3 * 3 + 2,
            height = 3 * n_distinct(grid_a$L) * n_distinct(grid_a$lambda))

# plot change in posterior
get_diff_post <- function(out_dir, grid, lag, c) {

}
pa_data$diff_post <- map(
    seq_len(nrow(pa_data)),
    ~get_diff_post(make_file_name_(out_dir, pa_data[., ], pa_data$lag[.],
                                   pa_data$c[.], "", "dist"))
)

fig_post_data <- pa_data %>%
    unnest(c(pa_xy, dist)) %>%
    pivot_longer(starts_with("m_"), "move", names_prefix = "m_",
                 values_to = "pa") %>%
    filter(!is.na(pa)) %>%
    filter(dist != 0) %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval, lag,
             c) %>%
    slice_tail(prop = 0.5) %>%
    ungroup()
fig_post_data$move <- move_names(as.numeric(fig_post_data$move))
fig_post_data$pa <- ifelse(fig_post_data$pa == 1, "both", "only one")

fig_post <- fig_post_data %>%
    ggplot(aes(x = move, y = dist, colour = pa)) +
    geom_violin(alpha = 0.5) +
    labs(title = sprintf("change in score distance by move across n = %d coupled chains",
                         n_distinct(grid_a$c)),
         subtitle = "first half of chain ignored",
         x = "move",
         y = "change in score distance",
         colour = element_blank()) +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
     facet_wrap(~ L + lambda + lag, ncol = n_distinct(grid_a$lag),
                scales = "free", labeller = "label_both") +
     ggsave(sprintf(fig_template, "dist-violin"), width = 3 * 3 + 2,
            height = 3 * n_distinct(grid_a$L) * n_distinct(grid_a$lambda))
