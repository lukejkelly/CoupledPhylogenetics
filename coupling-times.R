library("ape")
library("tidyverse")
library("fs")

source("estimators.R")
source("coupling-times-functions.R")

# Set target, e.g. target <- "20201217"
if (!exists("target")) {
    target <- readline("target directory name = ")
}

# Make grids of config and run settings
config_file <- file.path("..", target, "config.R")
grids <- make_grid(config_file)
# Coupled shorter runs is a, longer run is b, c is indices of coupled runs
grid_a <- grids$grid_a
grid_b <- grids$grid_b
grid_c <- grids$grid_c

# Target figure and output templates and directories
fig_dir <- file.path("..", target, "figs")
dir_create(fig_dir)
fig_template <- sprintf("%s/%%s.pdf", fig_dir)

out_dir <- file.path("..", target, "output")

# Parameters for unbiased estimators
k <- 2e2
m <- 1e3

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

################################################################################
# Marginal histograms
make_marginal_hist(out_dir, grid_a, grid_b, grid_c, k, m, "integrated_llkd",
                   "llkd")
make_marginal_hist(out_dir, grid_a, grid_b, grid_c, k, m, "root_time", "root")

################################################################################
# Estimators

# Root time
make_estimator_hist(out_dir, grid_a, grid_b, grid_c, k, m, "root_time", "root")
make_estimator_bias(out_dir, grid_a, grid_c, k, m, "root_time", "root")
make_estimator_mse(out_dir, grid_a, grid_b, grid_c, k, m, "root_time", "root")

# Clade support
grid_clade <- list(c(5, 6, 7, 8), c(6, 7), c(3, 4), c(3, 4, 7, 8)) %>%
    map(as.character) %>%
    tibble(L = c(10, 12, 6, 8), cl = .)
make_estimator_hist(out_dir, grid_a, grid_b, grid_c, k, m, "clade support",
                    "clade", grid_clade)
make_estimator_bias(out_dir, grid_a, grid_c, k, m, "clade support", "clade",
                    grid_clade)
make_estimator_mse(out_dir, grid_a, grid_b, grid_c, k, m, "clade support",
                   "clade", grid_clade)
