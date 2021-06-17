# Checking how long xi terms take to couple, originally written for
#    20210615a, maximal coupling of scaling one xi and CRN scaling all
#    20210615b, CRN coupling of sampling uniformly for one or all xi
#  with the topology and everything else but the node times fixed
# There was no discernible difference in the overall coupling times so we're
# checking if xi coupled earlier and whether there was a difference

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

compare_xi_times <- function(target) {
    target_dir <- file.path("..", target)

    # Make grids of config and run settings
    config_file <- file.path(target_dir, "config.R")
    grid_a <- make_grid(config_file)$grid_a
    out_dir <- file.path(target_dir, "output")

    d_xi <- compute_xi_tau(out_dir, grid_a)
    d_bl <- compute_tree_score_tau(out_dir, grid_a)

    return(list(xi = d_xi, bl = d_bl))
}

compute_xi_tau <- function(out_dir, grid_a) {
    d <- double(nrow(grid_a))
    for (i in seq_len(nrow(grid_a))) {
        svMisc::progress(i, nrow(grid_a))
        grid_a_i <- grid_a[i, ]

        x <- get_xi_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_x")
        y <- get_xi_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_y")
        lag_offset <- grid_a_i$lag / grid_a_i$sample_interval

        z <- numeric(grid_a_i$L)
        for (j in seq_along(z)) {
            z[j] <- get_tau(pull(x[, j]), pull(y[, j]), lag_offset)
        }
        d[i] <- max(z)
    }
    message("xi scores computed")
    return(d)
}

get_xi_ <- function(file_dir, grid, lag, c, d) {
    file_name <- make_file_name_(file_dir, grid, lag, c, paste0(d, "XI"), "txt")
    xi <- readr::read_table(file_name, FALSE,
                            sprintf("-%s-", paste0(rep("d", grid$L),
                                                   collapse = "")),
                            skip = 6)
    return(xi)
}

compute_tree_score_tau <- function(out_dir, grid_a) {
    # Return coupling time estimate from topological score distance
    d <- double(nrow(grid_a))
    for (i in seq_len(nrow(grid_a))) {
        svMisc::progress(i, nrow(grid_a))
        grid_a_i <- grid_a[i, ]

        x <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_x")
        y <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_y")

        lag_offset <- grid_a_i$lag / grid_a_i$sample_interval

        d[i] <- purrr::map2_dbl(x[-seq_len(lag_offset)], y, ape::dist.topo,
                                "score") %>%
                purrr::map_lgl(~ . != 0) %>%
                get_tau_(lag_offset)
    }
    message("tree distances computed")
    return(d)
}

# compute_tree_match_tau <- function(out_dir, grid_a) {
#     # Compute proportion of matching branch lengths assuming topology fixed
#     # then return marginal coupling time
#     d <- double(nrow(grid_a))
#     for (i in seq_len(nrow(grid_a))) {
#         svMisc::progress(i, nrow(grid_a))
#         grid_a_i <- grid_a[i, ]
#
#         x <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_x")
#         y <- get_trees_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_y")
#
#         lag_offset <- grid_a_i$lag / grid_a_i$sample_interval
#
#         d[i] <- purrr::map2_dbl(x[-seq_len(lag_offset)], y,
#                                 ~ mean(.x$edge.length == .y$edge.length)) %>%
#                 purrr::map_lgl(~ . != 1) %>%
#                 get_tau_(lag_offset)
#     }
#     message("tree distances computed")
#     return(d)
# }

##
targets <- c("20210615a", "20210615b")
out <- map(targets, compare_xi_times)

fig_data <- bind_rows(
    bind_cols(i = "0", xi = out[[1]]$xi, bl = out[[1]]$bl),
    bind_cols(i = "1", xi = out[[2]]$xi, bl = out[[2]]$bl)
)
fig_data$max <- pmax(fig_data$xi, fig_data$bl)

fig <- fig_data %>%
    ggplot(aes(x = xi, y = max, colour = i)) +
    geom_point() +
    geom_abline() +
    ggsave("../figs/xi.pdf")
# rsync -ahHSv ceremade:~/CoupledPhylogenetics/figs/xi.pdf ~/Workspace/CoupledPhylogenetics/figs
