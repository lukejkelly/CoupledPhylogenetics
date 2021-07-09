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

target <- "20210624d"
target_dir <- file.path("..", target)

# Make grids of config and run settings
config_file <- file.path(target_dir, "config.R")
grids <- make_grid(config_file)
grid_a <- grids$grid_a

# Target figure and output templates and directories
fig_dir <- file.path(target_dir, "figs0")
# fig_dir <- file.path(target_dir, "figs")
dir_create(fig_dir)
fig_template <- sprintf("%s/%%s.pdf", fig_dir)

out_dir <- file.path(target_dir, "output0")
# out_dir <- file.path(target_dir, "output")

# Some useful constants
n_L <- n_distinct(grid_a$L)
n_lambda <- n_distinct(grid_a$lambda)
n_c <- n_distinct(grid_a$c)

rl_a <- grid_a$run_length[1]
si_a <- grid_a$sample_interval[1]

################################################################################
# Who did (s = 1) or didn't (s = 0) couple
grid_a$s <- "y"

grid_a$s[grid_a$lag == 1e3 & grid_a$c %in% c(4, 12, 17)] <- "n"
grid_a$s[grid_a$lag == 1e4 & grid_a$c %in% c(8, 9, 11, 12, 14, 18)] <- "n"

################################################################################
# Trace plots
out1 <- grid_a %>%
    filter(lag == 1e4) %>%
    nest(x = -everything(), y = -everything(), z =  -everything())
for (i in seq_len(nrow(out1))) {
    svMisc::progress(i, nrow(out1))
    out_i <- out1[i, ]

    lag_offset <- out_i$lag / out_i$sample_interval
    x <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_x")$ncat
    y <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_y")$ncat

    z <- seq.int(0, length(x) - 1)
    out1$x[[i]] <- x
    out1$y[[i]] <- c(rep(NA, lag_offset), y)
    out1$z[[i]] <- z
}

fig1 <- out1 %>%
    unnest(c(x, y, z)) %>%
    pivot_longer(c(x, y)) %>%
    ggplot(aes(x = z, y = value, group = name, colour = s)) +
    geom_line(alpha = 0.5) +
    labs(x = "sample index / 10", y = "ncat", colour = "coupled?") +
    facet_wrap(~ c, ncol = 5, scales = "free")
    ggsave(sprintf(fig_template, "ncat"), width = 10, height = 8)

##
out2 <- grid_a %>%
    filter(lag == 1e4) %>%
    nest(x = -everything(), y = -everything(), z =  -everything())
for (i in seq_len(nrow(out2))) {
    svMisc::progress(i, nrow(out2))
    out_i <- out2[i, ]

    lag_offset <- out_i$lag / out_i$sample_interval
    x1 <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_x")$log_prior
    x2 <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_x")$integrated_llkd
    y1 <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_y")$log_prior
    y2 <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_y")$integrated_llkd

    z <- seq.int(0, length(x1) - 1)
    out2$x[[i]] <- (x1 + x2) %>% tail(1e3)
    out2$y[[i]] <- c(rep(NA, lag_offset), y1 + y2) %>% tail(1e3)
    out2$z[[i]] <- z %>% tail(1e3)
}

fig2 <- out2 %>%
    unnest(c(x, y, z)) %>%
    pivot_longer(c(x, y)) %>%
    ggplot(aes(x = z, y = value, group = name, colour = as.factor(s))) +
    geom_line(alpha = 0.5) +
    labs(x = "sample index / 10", y = "log_post", colour = "coupled?") +
    facet_wrap(~ c, ncol = 5, scales = "free")
    ggsave(sprintf(fig_template, "post"), width = 10, height = 8)
