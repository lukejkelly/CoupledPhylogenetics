library("ape")
library("tidyverse")
library("fs")

source("estimators.R")
source("coupling-times-functions")

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
k <- 1e2
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

make_estimator_hist(out_dir, grid_a, grid_b, grid_c, k, m, "root_time", "root")
make_estimator_bias(out_dir, grid_a, grid_c, k, m, "root_time", "root")
make_estimator_mse(out_dir, grid_a, grid_b, grid_c, k, m, "root_time", "root")


########
# Mean squared errors
out <- expand_grid(L = list_L, lambda = list_lambda, mc_1 = NA_real_,
                   mc_n = NA_real_, ue_n = NA_real_)
for (i in seq_len(nrow(out))) {
    L_i <- out$L[i]
    lambda_i <- out$lambda[i]
    u <- read.table(sprintf(out_template, L_i, lambda_i, grid_run_length[2],
                            grid_sample_interval[2], 0, "_x", "txt"),
                    FALSE, skip = 3)$V4
    v <- mean(u[-seq_len(k)])

    w_mc_1 <- double(length(grid_c))
    w_ue_1 <- double(length(grid_c))
    for (j in grid_c) {
        c_j <- grid_c[j]
        x <- read.table(sprintf(out_template, L_i, lambda_i, grid_run_length[1],
                                grid_sample_interval[1], c_j, "_x", "txt"),
                        FALSE, skip = 3)$V4
        y <- read.table(sprintf(out_template, L_i, lambda_i, grid_run_length[1],
                                grid_sample_interval[1], c_j, "_y", "txt"),
                        FALSE, skip = 3)$V4

        t <- (last(which(x[-1] != y)) + 1) %>% ifelse(is.na(.), 1, .)
        f <- unbiased_estimator(x, y, k, m, t)

        w_mc_1[j] <- f["mc"]
        w_ue_1[j] <- f["ue"]
    }
    w_mc_n <- mean(w_mc_1)
    w_ue_n <- mean(w_ue_1)

    out[i, "mc_1"] <- (v - w_mc_1)^2 %>% mean()
    out[i, "mc_n"] <- (v - w_mc_n)^2 %>% mean()
    out[i, "ue_n"] <- (v - w_ue_n)^2 %>% mean()
}
fig_data <- out %>%
    pivot_longer(-c(L, lambda), "estimator", values_to = "mse")

fig <- fig_data %>%
    ggplot(aes(x = estimator, y = mse)) +
    geom_col() +
    labs(title = "Monte Carlo (mc) and unbiased (ue) estimators of root time ",
         subtitle = paste(
             sprintf("%d coupled chains, samples %g to %g (interval %d);",
                     length(grid_c), k, m, grid_sample_interval[1]),
             "MSE uses MC average from 10x longer run;",
             "*_n is average of *_1 which uses a single chain (pair)"))
fig +
    facet_wrap(~ L + lambda, ncol = 3, scales = "free",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "root-mse_axes-free_scale-linear"), width = 12,
           height = 12, limitsize = FALSE)
fig +
    facet_wrap(~ L + lambda, ncol = 3, scales = "free",
               labeller = "label_both") +
    scale_y_continuous(trans = "log1p") +
    ggsave(sprintf(fig_template, "root-mse_axes-free_scale-log"), width = 12,
           height = 12, limitsize = FALSE)


################################################################################
# Unbiased estimators of clade probabilities

out_template <- sprintf(
    "%s/L%%d_r%.0e_l%%.0e_m%.0e_b%.0e_n%.0e_s%.0e-%%d%%s.%%s",
    out_dir, root_time, mu,  beta,
    grid_run_length[1], grid_sample_interval[1])
val_template <- sprintf(
    "%s/L%%d_r%.0e_l%%.0e_m%.0e_b%.0e_n%.0e_s%.0e-%d%%s.%%s",
    out_dir, root_time, mu,  beta,
    grid_run_length[2], grid_sample_interval[2], 0)

list_L <- c(10, 12, 6, 8)
list_lambda <- c(1e-01, 2e-1, 5e-2)
grid_clade <- list(c(5, 6, 7, 8),
                   c(6, 7),
                   c(3, 4),
                   c(3, 4, 7, 8)) %>%
              map(as.character)
grid_c <- seq.int(1, 100)

val <- expand_grid(L = list_L, lambda = list_lambda, mc = NA, bc = NA, ue = NA)
for (i in seq_len(nrow(val))) {
    print(i / nrow(val), digits = 3)
    t_x <- read.tree(sprintf(val_template, val$L[i], val$lambda[i], "_xcat",
                             "nex"),
                     skip = 8, comment.char = "#")
    v <- pluck(grid_clade, which(val$L[i] == list_L))
    x <- map_lgl(t_x, is.monophyletic, v)
    val$mc[i] <- mean(x[seq.int(k + 1, m + 1)])
    val$ue[i] <- val$mc[i]
}
val_data <- val %>%
    pivot_longer(-c(L, lambda))

out <- expand_grid(L = list_L, lambda = list_lambda, c = grid_c,
                   mc = NA, bc = NA, ue = NA)
for (i in seq_len(nrow(out))) {
    print(i / nrow(out), digits = 3)
    t_x <- read.tree(sprintf(out_template, out$L[i], out$lambda[i], out$c[i],
                             "_xcat", "nex"),
                     skip = 8, comment.char = "#")
    t_y <- read.tree(sprintf(out_template, out$L[i], out$lambda[i], out$c[i],
                             "_ycat", "nex"),
                     skip = 8, comment.char = "#")

    v <- pluck(grid_clade, which(out$L[i] == list_L))
    x <- map_lgl(t_x, is.monophyletic, v)
    y <- map_lgl(t_y, is.monophyletic, v)
    t <- (last(which(x[-1] != y)) + 1) %>% ifelse(is.na(.), 1, .)

    f <- unbiased_estimator(x, y, k, m, t)

    out$mc[i] <- f["mc"]
    out$bc[i] <- f["bc"]
    out$ue[i] <- f["ue"]
}

fig_data <- out %>%
    pivot_longer(-c(L, lambda, c))
group_data <- fig_data %>%
    group_by(L, lambda, name) %>%
    summarise(x_mean = mean(value), .groups = "drop")

fig <- fig_data %>%
    ggplot(aes(x = value, fill = as.factor(lambda), colour = NULL)) +
    stat_bin(aes(y = ..density.., group = as.factor(lambda)),
             bins = 20, position = "dodge", alpha = 0.5) +
    geom_vline(mapping = aes(xintercept = x_mean, colour = as.factor(lambda)),
               alpha = 0.8, data = group_data, size = 1, linetype = 2,
               show.legend = FALSE) +
    geom_vline(mapping = aes(xintercept = value), colour = "black",
               data = val_data, size = 1.5, linetype = 3) +
    guides(fill = guide_legend(title = "lambda"), colour = FALSE) +
    labs(title = sprintf(paste("Unbiased estimators of clade support",
                               "(%d paired chains)"),
                         length(grid_c)),
         subtitle = paste(sprintf("Samples %.02e to %.02e (interval %d),",
                                  k, m, grid_sample_interval[1]),
                          "MC average from 10x longer run in black,",
                          "Sample average dashed"))
fig +
    facet_wrap(~ L + lambda + name, ncol = 3, scales = "free",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "clades_axes-free"), width = 10,
           height = 30, limitsize = FALSE)

fig +
    facet_wrap(~ L + lambda + name, ncol = 3, scales = "fixed",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "clades_axes-fixed"), width = 10,
           height = 30, limitsize = FALSE)

################
# Mean squared errors

list_Lc <- c(10, 12, 6, 8)
list_lambdac <- c(1e-01, 2e-1, 5e-2)
grid_clade <- list(c(5, 6, 7, 8),
                   c(6, 7),
                   c(3, 4),
                   c(3, 4, 7, 8)) %>%
              map(as.character)
out <- expand_grid(L = list_Lc, lambda = list_lambdac, mc_1 = NA_real_,
                   mc_n = NA_real_, ue_n = NA_real_)
for (i in seq_len(nrow(out))) {
    L_i <- out$L[i]
    lambda_i <- out$lambda[i]
    cl_i <- pluck(grid_clade, which(L_i == list_Lc))

    u <- read.tree(sprintf(out_template, L_i, lambda_i, grid_run_length[2],
                            grid_sample_interval[2], 0, "_xcat", "nex"),
                   skip = 8, comment.char = "#")
    v <- u[-seq_len(k)] %>%
        map_lgl(is.monophyletic, cl_i) %>%
        mean()

    w_mc_1 <- double(length(grid_c))
    w_ue_1 <- double(length(grid_c))
    for (j in grid_c) {
        c_j <- grid_c[j]
        t_x <- read.tree(sprintf(out_template, L_i, lambda_i,
                                 grid_run_length[1], grid_sample_interval[1],
                                 c_j, "_xcat", "nex"),
                         skip = 8, comment.char = "#")
        t_y <- read.tree(sprintf(out_template, L_i, lambda_i,
                                 grid_run_length[1], grid_sample_interval[1],
                                 c_j, "_ycat", "nex"),
                         skip = 8, comment.char = "#")


        x <- map_lgl(t_x, is.monophyletic, cl_i)
        y <- map_lgl(t_y, is.monophyletic, cl_i)

        t <- (last(which(x[-1] != y)) + 1) %>% ifelse(is.na(.), 1, .)
        f <- unbiased_estimator(x, y, k, m, t)

        w_mc_1[j] <- f["mc"]
        w_ue_1[j] <- f["ue"]
    }
    w_mc_n <- mean(w_mc_1)
    w_ue_n <- mean(w_ue_1)

    out[i, "mc_1"] <- (v - w_mc_1)^2 %>% mean()
    out[i, "mc_n"] <- (v - w_mc_n)^2 %>% mean()
    out[i, "ue_n"] <- (v - w_ue_n)^2 %>% mean()
}
fig_data <- out %>%
    pivot_longer(-c(L, lambda), "estimator", values_to = "mse")

fig <- fig_data %>%
    ggplot(aes(x = estimator, y = mse)) +
    geom_col() +
    labs(title = "Monte Carlo (mc) and unbiased (ue) estimators of clade probabilities",
         subtitle = paste(
             sprintf("%d coupled chains, samples %g to %g (interval %d);",
                     length(grid_c), k, m, grid_sample_interval[1]),
             "MSE uses MC average from 10x longer run;",
             "*_n is average of *_1 which uses a single chain (pair)"))
fig +
    facet_wrap(~ L + lambda, ncol = 3, scales = "free",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "clade-mse_axes-free_scale-linear"),
           width = 12, height = 12, limitsize = FALSE)
fig +
    facet_wrap(~ L + lambda, ncol = 3, scales = "free",
               labeller = "label_both") +
    scale_y_continuous(trans = "log1p") +
    ggsave(sprintf(fig_template, "clade-mse_axes-free_scale-log"),
           width = 12, height = 12, limitsize = FALSE)
