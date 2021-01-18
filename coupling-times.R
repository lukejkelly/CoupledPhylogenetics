library("ape")
library("tidyverse")
library("fs")

source("unbiased-estimator.R")

target <- "20201217"
source(file.path("..", target, "config.R"))

# TODO: Allow these to vary
root_time <- list_root_time
mu <- list_mu
beta <- list_beta

# Target templates and directories
fig_dir <- file.path("..", target, "figs")
dir_create(fig_dir)
fig_template <- sprintf("%s/%%s.pdf", fig_dir)

out_dir <- file.path("..", target, "output")
out_template <- sprintf(
    "%s/L%%d_r%.0e_l%%.0e_m%.0e_b%.0e_n%%.0e_s%%.0e-%%d%%s.%%s",
    out_dir, root_time, mu,  beta)

list_c <- seq_len(100)

#################
# Coupling times
tau <- expand.grid(L = list_L, lambda = list_lambda, c = list_c, tau = NA)

for (i in seq_len(nrow(tau))) {
    print(i / nrow(tau), digits = 3)
    # x and y sampled
    x <- read.table(
        sprintf(out_template, tau$L[i], tau$lambda[i], grid_run_length[1],
                grid_sample_interval[1], tau$c[i], "_x", "txt"),
        skip = 3)
    y <- read.table(
        sprintf(out_template, tau$L[i], tau$lambda[i], grid_run_length[1],
                grid_sample_interval[1], tau$c[i], "_y", "txt"),
        skip = 3)
    # pre-computed tree distances
    z <- read.table(
        sprintf(out_template, tau$L[i], tau$lambda[i], grid_run_length[1],
                grid_sample_interval[1], tau$c[i], "", "dist"))

    tau_p <- (last(which(x$V2[-1] != y$V2)) + 1) %>% ifelse(is.na(.), 1, .)
    tau_l <- (last(which(x$V3[-1] != y$V3)) + 1) %>% ifelse(is.na(.), 1, .)
    tau_r <- (last(which(x$V4[-1] != y$V4)) + 1) %>% ifelse(is.na(.), 1, .)
    tau_t <- (last(which(z != 0)) + 1) %>% ifelse(is.na(.), 1, .)
    tau$tau[i] <- max(c(tau_p, tau_l, tau_r, tau_t))
}

fig <- tau %>%
    ggplot(aes(x = tau, fill = as.factor(lambda), colour = NULL)) +
    stat_ecdf(aes(colour = as.factor(lambda)), size = 1.5, pad = FALSE,
              alpha = 0.75) +
    guides(colour = guide_legend(title = "lambda")) +
    # stat_bin(aes(y = ..density.., group = as.factor(lambda)), bins = 10,
    #            position = "dodge") +
    # guides(fill = guide_legend(title = "lambda"), colour = FALSE) +
    labs(title = paste("ECDF of coupling time tau as number of leaves L and",
                        "birth rate lambda vary"),
         subtitle = sprintf(paste("%.02e iterations, coupling lag = %.02e,",
                                  "replications = %d, mu = %.02e,",
                                  "beta = %.02e"),
                            grid_run_length[1], grid_sample_interval[1],
                            length(list_c), mu, beta),
         x = "tau", y = NULL)
fig +
    facet_wrap(~ L, ncol = 1, scales = "free", labeller = "label_both") +
    ggsave(sprintf(fig_template, "coupling_axes-free"), width = 8, height = 10)

fig +
    facet_wrap(~ L, ncol = 1, scales = "fixed", labeller = "label_both") +
    ggsave(sprintf(fig_template, "coupling_axes-fixed"), width = 8, height = 10)

################################################################################
# Marginal histograms
k <- 1e2
m <- 1e3
inds <- seq.int(k, m)

out <- expand.grid(L = list_L, lambda = list_lambda, c = c(0, list_c))
out_l <- matrix(nrow = nrow(out), ncol = length(inds))
out_r <- matrix(nrow = nrow(out), ncol = length(inds))

for (i in seq_len(nrow(out))) {
    print(i / nrow(out), digits = 3)
    grid_ind <- ifelse(out$c[i] == 0, 2, 1)
    x <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i],
                grid_run_length[grid_ind], grid_sample_interval[grid_ind],
                out$c[i], "_x", "txt"),
        skip = 3)
    out_l[i, ] <- x$V3[inds]
    out_r[i, ] <- x$V4[inds]
}

out_lt <- data.frame(out, logl = out_l) %>%
    pivot_longer(-c(L, lambda, c), "samp", names_prefix = "logl.",
                 values_to = "logl")
out_rt <- data.frame(out, root = out_r) %>%
    pivot_longer(-c(L, lambda, c), "samp", names_prefix = "root.",
                 values_to = "root")
out <- full_join(out_lt, out_rt) %>%
    pivot_longer(c(logl, root))

fig <- out %>%
    ggplot(aes(x = value, fill = as.factor(c), colour = NULL)) +
    stat_bin(aes(y = ..density.., group = as.factor(c)), bins = 50,
               position = "dodge", show.legend = FALSE) +
    # guides(fill = guide_legend(title = "lambda"), colour = FALSE) +
    labs(title = "Marginal distributions",
         subtitle = sprintf("x-chain samples %.02e to %.02e (interval %d)",
                            k, m, grid_sample_interval[1]))
fig +
    facet_wrap(~ L + lambda + name, ncol = 6, scales = "free",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "marginal_axes-free"), width = 60,
           height = 50, limitsize = FALSE)
fig +
    facet_wrap(~ L + lambda + name, ncol = 6, scales = "fixed",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "marginal_axes-fixed"), width = 60,
           height = 50, limitsize = FALSE)

################################################################################
# Unbiased estimators of root time
k <- 1e2
m <- 1e3

val <- expand.grid(L = list_L, lambda = list_lambda, mc = NA, bc = NA, ue = NA)
for (i in seq_len(nrow(val))) {
    print(i / nrow(val), digits = 3)
    x <- read.table(
        sprintf(out_template, val$L[i], val$lambda[i], grid_run_length[2],
                grid_sample_interval[2], 0, "_x", "txt"),
        FALSE, skip = 3)$V4

    val$mc[i] <- mean(x[seq.int(k, m) + 1])
    val$ue[i] <- val$mc[i]
}
val_data <- val %>%
    pivot_longer(-c(L, lambda))

out <- expand.grid(L = list_L, lambda = list_lambda, c = list_c, t = NA,
                   mc = NA, bc = NA, ue = NA)
for (i in seq_len(nrow(out))) {
    print(i / nrow(out), digits = 3)
    x <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i], grid_run_length[1],
                grid_sample_interval[1], out$c[i], "_x", "txt"),
        FALSE, skip = 3)$V4
    y <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i], grid_run_length[1],
                grid_sample_interval[1], out$c[i], "_y", "txt"),
        FALSE, skip = 3)$V4

    out$t[i] <- (last(which(x[-1] != y)) + 1) %>% ifelse(is.na(.), 1, .)

    f <- unbiased_estimator(x, y, k, m, out$t[i])
    out$mc[i] <- f["mc"]
    out$bc[i] <- f["bc"]
    out$ue[i] <- f["ue"]
}
fig_data <- out %>%
    select(-t) %>%
    pivot_longer(-c(L, lambda, c))
group_data <- fig_data %>%
    group_by(L, lambda, name) %>%
    summarise(x_mean = mean(value), .groups = "drop")

fig <- fig_data %>%
    ggplot(aes(x = value, fill = as.factor(lambda), colour = NULL)) +
    stat_bin(aes(y = ..density.., group = as.factor(lambda)),
             bins = 20, position = "dodge", alpha = 0.5) +
    geom_vline(mapping = aes(xintercept = x_mean, colour = as.factor(lambda),
                             linetype = as.factor(lambda)),
               alpha = 0.8, data = group_data, size = 1, linetype = 2,
               show.legend = FALSE) +
    geom_vline(mapping = aes(xintercept = value), colour = "black",
              data = val_data, size = 1.5, linetype = 3) +
    guides(fill = guide_legend(title = "lambda"), colour = FALSE) +
    labs(title = sprintf("Unbiased estimators of root time (%d paired chains)",
                         length(list_c)),
         subtitle = paste(sprintf("Samples %.02e to %.02e (interval %d),",
                                  k, m, grid_sample_interval[1]),
                          "MC average from 10x longer run"))
fig +
    facet_wrap(~ L + lambda + name, ncol = 3, scales = "free",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "root_axes-free"), width = 10,
           height = 30, limitsize = FALSE)
fig +
    facet_wrap(~ L + lambda + name, ncol = 3, scales = "fixed",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "root_axes-fixed"), width = 10,
           height = 30, limitsize = FALSE)

########
# Mean squared errors
out <- expand_grid(L = list_L, lambda = list_lambda, mc_1 = NA_real_,
                   mc_n = NA_real_, ue_1 = NA_real_, ue_n = NA_real_)
for (i in seq_len(nrow(out))) {
    L_i <- out$L[i]
    lambda_i <- out$lambda[i]
    u <- read.table(sprintf(out_template, L_i, lambda_i, grid_run_length[2],
                            grid_sample_interval[2], 0, "_x", "txt"),
                    FALSE, skip = 3)$V4
    v <- mean(u[-seq_len(k)])

    w_mc_1 <- double(length(list_c))
    w_ue_1 <- double(length(list_c))
    for (j in list_c) {
        c_j <- list_c[j]
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
    out[i, "ue_1"] <- (v - w_ue_1)^2 %>% mean()
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
                     length(list_c), k, m, grid_sample_interval[1]),
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
list_clade <- list(c(5, 6, 7, 8),
                   c(6, 7),
                   c(3, 4),
                   c(3, 4, 7, 8)) %>%
              map(as.character)
list_c <- seq.int(1, 100)

val <- expand_grid(L = list_L, lambda = list_lambda, mc = NA, bc = NA, ue = NA)
for (i in seq_len(nrow(val))) {
    print(i / nrow(val), digits = 3)
    t_x <- read.tree(sprintf(val_template, val$L[i], val$lambda[i], "_xcat",
                             "nex"),
                     skip = 8, comment.char = "#")
    v <- pluck(list_clade, which(val$L[i] == list_L))
    x <- map_lgl(t_x, is.monophyletic, v)
    val$mc[i] <- mean(x[seq.int(k + 1, m + 1)])
    val$ue[i] <- val$mc[i]
}
val_data <- val %>%
    pivot_longer(-c(L, lambda))

out <- expand_grid(L = list_L, lambda = list_lambda, c = list_c,
                   mc = NA, bc = NA, ue = NA)
for (i in seq_len(nrow(out))) {
    print(i / nrow(out), digits = 3)
    t_x <- read.tree(sprintf(out_template, out$L[i], out$lambda[i], out$c[i],
                             "_xcat", "nex"),
                     skip = 8, comment.char = "#")
    t_y <- read.tree(sprintf(out_template, out$L[i], out$lambda[i], out$c[i],
                             "_ycat", "nex"),
                     skip = 8, comment.char = "#")

    v <- pluck(list_clade, which(out$L[i] == list_L))
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
                         length(list_c)),
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
list_clade <- list(c(5, 6, 7, 8),
                   c(6, 7),
                   c(3, 4),
                   c(3, 4, 7, 8)) %>%
              map(as.character)
out <- expand_grid(L = list_Lc, lambda = list_lambdac, mc_1 = NA_real_,
                   mc_n = NA_real_, ue_1 = NA_real_, ue_n = NA_real_)
for (i in seq_len(nrow(out))) {
    L_i <- out$L[i]
    lambda_i <- out$lambda[i]
    cl_i <- pluck(list_clade, which(L_i == list_Lc))

    u <- read.tree(sprintf(out_template, L_i, lambda_i, grid_run_length[2],
                            grid_sample_interval[2], 0, "_xcat", "nex"),
                   skip = 8, comment.char = "#")
    v <- u[-seq_len(k)] %>%
        map_lgl(is.monophyletic, cl_i) %>%
        mean()

    w_mc_1 <- double(length(list_c))
    w_ue_1 <- double(length(list_c))
    for (j in list_c) {
        c_j <- list_c[j]
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
    out[i, "ue_1"] <- (v - w_ue_1)^2 %>% mean()
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
                     length(list_c), k, m, grid_sample_interval[1]),
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
