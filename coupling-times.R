library("ape")
library("tidyverse")
library("fs")

source("unbiased-estimator.R")

target <- "20201214"

root_time <- 1e3
mu <- 5e-4
beta <- 0

fig_dir <- file.path("..", target, "figs")
dir_create(fig_dir)
fig_template <- sprintf("%s/%%s.pdf", fig_dir)

out_dir <- file.path("..", target, "output")

out_template <- sprintf(
    "%s/L%%d_r%.0e_l%%.0e_m%.0e_b%.0e_n%%.0e_s%%.0e-%%d_%%s.%%s",
    out_dir, root_time, mu,  beta)

list_L <- seq.int(4, 10, 2)
list_lambda <- c(5e-02, 1e-01, 2e-01)
list_c <- seq_len(100)
grid_run_length <- c(2.5e5, 2.5e6)
grid_sample_interval <- c(1e1, 1e2)

out <- data.frame(
    expand.grid(L = list_L, lambda = list_lambda, c = list_c),
    tau_p = NA, tau_l = NA, tau_r = NA, # tau_t = NA, tau_o = NA,
    tau = NA
)

# c("sample", "log prior", "integrated llkd", "root time", "mu", "p", "lambda",
#   "kappa", "rho", "ncat", "log likelihood", "beta")

for (i in seq_len(nrow(out))) {
    print(i / nrow(out), digits = 2)
    p_x <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i], grid_run_length[1],
                grid_sample_interval[1], out$c[i], "x", "txt"),
        FALSE, skip = 3)
    p_y <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i], grid_run_length[1],
                grid_sample_interval[1], out$c[i], "y", "txt"),
        FALSE, skip = 3)

    out$tau_p[i] <- last(which(p_x$V2[-1] != p_y$V2)) + 1
    out$tau_l[i] <- last(which(p_x$V3[-1] != p_y$V3)) + 1
    out$tau_r[i] <- last(which(p_x$V4[-1] != p_y$V4)) + 1
    tau_max <- max(c(out$tau_p[i], out$tau_l[i], out$tau_r[i]))
    out$tau[i] <- ifelse(tau_max == m_inds, NA, tau_max)
}

fig <- out %>%
    ggplot(aes(x = tau, fill = as.factor(lambda), colour = NULL)) +
    stat_ecdf(aes(colour = as.factor(lambda)), size = 2, pad = FALSE) +
    guides(colour = guide_legend(title = "lambda")) +
    # stat_bin(aes(y = ..density.., group = as.factor(lambda)), bins = 10,
    #            position = "dodge") +
    # guides(fill = guide_legend(title = "lambda"), colour = FALSE) +
    labs(title = paste("Coupling time tau as number of leaves L and birth rate",
                       "lambda vary"),
         subtitle = sprintf(paste("%.02e samples, coupling lag = %.02e,",
                                  "replications = %d, mu = %.02e,",
                                  "beta = %.02e"),
                            grid_run_length[1], grid_sample_interval[1],
                            length(list_c), mu, beta),
         x = "tau")
fig +
    facet_wrap(~ L, ncol = 1, scales = "free", labeller = "label_both") +
    ggsave(sprintf(fig_template, "coupling-axes-free"), width = 8, height = 10)

fig +
    facet_wrap(~ L, ncol = 1, scales = "fixed", labeller = "label_both") +
    ggsave(sprintf(fig_template, "coupling-axes-fixed"), width = 8, height = 10)

################################################################################
# Marginal histograms
k <- 1e4
m <- 2e4
inds <- seq.int(k, m, 100)

out <- expand.grid(L = list_L, lambda = list_lambda, c = list_c)
out_l <- matrix(nrow = nrow(out), ncol = length(inds))
out_r <- matrix(nrow = nrow(out), ncol = length(inds))

for (i in seq_len(nrow(out))) {
    print(i / nrow(out), digits = 2)
    p_x <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i], grid_run_length[1],
                grid_sample_interval[1], out$c[i], "x", "txt"),
        FALSE, skip = 3)
    out_l[i, ] <- p_x$V3[inds]
    out_r[i, ] <- p_x$V4[inds]
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
         subtitle = sprintf("x-chain samples %.02e to %.02e, at interval %d",
                            k, m, grid_sample_interval[1]))
fig +
    facet_wrap(~ L + lambda + name, ncol = 6, scales = "free",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "marginal-axes-free"), width = 60,
           height = 50, limitsize = FALSE)

fig +
    facet_wrap(~ L + lambda + name, ncol = 6, scales = "fixed",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "marginal-axes-fixed"), width = 60,
           height = 50, limitsize = FALSE)

################################################################################
# Unbiased estimators
k <- 1e4
m <- 2e4
out <- expand.grid(L = list_L, lambda = list_lambda, c = list_c,
                   mc = NA, bc = NA, ue = NA)

for (i in seq_len(nrow(out))) {
    print(i / nrow(out), digits = 3)
    x <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i], grid_run_length[1],
                grid_sample_interval[1], out$c[i], "x", "txt"),
        FALSE, skip = 3)$V4
    y <- read.table(
        sprintf(out_template, out$L[i], out$lambda[i], grid_run_length[1],
                grid_sample_interval[1], out$c[i], "y", "txt"),
        FALSE, skip = 3)$V4

    t <- last(which(x[-1] != y)) + 1

    f <- unbiased_estimator(x, y, k, m, t)
    out$mc[i] <- f["mc"]
    out$bc[i] <- f["bc"]
    out$ue[i] <- f["ue"]
}

val <- expand.grid(L = list_L, lambda = list_lambda, c = 0, mc = NA)
for (i in seq_len(nrow(val))) {
    print(i / nrow(val), digits = 3)
    x <- read.table(
        sprintf(out_template, val$L[i], val$lambda[i], grid_run_length[2],
                grid_sample_interval[2], val$c[i], "x", "txt"),
        FALSE, skip = 3)$V4

    val$mc[i] <- mean(x[seq.int(k, m) + 1])
}

for (i in seq_len(nrow(out))) {
    v_i <- filter(val, L == out$L[i], lambda == out$lambda[i])$mc
    out$mc[i] <- out$mc[i] - v_i
    out$ue[i] <- out$ue[i] - v_i
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
    geom_vline(mapping = aes(xintercept = x_mean, colour = as.factor(lambda),
                             linetype = as.factor(lambda)),
               alpha = 0.8, data = group_data, size = 1, linetype = 2,
               show.legend = FALSE) +
    guides(fill = guide_legend(title = "lambda"), colour = FALSE) +
    labs(title = "Unbiased estimators",
         subtitle = sprintf("Samples %.02e to %.02e, at interval %d",
                            k, m, grid_sample_interval[1]))
fig +
    facet_wrap(~ L + name, ncol = 3, scales = "free",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "unbiased-axes-free"), width = 10,
           height = 8, limitsize = FALSE)

fig +
    facet_wrap(~ L + name, ncol = 3, scales = "fixed",
               labeller = "label_both") +
    ggsave(sprintf(fig_template, "unbiased-axes-fixed"), width = 10,
           height = 8, limitsize = FALSE)
