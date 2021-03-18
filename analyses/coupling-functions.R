library("tidyverse")
library("ape")

# Functions to construct file names and paths
get_file_template_a <- function() {
    file_template_a <- "L%d_r%e_l%e_m%e_b%e_l%e-%d%s.%s"
    return(file_template_a)
}
get_file_template_b <- function() {
    file_template_b <- "L%d_r%e_l%e_m%e_b%e.%s"
    return(file_template_b)
}

make_file_name <- function(file_dir, L, root_time, lambda, mu, beta,
                           lag = NULL, c, d, e) {
    # Run index is c, decorator d = "x" or "y" and ending e
    if (is.null(lag)) {
        file_name <- sprintf(get_file_template_b(),
                             L, root_time, lambda, mu, beta, e)
    } else {
        file_name <- sprintf(get_file_template_a(),
                             L, root_time, lambda, mu, beta, lag, c, d, e)
    }
    return(file.path(file_dir, file_name))
}

make_file_name_ <- function(file_dir, grid, lag, c, d, e) {
    file_name <- make_file_name(file_dir, grid$L, grid$root_time, grid$lambda,
                                grid$mu, grid$beta, lag, c, d, e)
    return(file_name)
}

read_nexus_file <- function(target_dir, L, root_time, lambda, mu, beta) {
    file_name <- sprintf(
        file.path(target_dir, "data", "L%d_r%e_l%e_m%e_b%e.nex"),
        L, root_time, lambda, mu, beta
    )
    f <- ape::read.nexus(file_name)
    return(f)
}

read_nexus_file_ <- function(target_dir, grid) {
    f <- read_nexus_file(target_dir, grid$L, grid$root_time, grid$lambda,
                         grid$mu, grid$beta)
    return(f)
}

# Reading tree samples
get_trees <- function(tree_file) {
    trees <- ape::read.tree(tree_file, skip = 7, comment.char = "#")
    return(trees)
}

get_trees_ <- function(file_dir, grid, lag, c, d) {
    file_name <- make_file_name_(file_dir, grid, lag, c, d, "nex")
    trees <- get_trees(file_name)
    return(trees)
}

# Reading parameter samples
get_par_names <- function() {
    par_names <- c("sample", "log_prior", "integrated_llkd", "root_time", "mu",
                   "p", "lambda", "kappa", "rho", "ncat", "log_likelihood",
                   "beta")
    return(par_names)
}

get_pars <- function(par_file) {
    pars <- readr::read_table(par_file, get_par_names(), readr::cols(),
                              skip = 3)
    return(pars)
}

get_pars_ <- function(file_dir, grid, lag, c, d) {
    file_name <- make_file_name_(file_dir, grid, lag, c, d, "txt")
    pars <- get_pars(file_name)
    return(pars)
}

# Read config file and make settings grids
make_grid <- function(config_file) {
    source(config_file, TRUE)
    grid_config <- expand_grid(L = list_L, root_time = list_root_time,
                               lambda = list_lambda, mu = list_mu,
                               beta = list_beta)
    grid_a <- expand_grid(grid_config, run_length = list_run_length$coupled,
                          sample_interval = list_sample_interval$coupled,
                          lag = list_lag, c = list_c)
    grid_b <- tibble(grid_config, run_length = list_run_length$marginal,
                     sample_interval = list_sample_interval$marginal)
    return(list(grid_a = grid_a, grid_b = grid_b))
}

# Get coupling time
get_tau <- function(x, y, lag_offset) {
    # Assume x offset from y by by lag_offset indices
    tau <- get_tau_(x[-seq_len(lag_offset)] != y, lag_offset)
    return(tau)
}

get_tau_ <- function(z, lag_offset) {
    # z[i] <- x[i + lag_offset] != y[i]
    # last time chains differed since marginally they may couple and decouple
    t_off <- dplyr::last(which(z))
    if (is.na(t_off)) {
        t_off <- 0
    } else if (t_off == length(z)) {
        stop("chains did not couple")
    }
    tau <- t_off + lag_offset
    return(tau)
}

# Making coupling time figures
get_coupling_times <- function(out_dir, grid_a) {
    tau <- matrix(NA_integer_, nrow(grid_a))
    for (i in seq_along(tau)) {
        svMisc::progress(i, length(tau))
        grid_a_i <- grid_a[i, ]

        # sampled x and y
        x <- get_pars_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_x")
        y <- get_pars_(out_dir, grid_a_i, grid_a_i$lag, grid_a_i$c, "_y")
        # pre-computed tree distances
        z <- read.table(make_file_name_(out_dir, grid_a_i, grid_a_i$lag,
                                        grid_a_i$c, "", "dist"))

        lag_offset <- grid_a_i$lag / grid_a_i$sample_interval

        tau_p <- get_tau(x$log_prior, y$log_prior, lag_offset)
        tau_l <- get_tau(x$integrated_llkd, y$integrated_llkd, lag_offset)
        tau_r <- get_tau(x$root_time, y$root_time, lag_offset)
        tau_t <- get_tau_(z != 0, lag_offset)
        tau[i] <- max(c(tau_p, tau_l, tau_r, tau_t))
    }
    message("coupling times computed")
    return(tau)
}

make_w1_figure <- function(out_dir, fig_tau_data, grid_d, lag, iters) {

    w1_data <- expand_grid(fig_tau_data, iter = iters, w1_root = NA_real_,
                           w1_clade = NA_real_, w1_tree = NA_real_)  %>%
         group_by(L, root_time, lambda, mu, beta, run_length, sample_interval,
                  samp, tau)
    w1_keys <- group_keys(w1_data)
    w1_inds <- group_rows(w1_data)

    for (i in seq_len(nrow(w1_keys))) {
        svMisc::progress(i, nrow(w1_keys))
        grid_a_i <- w1_keys[i, 1:7]
        grid_c_i <- strtoi(w1_keys[i, 8])
        inds_i <- w1_inds[[i]]

        # Parameters
        p_x <- get_pars_(out_dir, grid_a_i, grid_c_i, "_x")$root_time
        p_y <- get_pars_(out_dir, grid_a_i, grid_c_i, "_y")$root_time
        tau_i <- w1_keys$tau[i]
        w1_data$w1_root[inds_i] <- w1_bound_estimators(tau_i, lag, iters, p_x,
                                                       p_y)

        # Clades and trees
        d_i <- plyr::match_df(grid_d, w1_keys[i, ],
                              c("L", "root_time", "lambda", "mu", "beta",
                                "run_length", "sample_interval"))
        t_x <- get_trees_(out_dir, grid_a_i, grid_c_i, "_x")
        t_y <- get_trees_(out_dir, grid_a_i, grid_c_i, "_y")

        c_x <- map_lgl(t_x, is.monophyletic, d_i$cl[[1]])
        c_y <- map_lgl(t_y, is.monophyletic, d_i$cl[[1]])
        w1_data$w1_clade[inds_i] <- w1_bound_estimators(tau_i, lag, iters, c_x,
                                                        c_y)
        g_x <- map_lgl(t_x, all.equal, d_i$tr[[1]], FALSE)
        g_y <- map_lgl(t_y, all.equal, d_i$tr[[1]], FALSE)
        w1_data$w1_tree[inds_i] <- w1_bound_estimators(tau_i, lag, iters, g_x,
                                                        g_y)
    }
    message("w1 bounds computed")

    fig_w1_data <- w1_data %>%
        ungroup(samp, tau) %>%
        select(-c(samp, tau)) %>%
        group_by(iter, .add = TRUE) %>%
        summarise(w1_root = mean(w1_root),
                  w1_clade = mean(w1_clade),
                  w1_tree = mean(w1_tree),
                  .groups = "drop") %>%
        pivot_longer(starts_with("w1_"), "stat", names_prefix = "w1_",
                     values_to = "w1") %>%
        mutate(across("stat", factor, c("root", "clade", "tree", "sum")))
    fig_w1 <- fig_w1_data %>%
        ggplot(aes(x = iter, y = w1, colour = as.factor(lambda))) +
        geom_line(size = 1.5, alpha = 0.75) +
        guides(colour = guide_legend(title = "lambda")) +
        labs(title = "W1 bound as leaf count L and birth rate lambda vary",
             subtitle = sprintf(
                 "%.02e iterations, coupling lag = %.02e, replications = %d",
                 w1_keys$run_length[1], w1_keys$sample_interval[1],
                 n_distinct(w1_keys$samp)),
             x = "iteration",
             y = "d_w1")
    for (scales in c("free", "fixed")) {
        fig_w1 +
            facet_wrap(~ L + stat, ncol = 3, scales = scales,
                       labeller = "label_both") +
            ggsave(sprintf(fig_template, sprintf("w1_axes-%s", scales)),
                   width = 10, height = 10)
    }
    fig_w1 +
        ylim(0, 2) +
        facet_wrap(~ L + stat, ncol = 3) +
        ggsave(sprintf(fig_template, "w1_axes-clipped"), width = 10,
               height = 10)
}

# Make marginal histograms
get_marginal_data <- function(out_dir, grid_a, grid_b, grid_c, par_name, k, m) {

    out_a <- tidyr::expand_grid(grid_a, c = grid_c)
    out_b <- tidyr::expand_grid(grid_b, c = 0)
    out_c <- bind_rows(out_a, out_b) %>%
        nest(x = -everything())
    for (i in seq_len(nrow(out_c))) {
        svMisc::progress(i, nrow(out_c))
        if (out_c$c[i] == 0) {
            x <- get_pars_(out_dir, out_c[i, 1:7], out_c$c[i], "")
            out_c$x[[i]] <- x[[par_name]][ind_x(k, grid_b$run_length[1] / grid_b$sample_interval[1])]

        } else {
            x <- get_pars_(out_dir, out_c[i, 1:7], out_c$c[i], "_x")
            out_c$x[[i]] <- x[[par_name]][ind_x(k, m)]
        }
    }
    message(sprintf("%s marginal data read", par_name))
    out <- tidyr::unnest(out_c, x)
    return(out)
}

make_marginal_hist <- function(out_dir, grid_a, grid_b, grid_c, par_name,
                               par_label, k = NULL, m = NULL) {
    if (is.null(m)) {
        m <- floor(grid_a$run_length[1] / grid_a$sample_interval[1])
    }
    if (is.null(k)) {
        k <- floor(m / 10)
    }
    out <- get_marginal_data(out_dir, grid_a, grid_b, grid_c, par_name, k, m)
    fig <- out %>%
        mutate(type = ifelse(c > 0, "a", "b")) %>%
        ggplot(aes(x = x, fill = as.factor(type), colour = NULL)) +
        stat_bin(aes(y = ..density..), bins = 30, position = "dodge") +
        labs(title = "Marginal distributions",
             subtitle = sprintf("x-chain samples %.02e to %.02e; (a) %d coupled chains and (b) 10x-longer chain", k, m, length(grid_c)),
             fill = NULL)
    for (scales in c("free", "fixed")) {
        fig +
        facet_wrap(~ L + lambda, ncol = 2, scales = scales,
                   labeller = "label_both") +
        ggsave(sprintf(fig_template,
                       sprintf("%s-hist_axes-%s", par_label, scales)),
               width = 8, height = 10, limitsize = FALSE)
    }
}

# Plot estimators
get_estimators <- function(out_dir, grid_a, grid_c, grid_d, par_name, k = NULL,
                           m = NULL) {
    if (is.null(m)) {
        m <- floor(grid_a$run_length[1] / grid_a$sample_interval[1])
    }
    if (is.null(k)) {
        k <- floor(m * c(1, 2, 5) / 10)
    }
    if (is.null(grid_d)) {
        out_l <- grid_a
    } else {
        out_l <- grid_d
    }
    out <- tidyr::expand_grid(out_l, c = grid_c, t = NA_integer_, k = k,
                              mc = NA_real_, bc = NA_real_, ue = NA_real_) %>%
        nest(k = k, s = c(mc, bc, ue))
    for (i in seq_len(nrow(out))) {
        svMisc::progress(i, nrow(out))

        xy <- get_samples_for_estimators(out_dir, out[i, ], par_name)
        x <- xy$x
        y <- xy$y

        out$t[i] <- get_tau(x, y)
        out$s[[i]] <- map_dfr(k, ~unbiased_estimator(x, y, ., m, out$t[i]))
    }
    message(sprintf("%s estimators computed", par_name))
    out <- unnest(out, cols = c(k, s))
    return(out)
}

get_samples_for_estimators <- function(out_dir, out_i, par_name) {
    if (par_name == "root_time") {
        x <- get_pars_(out_dir, out_i[1:7], out_i$c, "_x")[[par_name]]
        y <- get_pars_(out_dir, out_i[1:7], out_i$c, "_y")[[par_name]]
    } else {
        t_x <- get_trees_(out_dir, out_i[1:7], out_i$c, "_x")
        t_y <- get_trees_(out_dir, out_i[1:7], out_i$c, "_y")
        if (par_name == "topology support") {
            # TODO: check this
            v <- out_i$tr[[1]]
            x <- purrr::map_lgl(t_x, all.equal, v, FALSE)
            y <- purrr::map_lgl(t_y, all.equal, v, FALSE)
        } else if (par_name == "clade support") {
            # TODO: check this
            v <- out_i$cl[[1]]
            x <- purrr::map_lgl(t_x, ape::is.monophyletic, v)
            y <- purrr::map_lgl(t_y, ape::is.monophyletic, v)
        } else {
            stop("Incorrect par_name")
        }
    }
    return(list(x = x, y = y))
}

estimate_ground_truth <- function(out_dir, grid_b, grid_d, par_name) {
    m <- floor(grid_b$run_length[1] / grid_b$sample_interval[1])
    k <- floor(m / 2)

    out <- tibble(grid_b, mc = NA_real_)
    for (i in seq_len(nrow(out))) {
        svMisc::progress(i, nrow(out))
        if (is.null(grid_d)) {
            x <- get_pars_(out_dir, out[i, 1:7], 0, "")[[par_name]]
        } else {
            t <- get_trees_(out_dir, out[i, 1:7], 0, "")
            if (par_name == "topology support") {
                v <- grid_d$tr[[i]]
                x <- map_lgl(t, all.equal, v, FALSE)
            } else if (par_name == "clade support") {
                v <- grid_d$cl[[i]]
                x <- map_lgl(t, is.monophyletic, v)
            } else {
                stop("Incorrect par_name")
            }
        }
        out$mc[[i]] <- monte_carlo_estimator(x, k, m)
    }
    message(sprintf("%s ground truth computed", par_name))
    return(out)
}

make_estimator_figs <- function(out_dir, grid_a, grid_b, grid_c, grid_d,
                                par_name, par_label, k = NULL, m = NULL) {

    out_a <- get_estimators(out_dir, grid_a, grid_c, grid_d, par_name, k, m)
    out_b <- estimate_ground_truth(out_dir, grid_b, grid_d, par_name)

    make_estimator_hist(out_a, out_b, par_name, par_label)
    make_estimator_bias(out_a, par_name, par_label)
    make_estimator_mse(out_a, out_b, par_name, par_label)
}

make_estimator_hist <- function(out_a, out_b, par_name, par_label) {
    # TODO: add m to title
    # Coupled estimators
    fig_data_a <- out_a %>%
        pivot_longer(c(mc, bc, ue))
    # Average across chains
    fig_data_g <- fig_data_a %>%
        select(-t) %>%
        group_by(L, root_time, lambda, mu, beta, run_length, sample_interval, k,
                 name) %>%
        summarise(value = mean(value), .groups = "drop")
    # Ground truth estimates
    fig_data_b <- out_b %>%
        mutate(ue = mc) %>%
        pivot_longer(c(mc, ue))
    # Merge overall estimates
    fig_data_c <- bind_rows(bind_cols(fig_data_g, estimate = "average"),
                            bind_cols(fig_data_b, estimate = "ground truth"))

    fig <- fig_data_a %>%
        ggplot(aes(x = value)) +
        geom_density(aes(colour = as.factor(k)), alpha = 0.5, fill = NA) +
        geom_vline(aes(xintercept = value, colour = as.factor(k),
                       linetype = estimate),
                   fig_data_c, alpha = 0.8, size = 1) +
        guides(colour = guide_legend(title = "k"), fill = FALSE) +
        facet_wrap(~ L + lambda + name, ncol = 3, scales = "free",
                  labeller = "label_both") +
        labs(title = sprintf("Monte Carlo (mc) and bias-corrected (bc) unbiased (ue) estimators of %s",
                             par_name),
             subtitle = sprintf("Samples k to m\n%d pairs of coupled chains, length m = %.02e and subsample %.02e\n1 ground truth chain, length %.02e at subsample %.02e",
                                n_distinct(out_a$c),
                                out_a$run_length[1], out_a$sample_interval[1],
                                out_b$run_length[1], out_b$sample_interval[1]))
    for (scales in c("free", "fixed")) {
        fig +
        facet_wrap(~ L + lambda + name, ncol = 3, scales = scales,
                   labeller = "label_both") +
        ggsave(sprintf(fig_template,
                       sprintf("%s-est-hist_axes-%s", par_label, scales)),
               width = 8, height = 12, limitsize = FALSE)
    }
}

make_estimator_bias <- function(out_a, par_name, par_label) {
    # TODO: add m to title
    fig <- out_a %>%
        ggplot(aes(x = t, y = bc, fill = as.factor(k), colour = as.factor(k))) +
        geom_point(alpha = 0.6) +
        labs(title = sprintf("Coupling time tau and bias correction in unbiased estimators of %s",
                             par_name),
             subtitle = sprintf("Samples k to m; %d coupled chains of length %.02e at sampling interval %.02e",
                                n_distinct(out_a$c), out_a$run_length[1],
                                out_a$sample_interval[1]),
             x = "tau",
             y = "bias correction") +
        guides(colour = guide_legend(title = "k"), fill = FALSE)
    for (scales in c("free", "fixed")) {
        fig +
        facet_wrap(~ L + lambda, ncol = 2, scales = scales,
                   labeller = "label_both") +
        ggsave(sprintf(fig_template,
                       sprintf("%s-bc-tau_axes-%s", par_label, scales)),
               width = 8, height = 8, limitsize = FALSE)
    }
}

get_estimator_mse <- function(v, est) {
    m <- mean((v - est)^2)
    return(m)
}

make_estimator_mse <- function(out_a, out_b, par_name, par_label) {
    # TODO: add m to title
    if (par_name == "root_time") {
        out_f <- out_a %>%
            select(-c(c, t, bc))
    } else {
        out_f <- out_a %>%
            select(-c(cl, tr, c, t, bc))
    }
    fig_data <- out_f %>%
        nest(s = c(mc, ue)) %>%
        tibble(mse_mc_1 = NA_real_, mse_mc_n = NA_real_, mse_ue_n = NA_real_)

    for (i in seq_len(nrow(fig_data))) {
        mc_1 <- fig_data$s[[i]]$mc
        mc_n <- mean(mc_1)
        ue_n <- mean(fig_data$s[[i]]$ue)

        v <- out_b %>%
            semi_join(
                fig_data[i, ],
                c("L", "root_time", "lambda", "mu", "beta")
            ) %>%
            pull(mc)
        fig_data$mse_mc_1[i] <- get_estimator_mse(v, mc_1)
        fig_data$mse_mc_n[i] <- get_estimator_mse(v, mc_n)
        fig_data$mse_ue_n[i] <- get_estimator_mse(v, ue_n)
    }
    label_data <- out_a %>%
        select(L, root_time, lambda, mu, beta, t, k) %>%
        nest(t = t)
    label_data$mt <- map2_dbl(label_data$k, label_data$t, ~mean(.x >= .y$t))

    fig <- fig_data %>%
        pivot_longer(c(mse_mc_1, mse_mc_n, mse_ue_n), "estimator",
                     names_prefix = "mse_", values_to = "mse") %>%
        ggplot(aes(x = as.factor(k), y = mse, fill = estimator)) +
        geom_col(position = "dodge") +
        geom_label(aes(x = as.factor(k), y = 0, label = mt), label_data,
                   fill = NA, position = position_dodge(width = 1),
                   inherit.aes = FALSE) +
        labs(title = sprintf("MSE of Monte Carlo (mc) and unbiased (ue) estimators of %s",
                             par_name),
             subtitle = sprintf("Samples k to m in %d coupled chains; *_n is average of single chain/pair estimators *_1",
                                length(grid_c)),
             x = "k")
     for (scales in c("free", "fixed")) {
         fig +
         facet_wrap(~ L + lambda, ncol = 2, scales = scales,
                    labeller = "label_both") +
         ggsave(sprintf(fig_template,
                        sprintf("%s-mse_axes-%s", par_label, scales)),
                width = 8, height = 8)
     }
}

trace_estimator <- function(out_dir, grid_a, grid_b, grid_c, grid_d,
                                par_name, par_label) {

    m <- grid_a$run_length[1] / grid_a$sample_interval[1]
    k <- floor(m * c(1 / 10, 1 / 5, 2 / 5))

    if (is.null(grid_d)) {
        out_l <- grid_a
    } else {
        out_l <- grid_d
    }
    out_r <- tibble(k = k, m = map(k, ~seq.int(. + 100, m))) %>%
        unnest(m) %>%
        tibble(mc = NA_real_, bc = NA_real_, ue = NA_real_)
    fig_data <- tidyr::expand_grid(out_l, c = grid_c, out_r) %>%
        nest(km = k:m, s = mc:ue)

    # get estimators for each combination of k and m
    for (i in seq_len(nrow(fig_data))) {
        svMisc::progress(i, nrow(fig_data))

        xy <- get_samples_for_estimators(out_dir, fig_data[i, ], par_name)
        x <- xy$x
        y <- xy$y
        t <- get_tau(x, y)

        fig_data$s[[i]] <- map2_dfr(
            fig_data$km[[i]]$k,
            fig_data$km[[i]]$m,
            ~unbiased_estimator(x, y, .x, .y, t)
        )
    }
    message("estimator terms computed")

    # get mse for each estimator and sd of bias correction
    fig_data <- fig_data %>%
        unnest(c(km, s)) %>%
        nest(s = c(c, mc:ue)) %>%
        tibble(mse_mc_1 = NA_real_, mse_mc_n = NA_real_, mse_ue_n = NA_real_)
    out_b <- estimate_ground_truth(out_dir, grid_b, grid_d, par_name)

    for (i in seq_len(nrow(fig_data))) {
        svMisc::progress(i, nrow(fig_data))

        mc_1 <- fig_data$s[[i]]$mc
        mc_n <- mean(mc_1)
        ue_n <- mean(fig_data$s[[i]]$ue)

        v <- out_b %>%
            semi_join(
                fig_data[i, ],
                c("L", "root_time", "lambda", "mu", "beta")
            ) %>%
            pull(mc)
        fig_data$mse_mc_1[i] <- get_estimator_mse(v, mc_1)
        fig_data$mse_mc_n[i] <- get_estimator_mse(v, mc_n)
        fig_data$mse_ue_n[i] <- get_estimator_mse(v, ue_n)
    }
    message("mse terms computed")

    fig1 <- fig_data %>%
        pivot_longer(c(mse_mc_1, mse_mc_n, mse_ue_n), "estimator",
                     names_prefix = "mse_", values_to = "mse") %>%
        ggplot(aes(x = m, y = mse, colour = as.factor(k),
                   linetype = estimator)) +
        geom_line(alpha = 0.75) +
        labs(title = sprintf("MSE of Monte Carlo (mc) and unbiased (ue) estimators of %s as k and m vary",
                             par_name),
             subtitle = sprintf("Samples k to m in n = %d coupled chains",
                                length(grid_c)),
             x = "m",
             colour = "k") +
         scale_y_continuous(trans = "log1p") +
         facet_wrap(~ L + lambda, ncol = 2, scales = "free",
                    labeller = "label_both") +
         ggsave(sprintf(fig_template, sprintf("%s-mse-trace", par_label)),
                width = 8, height = 8)

    # plot sd of bc
    fig_data$sd_bc <- map_dbl(fig_data$s, ~sd(.$bc))
    fig <- fig_data %>%
        ggplot(aes(x = m, y = sd_bc, colour = as.factor(k))) +
        geom_line(alpha = 0.75) +
        labs(title = sprintf("SD of bias correction (bc) in unbiased estimators of %s as k and m vary",
                             par_name),
             subtitle = sprintf("Samples k to m in n = %d coupled chains",
                                length(grid_c)),
             x = "m",
             y = "sd(bc)",
             colour = "k") +
         scale_y_continuous(trans = "log1p") +
         facet_wrap(~ L + lambda, ncol = 2, scales = "free",
                    labeller = "label_both") +
         ggsave(sprintf(fig_template, sprintf("%s-bc-trace", par_label)),
                width = 8, height = 8)
}
