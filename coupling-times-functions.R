# Functions to construct file names and paths
get_file_template <- function() {
    # file_template <- "L%d_r%.0e_l%.0e_m%.0e_b%.0e_n%.0e_s%.0e-%d%s.%s"
    file_template <- "L%d_r%e_l%e_m%e_b%e_n%e_s%e-%d%s.%s"
    return(file_template)
}

make_file_name <- function(file_dir, L, root_time, lambda, mu, beta, run_length,
                           sample_interval, c, d, e) {
    # Run index is c, decorator d and ending e
    file_name <- file.path(
        file_dir,
        sprintf(get_file_template(), L, root_time, lambda,
                mu, beta, run_length, sample_interval, c, d, e))
    return(file_name)
}

make_file_name_ <- function(file_dir, grid, c, d, e) {
    file_name <- make_file_name(file_dir, grid$L, grid$root_time, grid$lambda,
                                grid$mu, grid$beta, grid$run_length,
                                grid$sample_interval, c, d, e)
    return(file_name)
}

# Reading tree samples
get_trees <- function(tree_file) {
    trees <- ape::read.tree(tree_file, skip = 8, comment.char = "#")
    return(trees)
}

get_trees_ <- function(file_dir, grid, c, d) {
    file_name <- make_file_name_(file_dir, grid, c, d, "nex")
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

get_pars_ <- function(file_dir, grid, c, d) {
    file_name <- make_file_name_(file_dir, grid, c, d, "txt")
    pars <- get_pars(file_name)
    return(pars)
}

# Read config file and make settings grids
make_grid <- function(config_file) {
    source(config_file, TRUE)
    grid_config <- expand_grid(L = list_L, root_time = list_root_time,
                               lambda = list_lambda, mu = list_mu,
                               beta = list_beta)
    grid_a <- expand_grid(grid_config, run_length = grid_run_length[1],
                          sample_interval = grid_sample_interval[1])
    grid_b <- expand_grid(grid_config, run_length = grid_run_length[2],
                          sample_interval = grid_sample_interval[2])
    return(list(grid_a = grid_a, grid_b = grid_b, grid_c = grid_c))
}

# Get coupling time
get_tau <- function(x, y) {
    # Assume x offset by 1 index
    tau <- get_tau_(x[-1] != y)
    return(tau)
}

get_tau_ <- function(z) {
    tau <- dplyr::last(which(z)) + 1
    if (is.na(tau)) {
        tau <- 1
    }
    return(tau)
}

# Making coupling time figures
get_coupling_times <- function(out_dir, grid_a, grid_c) {
    tau <- matrix(NA_integer_, nrow(grid_a), length(grid_c))
    for (i in seq_len(nrow(tau))) {
        print(i / nrow(tau), digits = 3)
        grid_a_i <- grid_a[i, ]
        for (j in seq_along(grid_c)) {
            grid_c_j <- grid_c[j]
            # sampled x and y
            x <- get_pars_(out_dir, grid_a_i, grid_c_j, "_x")
            y <- get_pars_(out_dir, grid_a_i, grid_c_j, "_y")
            # pre-computed tree distances
            z <- read.table(make_file_name_(out_dir, grid_a_i, grid_c_j, "",
                            "dist"))

            tau_p <- get_tau(x$log_prior, y$log_prior)
            tau_l <- get_tau(x$integrated_llkd, y$integrated_llkd)
            tau_r <- get_tau(x$root_time, y$root_time)
            tau_t <- get_tau_(z != 0)

            tau[i, j] <- max(c(tau_p, tau_l, tau_r, tau_t))
        }
    }
    out <- cbind(grid_a, tau = tau) %>%
        pivot_longer(starts_with("tau"), "samp", names_prefix = "tau.",
                     values_to = "tau")
    return(out)
}

# Make marginal histograms
get_marginal_data <- function(out_dir, grid_a, grid_b, grid_c, k, m, par_name) {
    inds <- seq.int(k, m)

    out_a <- tidyr::expand_grid(grid_a, c = grid_c)
    out_b <- tidyr::expand_grid(grid_b, c = 0)
    out_c <- rbind(out_a, out_b) %>%
        nest(x = -everything())
    for (i in seq_len(nrow(out_c))) {
        print(i / nrow(out_c), digits = 3)
        x <- get_pars_(out_dir, out_c[i, 1:7], out_c$c[i], "_x")
        out_c$x[[i]] <- x[[par_name]][inds]
    }
    out <- tidyr::unnest(out_c, x)
    return(out)
}

make_marginal_hist <- function(out_dir, grid_a, grid_b, grid_c, k, m,
                                 par_name, par_label) {
    out <- get_marginal_data(out_dir, grid_a, grid_b, grid_c, k, m, par_name)
    fig <- out %>%
        mutate(type = ifelse(c > 0, "a", "b")) %>%
        ggplot(aes(x = x, fill = as.factor(type), colour = NULL)) +
        stat_bin(aes(y = ..density..), bins = 50, position = "dodge") +
        labs(title = "Marginal distributions",
             subtitle = sprintf("x-chain samples %.02e to %.02e; (a) %d coupled chains and (b) 10x-longer chain", k, m, length(grid_c)),
             fill = NULL)
    for (scales in c("free", "fixed")) {
        fig +
        facet_wrap(~ L + lambda, ncol = 3, scales = scales,
                   labeller = "label_both") +
        ggsave(sprintf(fig_template,
                       sprintf("%s-hist_axes-%s", par_label, scales)),
               width = 8, height = 10, limitsize = FALSE)
    }
}

# Plot estimators
get_estimators <- function(out_dir, grid_a, grid_c, k, m, par_name,
                           grid_clade) {
    out <- tidyr::expand_grid(grid_a, c = grid_c, t = NA_integer_,
                              mc = NA_real_, bc = NA_real_, ue = NA_real_)
    for (i in seq_len(nrow(out))) {
        print(i / nrow(out), digits = 3)
        if (is.null(grid_clade)) {
            x <- get_pars_(out_dir, out[i, 1:7], out$c[i], "_x")[[par_name]]
            y <- get_pars_(out_dir, out[i, 1:7], out$c[i], "_y")[[par_name]]
        } else {
            v <- purrr::pluck(grid_clade$cl, which(out$L[i] == grid_clade$L))
            x <- get_trees_(out_dir, out[i, 1:7], out$c[i], "_x") %>%
                map_lgl(is.monophyletic, v)
            y <- get_trees_(out_dir, out[i, 1:7], out$c[i], "_y") %>%
                map_lgl(is.monophyletic, v)
        }

        out$t[i] <- get_tau(x, y)
        f <- unbiased_estimator(x, y, k, m, out$t[i])
        out$mc[i] <- f["mc"]
        out$bc[i] <- f["bc"]
        out$ue[i] <- f["ue"]
    }
    return(out)
}

estimate_ground_truth <- function(out_dir, grid_b, k, m, par_name, grid_clade) {
    out <- tidyr::expand_grid(grid_b, mc = NA_real_)
    for (i in seq_len(nrow(out))) {
        print(i / nrow(out), digits = 3)
        if (is.null(grid_clade)) {
            x <- get_pars_(out_dir, out[i, 1:7], 0, "_x")[[par_name]]
        } else {
            v <- purrr::pluck(grid_clade$cl, which(out$L[i] == grid_clade$L))
            x <- get_trees_(out_dir, out[i, 1:7], 0, "_x") %>%
                map_lgl(is.monophyletic, v)
        }
        out$mc[i] <- monte_carlo_estimator(x, k, m)
    }
    return(out)
}

make_estimator_hist <- function(out_dir, grid_a, grid_b, grid_c, k, m,
                                par_name, par_label, grid_clade = NULL) {
    out_a <- get_estimators(out_dir, grid_a, grid_c, k, m, par_name, grid_clade)
    out_b <- estimate_ground_truth(out_dir, grid_b, k, m, par_name, grid_clade)

    # Histograms
    fig_data_a <- out_a %>%
        pivot_longer(c(mc, bc, ue))
    fig_data_g <- fig_data_a %>%
        select(-t) %>%
        group_by(L, root_time, lambda, mu, beta, run_length, sample_interval,
                 name) %>%
        summarise(x_mean = mean(value), .groups = "drop")
    fig_data_b <- out_b %>%
        mutate(ue = mc) %>%
        pivot_longer(c(mc, ue))

    fig <- fig_data_a %>%
        ggplot(aes(x = value)) +
        stat_bin(aes(y = ..density.., group = as.factor(lambda)),
                 bins = 20, position = "dodge", alpha = 0.5) +
        geom_vline(mapping = aes(xintercept = x_mean, colour = "blue",
                                 linetype = 3),
                   alpha = 0.8, data = fig_data_g, size = 1, linetype = 2,
                   show.legend = FALSE) +
        geom_vline(mapping = aes(xintercept = value), colour = "red",
                  data = fig_data_b, size = 1.5, linetype = 3) +
        guides(fill = guide_legend(title = "lambda"), colour = FALSE) +
        labs(title = sprintf("Monte Carlo (mc) and bias-corrected (bc) unbiased (ue) estimators of %s",
                             par_name),
             subtitle = sprintf("Samples %.02e to %.02e; histogram of %d coupled estimators, dotted line is their average;\ndashed line is MC estimator from 10x-longer chain",
                                k, m, length(grid_c)))
    for (scales in c("free", "fixed")) {
        fig +
        facet_wrap(~ L + lambda + name, ncol = 3, scales = scales,
                   labeller = "label_both") +
        ggsave(sprintf(fig_template,
                       sprintf("%s-est-hist_axes-%s", par_label, scales)),
               width = 8, height = 30, limitsize = FALSE)
    }
}

make_estimator_bias <- function(out_dir, grid_a, grid_c, k, m, par_name,
                               par_label, grid_clade = NULL) {
    out <- get_estimators(out_dir, grid_a, grid_c, k, m, par_name, grid_clade)

    fig <- out %>%
        ggplot(aes(x = t, y = bc)) +
        geom_point(alpha = 0.2) +
        # geom_density2d_filled(show.legend = c(fill = TRUE)) +
        labs(title = sprintf("Coupling time tau and bias correction in unbiased estimators of %s",
                             par_name),
             subtitle = sprintf("Samples %.02e to %.02e; %d coupled estimators",
                                k, m, length(grid_c)),
             x = "tau",
             y = "bias correction") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    for (scales in c("free", "fixed")) {
        fig +
        facet_wrap(~ L + lambda, ncol = 3, scales = scales,
                   labeller = "label_both") +
        ggsave(sprintf(fig_template,
                       sprintf("%s-bc-tau_axes-%s", par_label, scales)),
               width = 8, height = 10, limitsize = FALSE)
    }
}

get_estimator_mse <- function(v, est) {
    m <- mean((v - est)^2)
    return(m)
}

make_estimator_mse <- function(out_dir, grid_a, grid_b, grid_c, k, m,
                               par_name, par_label, grid_clade = NULL) {
    out_a <- get_estimators(out_dir, grid_a, grid_c, k, m, par_name, grid_clade)
    out_b <- estimate_ground_truth(out_dir, grid_b, k, m, par_name, grid_clade)

    # MSEs
    fig_data <- expand_grid(grid_b, mc_1 = NA_real_, mc_n = NA_real_,
                            ue_n = NA_real_)
    for (i in seq_len(nrow(fig_data))) {
        out_i <- filter(out_a, L == fig_data$L[i], lambda == fig_data$lambda[i])
        mc_1 <- out_i$mc
        mc_n <- mean(mc_1)
        ue_n <- mean(out_i$ue)

        v <- out_b$mc[i]
        fig_data$mc_1[i] <- get_estimator_mse(v, mc_1)
        fig_data$mc_n[i] <- get_estimator_mse(v, mc_n)
        fig_data$ue_n[i] <- get_estimator_mse(v, ue_n)
    }

    fig <- fig_data %>%
        pivot_longer(c(mc_1, mc_n, ue_n), "estimator", values_to = "mse") %>%
        ggplot(aes(x = estimator, y = mse)) +
        geom_col() +
        labs(title = sprintf("MSE of Monte Carlo (mc) and unbiased (ue) estimators of %s",
                             par_name),
             subtitle = sprintf("%d coupled chains, samples %g to %g; *_n is average of *_1 which uses a single chain (pair)",
                         length(grid_c), k, m))
     for (scales in c("free", "fixed")) {
         fig +
         facet_wrap(~ L + lambda, ncol = 3, scales = scales,
                    labeller = "label_both") +
         ggsave(sprintf(fig_template,
                        sprintf("%s-mse_axes-%s", par_label, scales)),
                width = 8, height = 10)
     }
}

make_w1_figure <- function(out_dir, fig_tau_data, grid_clade, lag, iters) {

    w1_data <- expand_grid(fig_tau_data, iter = iters, w1_root = NA_real_,
                           w1_clade = NA_real_)  %>%
         group_by(L, root_time, lambda, mu, beta, run_length, sample_interval,
                  samp, tau)
    w1_keys <- group_keys(w1_data)
    w1_inds <- group_rows(w1_data)

    for (i in seq_len(nrow(w1_keys))) {
        print(i / nrow(w1_keys), digits = 3)
        grid_a_i <- w1_keys[i, 1:7]
        grid_c_i <- strtoi(w1_keys[i, 8])
        inds_i <- w1_inds[[i]]

        # Parameters
        p_x <- get_pars_(out_dir, grid_a_i, grid_c_i, "_x")$root_time
        p_y <- get_pars_(out_dir, grid_a_i, grid_c_i, "_y")$root_time
        tau_i <- w1_keys$tau[i]
        w1_data$w1_root[inds_i] <- w1_bound_estimators(tau_i, lag, iters, p_x,
                                                       p_y)

        # Clades
        v <- purrr::pluck(grid_clade$cl, which(w1_keys$L[i] == grid_clade$L))
        t_x <- get_trees_(out_dir, grid_a_i, grid_c_i, "_x") %>%
            map_lgl(is.monophyletic, v)
        t_y <- get_trees_(out_dir, grid_a_i, grid_c_i, "_y") %>%
            map_lgl(is.monophyletic, v)
        w1_data$w1_clade[inds_i] <- w1_bound_estimators(tau_i, lag, iters, t_x,
                                                        t_y)
    }

    fig_w1_data <- w1_data %>%
        ungroup(samp, tau) %>%
        select(-c(samp, tau)) %>%
        group_by(iter, .add = TRUE) %>%
        summarise(w1_root = mean(w1_root), w1_clade = mean(w1_clade),
                  w1_sum = sum(w1_root + w1_clade), .groups = "drop") %>%
        pivot_longer(starts_with("w1_"), "stat", names_prefix = "w1_",
                     values_to = "w1")
    fig_w1 <- fig_w1_data %>%
        ggplot(aes(x = iter, y = w1, colour = as.factor(lambda))) +
        geom_line(size = 1.5, alpha = 0.75) +
        guides(colour = guide_legend(title = "lambda")) +
        labs(title = "W1 upper bound as number of leaves L and birth rate lambda vary",
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
