library("tidyr")
library("dplyr")
library("ape")
library("ggplot2")
library("purrr")

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
    trees <- ape::read.tree(tree_file, skip = 7, comment.char = "#", keep.multi = TRUE)
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
    pars <- readr::read_table2(par_file, get_par_names()[c(2:5, 8:12)],
                              "-dddd--ddidd-", skip = 3)
    return(pars)
}

get_pars_ <- function(file_dir, grid, lag, c, d) {
    file_name <- make_file_name_(file_dir, grid, lag, c, d, "txt")
    pars <- get_pars(file_name)
    return(pars)
}

get_pa_xy <- function(file_dir, grid, lag, c, moves) {
    file_name <- make_file_name_(file_dir, grid, lag, c, "_xy", "pa")
    pa_xy <- readr::read_csv(file_name, col_types = readr::cols()) %>%
        filter(t <= grid$tau)
    t <- pa_xy$t
    pa <- pa_xy %>%
        select(-t) %>%
        select(!!moves) %>%
        # mutate(across(everything(), ~replace(., is.nan(.), 0))) %>%
        rename_with(~ paste0("m_", .))
    return(tibble(t, pa))
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
get_coupling_times_ <- function(out_dir, grid_a) {
    # estimate coupling times from output
    tau <- rep(NA_integer_, nrow(grid_a))
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
        tau_b <- get_tau(x$beta, y$beta, lag_offset)
        tau_c <- get_tau(x$ncat, y$ncat, lag_offset)
        tau[i] <- max(c(tau_p, tau_l, tau_r, tau_t, tau_b, tau_c))
    }
    message("coupling times computed")
    return(tau)
}

get_coupling_times <- function(out_dir, grid_a) {
    # get coupling times written by TraitLab
    tau <- rep(NA_integer_, nrow(grid_a))
    for (i in seq_along(tau)) {
        svMisc::progress(i, length(tau))
        grid_a_i <- grid_a[i, ]
        file_name <- make_file_name_(out_dir, grid_a_i, grid_a_i$lag,
                                     grid_a_i$c, "_x", "tau")
        tau[i] <- scan(file_name)
    }
    message("coupling times computed")
    return(tau)
}

make_tau_ecdf <- function(grid_a) {
    fig_tau <- grid_a %>%
        ggplot(aes(
            x = tau - lag / sample_interval,
            colour = factor(sprintf("%.4g", lag),
                            levels = sprintf("%.4g", unique(sort(lag))))
        )) +
        stat_ecdf(pad = FALSE, alpha = 0.75) +
        labs(title = "ECDF of coupling time tau",
             subtitle = sprintf("minimum %.02e iterations, replications = %d",
                                grid_a$run_length[1], n_distinct(grid_a$c)),
             x = sprintf("(tau - lag) / %d", grid_a$sample_interval[1]),
             y = "ECDF", # latex2exp::TeX(r'($\hat{F}(\tau)$)'),
             colour = "lag")
    if (n_distinct(grid_a$L) > 1) {
        for (scales in c("free", "fixed")) {
            fig_p <- fig_tau +
                facet_wrap(~ L, ncol = n_distinct(grid_a$L), scales = scales,
                           labeller = "label_both")
            ggsave(sprintf(fig_template, sprintf("tau-ecdf_axes-%s", scales)),
                   fig_p,
                   width = 3 * n_distinct(grid_a$L) + 2,
                   height = 3)
        }
    } else {
        ggsave(sprintf(fig_template, "tau"), fig_tau, width = 5, height = 3)
    }
}

get_estimator_k <- function(grid_a) {
    grid_a <- grid_a %>%
        nest(gr = c(c, tau)) %>%
        mutate(k = map(gr, ~ quantile(.$tau, 0.99))) %>%
        unnest(cols = c(gr, k))
    return(grid_a)
}

make_tv_figure <- function(out_dir, grid_a, iters) {
    tv_data <- expand_grid(grid_a, iter = iters, tv = NA_real_)
    tv_data$tv <- tv_bound_estimator(tv_data$tau,
                                     tv_data$lag / tv_data$sample_interval,
                                     tv_data$iter)
    fig_tv_data <- tv_data %>%
        select(-tau) %>%
        nest(s = c(c, tv)) %>%
        mutate(tv = map_dbl(s, ~mean(.$tv)))
    fig_tv <- fig_tv_data %>%
        ggplot(aes(
            x = iter,
            y = tv,
            colour = factor(sprintf("%.4g", lag),
                            levels = sprintf("%.4g", unique(sort(lag))))
        )) +
        geom_line(size = 1, alpha = 0.75) +
        labs(title = "TV upper bound",
             subtitle = sprintf("minimum %.02e iterations, replications = %d",
                                grid_a$run_length[1], n_distinct(grid_a$c)),
             x = sprintf("iteration / %d", grid_a$sample_interval[1]),
             y = "TV upper bound",
             colour = "lag")
     if (n_distinct(grid_a$L) > 1) {
         for (scales in c("free", "fixed")) {
             fig_p <- fig_tv +
                 # scale_y_continuous(trans = "log1p") +
                 facet_wrap(~ L, ncol = n_L, scales = scales,
                            labeller = "label_both")
             ggsave(sprintf(fig_template, sprintf("tv_axes-%s", scales)),
                    fig_p,
                    width = 3 * n_distinct(grid_a$L) + 2,
                    height = 3)
         }
     } else {
         ggsave(sprintf(fig_template, "tv"), fig_tv, width = 5, height = 3)
     }
}

make_w1_figure <- function(out_dir, grid_a, grid_d, iters) {
    w1_data <- expand_grid(grid_a, iter = iters, w1_root = NA_real_,
                           w1_clade = NA_real_, w1_tree = NA_real_)  %>%
        nest(s = c(iter, w1_root, w1_clade, w1_tree))

    for (i in seq_len(nrow(w1_data))) {
        svMisc::progress(i, nrow(w1_data))

        w1_i <- w1_data[i, ]
        lag_i <- w1_i$lag
        c_i <- w1_i$c
        tau_i <- w1_i$tau
        si_i <-  w1_i$sample_interval

        # root time
        p_x <- get_pars_(out_dir, w1_i, lag_i, c_i, "_x")$root_time
        p_y <- get_pars_(out_dir, w1_i, lag_i, c_i, "_y")$root_time

        w1_data$s[[i]]$w1_root <- map_dbl(
            w1_data$s[[i]]$iter,
            ~w1_bound_estimator(p_x, p_y, ., tau_i, lag_i /  si_i)
        )

        # Clades and trees
        d_i <- plyr::match_df(grid_d, w1_i,
                              c("L", "root_time", "lambda", "mu", "beta",
                                "run_length", "sample_interval"))
        t_x <- get_trees_(out_dir, w1_i, lag_i, c_i, "_x")
        t_y <- get_trees_(out_dir, w1_i, lag_i, c_i, "_y")

        c_x <- map_lgl(t_x, is.monophyletic, d_i$cl[[1]])
        c_y <- map_lgl(t_y, is.monophyletic, d_i$cl[[1]])

        w1_data$s[[i]]$w1_clade <- map_dbl(
            w1_data$s[[i]]$iter,
            ~w1_bound_estimator(c_x, c_y, ., tau_i, lag_i / si_i)
        )

        g_x <- map_lgl(t_x, all.equal, d_i$tr[[1]], FALSE)
        g_y <- map_lgl(t_y, all.equal, d_i$tr[[1]], FALSE)

        w1_data$s[[i]]$w1_tree <- map_dbl(
            w1_data$s[[i]]$iter,
            ~w1_bound_estimator(g_x, g_y, ., tau_i, lag_i / si_i)
        )
    }
    message("w1 bounds computed")

    fig_w1_data <- w1_data %>%
        unnest(s) %>%
        select(-c(tau)) %>%
        nest(s = c(c, w1_root, w1_clade, w1_tree)) %>%
        mutate(
            w1_root = map_dbl(s, ~mean(.$w1_root)),
            w1_clade = map_dbl(s, ~mean(.$w1_clade)),
            w1_tree = map_dbl(s, ~mean(.$w1_tree))
        ) %>%
        pivot_longer(starts_with("w1_"), "stat", names_prefix = "w1_",
                     values_to = "w1") %>%
        mutate(across("stat", factor, c("root", "clade", "tree")))
    fig_w1 <- fig_w1_data %>%
        ggplot(aes(
            x = iter,
            y = w1,
            colour = factor(sprintf("%.4g", lag),
                                    levels = sprintf("%.4g", unique(sort(lag))))
        )) +
        geom_line(alpha = 0.75) +
        labs(title = "Estimated W1 bound",
             subtitle = sprintf("%.02e iterations, replications = %d",
                                w1_data$run_length[1], n_distinct(w1_data$c)),
             x = sprintf("iteration / %d", si_a),
             y = "d_w1",
             colour = "lag") +
        facet_wrap(~ L + lambda + stat, ncol = 3, scales = "free",
                   labeller = "label_both")
    ggsave(sprintf(fig_template, "w1_axes-free"),
           fig_w1,
           width = 3 * 3 + 2,
           height = prod(3, n_distinct(grid_a$L), n_distinct(grid_a$lambda)))
}

# Make marginal histograms
get_marginal_data <- function(out_dir, grid_a, grid_b, par_name, k, m) {
    if (all(grid_b$run_length == 0)) {
        out_r <- grid_a
    } else {
        out_r <- bind_rows(grid_a, tibble(grid_b(grid_b$run_length > 0),
                                          lag = NA, c = NA))
    }
    out <- out_r %>%
        nest(x = -everything())
    for (i in seq_len(nrow(out))) {
        svMisc::progress(i, nrow(out))
        out_i <- out[i, ]
        if (is.na(out_i$c)) {
            x <- get_pars_(out_dir, out_i, NULL)[[par_name]]
            # use second half of long chain samples
            m0 <- out_i$run_length / out_i$sample_interval
            k0 <- floor(m0 / 2)
            out$x[[i]] <- x[ind0(k0, m0)]
        } else {
            x <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_x")[[par_name]]
            out$x[[i]] <- x[ind0(k[i], m[i])]
        }
    }
    message(sprintf("%s marginal data read", par_name))
    out <- tidyr::unnest(out, x)
    return(out)
}

make_marginal_hist <- function(out_dir, grid_a, grid_b, par_name, par_label,
                               k = NULL, m = NULL) {
    if (is.null(m)) {
        if (all(grid_a$run_length > 0)) {
            m <- floor(grid_a$run_length / grid_a$sample_interval)
        } else {
            m <- floor(grid_a$lag / grid_a$sample_interval)
        }
    }
    if (is.null(k)) {
        k <- floor(m / 5)
    }
    out <- get_marginal_data(out_dir, grid_a, grid_b, par_name, k, m)
    fig_data <- out %>%
        mutate(type = ifelse(is.na(c), "ground truth", paste("lag", lag)))
    if (all(out$x %% 1 == 0)) {
        fig <- fig_data %>%
            ggplot(aes(x = x, colour = as.factor(type), fill = as.factor(type)),
                       alpha = 0.5) +
            geom_histogram(aes(), binwidth = 1, center = 0, position = "dodge")
    } else {
        fig <- fig_data %>%
            ggplot(aes(x = x, colour = as.factor(type)), fill = NULL,
                   alpha = 0.5) +
            geom_density(aes())
    }
    fig <- fig +
        labs(title = sprintf("Marginal distribution of %s", par_name),
             subtitle = sprintf("x-chain samples %.02e to %.02e in %d coupled chains",
                                k, m, n_distinct(grid_a$c)),
             colour = NULL,
             fill = NULL,
             x = par_label)
    if (n_distinct(grid_a$L) > 1 || n_distinct(grid_a$lambda) > 1) {
        for (scales in c("free", "fixed")) {
            fig_p <- fig +
            facet_wrap(~ L + lambda, ncol = 2, scales = scales,
                       labeller = "label_both")
            ggsave(sprintf(fig_template, sprintf("%s-hist_axes-%s", par_label,
                                                 scales)),
                   fig_p,
                   width = 3 * n_distinct(grid_a$lambda) + 2,
                   height = 3 * n_distinct(grid_a$L))
        }
    } else {
        ggsave(sprintf(fig_template, sprintf("%s-hist", par_label)),
               fig, width = 5, height = 3)
    }
}

# Plot estimators
make_estimator_figs <- function(out_dir, grid_a, grid_b, grid_d, par_name,
                                par_label, k = NULL, m = NULL) {

    out_a <- get_estimators(out_dir, grid_a, grid_d, par_name, k, m)
    out_b <- estimate_ground_truth(out_dir, grid_b, grid_d, par_name)

    make_estimator_hist(out_a, out_b, par_name, par_label)
    make_estimator_bias(out_a, par_name, par_label)
    make_estimator_mse(out_a, out_b, par_name, par_label)
}

get_estimators <- function(out_dir, grid_a, grid_d, par_name, k = NULL,
                           m = NULL) {
    if (is.null(m)) {
        m <- floor(grid_a$run_length[1] / grid_a$sample_interval[1])
    }
    if (is.null(k)) {
        k <- floor(m * c(1, 2) / 10)
    }
    if (is.null(grid_d)) {
        out_l <- grid_a
    } else {
        out_l <- left_join(grid_a,
                           select(grid_d, -c(run_length, sample_interval)),
                           c("L", "root_time", "lambda", "mu", "beta"))
    }
    out <- tidyr::expand_grid(out_l, k = k, mc = NA_real_, bc = NA_real_,
                              ue = NA_real_) %>%
        nest(k = k, s = c(mc, bc, ue))
    for (i in seq_len(nrow(out))) {
        svMisc::progress(i, nrow(out))
        out_i <- out[i, ]
        xy <- get_samples_for_estimators(out_dir, out_i, par_name)
        x <- xy$x
        y <- xy$y
        out$s[[i]] <- map_dfr(
            k,
            ~unbiased_estimator(x, y, ., m, out_i$tau,
                                out_i$lag / out_i$sample_interval)
        )
    }
    message(sprintf("%s estimators computed", par_name))
    out <- unnest(out, cols = c(k, s))
    return(out)
}

get_samples_for_estimators <- function(out_dir, out_i, par_name) {
    if (par_name == "root_time") {
        x <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_x")[[par_name]]
        y <- get_pars_(out_dir, out_i, out_i$lag, out_i$c, "_y")[[par_name]]
    } else {
        t_x <- get_trees_(out_dir, out_i, out_i$lag, out_i$c, "_x")
        t_y <- get_trees_(out_dir, out_i, out_i$lag, out_i$c, "_y")
        if (par_name == "topology support") {
            v <- out_i$tr[[1]]
            x <- purrr::map_lgl(t_x, all.equal, v, FALSE)
            y <- purrr::map_lgl(t_y, all.equal, v, FALSE)
        } else if (par_name == "clade support") {
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

    if (is.null(grid_d)) {
        out_l <- grid_b
    } else {
        out_l <- right_join(grid_b, grid_d,
                            c("L", "root_time", "lambda", "mu", "beta",
                              "run_length", "sample_interval"))
    }
    out <- tibble(out_l, mc = NA_real_)
    for (i in seq_len(nrow(out))) {
        svMisc::progress(i, nrow(out))
        out_i <- out[i, ]
        if (par_name == "root_time") {
            x <- get_pars_(out_dir, out_i, NULL)[[par_name]]
        } else {
            t <- get_trees_(out_dir, out_i, NULL)
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

make_estimator_hist <- function(out_a, out_b, par_name, par_label) {
    # Coupled estimators
    fig_data_a <- out_a %>%
        pivot_longer(c(mc, bc, ue))
    # Average across chains
    fig_data_g <- fig_data_a %>%
        select(-tau) %>%
        group_by(L, root_time, lambda, mu, beta, run_length, sample_interval,
                 lag, k, name) %>%
        summarise(value = mean(value), .groups = "drop")
    # Ground truth estimates
    fig_data_b <- out_b %>%
        tidyr::expand_grid(lag = unique(out_a$lag), k = unique(out_a$k)) %>%
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
             subtitle = sprintf("%d pairs of coupled chains, samples k to m = %.02e (interval %.02e)\n1 ground truth chain, length %.02e (interval %.02e)",
                                n_distinct(out_a$c),
                                out_a$run_length[1], out_a$sample_interval[1],
                                out_b$run_length[1], out_b$sample_interval[1])) +
             guides(colour = guide_legend(title = sprintf("k / %d", out_a$sample_interval)))
    for (scales in c("free", "fixed")) {
        fig_p <- fig +
            facet_wrap(~ L + lambda + lag + name, ncol = 3, scales = scales,
                       labeller = "label_both")
        ggsave(sprintf(fig_template,
                       sprintf("%s-est-hist_axes-%s", par_label, scales)),
               fig_p,
               width = 3 * 3 + 2,
               height = 3 * prod(n_distinct(out_a$L), n_distinct(out_a$lambda),
                                 n_distinct(out_a$lag)))
    }
}

make_estimator_bias <- function(out_a, par_name, par_label) {
    # TODO: add m to title
    fig <- out_a %>%
        ggplot(aes(x = tau, y = bc, fill = as.factor(k),
               colour = as.factor(k))) +
        geom_point(alpha = 0.6) +
        labs(title = sprintf("Coupling time tau and bias correction in unbiased estimators of %s",
                             par_name),
             subtitle = sprintf("Samples k to m = %.02e; %d pairs of coupled chains",
                                out_a$run_length[1], n_distinct(out_a$c)),
             x = sprintf("tau / %d", out_a$sample_interval[1]),
             y = "bias correction") +
        guides(colour = guide_legend(title = sprintf("k / %d",
                                                     out_a$sample_interval[1])),
               fill = FALSE)
    for (scales in c("free", "fixed")) {
        fig_p <- fig +
            facet_wrap(~ L + lambda + lag, ncol = n_distinct(out_a$lag),
                       scales = scales, labeller = "label_both")
        ggsave(sprintf(fig_template,
                       sprintf("%s-bc-tau_axes-%s", par_label, scales)),
               fig_p,
               width = 3 * n_distinct(out_a$lag) + 2,
               height = 3 * n_distinct(out_a$L) * n_distinct(out_a$lambda))
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
            select(-c(c, tau, bc))
    } else {
        out_f <- out_a %>%
            select(-c(c, tau, cl, tr, bc))
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
        select(L, root_time, lambda, mu, beta, lag, tau, k) %>%
        nest(tau = tau)
    label_data$mt <- map2_dbl(label_data$k, label_data$tau, ~mean(.x >= .y$tau))

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
             subtitle = sprintf("Using samples k to m/tau in %d coupled chains (*_n is average of single chain/pair estimators *_1)",
                                length(n_distinct(out_a$c))),
             x = "k")
     for (scales in c("free", "fixed")) {
         fig_p <- fig +
             facet_wrap(~ L + lambda + lag, ncol = n_distinct(out_a$lag),
                        scales = scales, labeller = "label_both")
         ggsave(sprintf(fig_template,
                        sprintf("%s-mse_axes-%s", par_label, scales)),
                fig_p,
                width = 3 * n_distinct(out_a$lag) + 2,
                height = 3 * n_distinct(out_a$L) * n_distinct(out_a$lambda))
     }
}

trace_estimator <- function(out_dir, grid_a, grid_b, grid_d, par_name,
                            par_label) {

    m <- grid_a$run_length[1] / grid_a$sample_interval[1]
    k <- floor(m * c(1, 2, 4) / 10)

    if (is.null(grid_d)) {
        out_l <- grid_a
    } else {
        out_l <- left_join(grid_a, grid_d,
                           c("L", "root_time", "lambda", "mu", "beta",
                             "run_length", "sample_interval"))
    }
    out_r <- tibble(k = k, m = map(k, ~seq.int(. + 100, m, 10))) %>%
        unnest(m) %>%
        tibble(mc = NA_real_, bc = NA_real_, ue = NA_real_)
    fig_data <- tidyr::expand_grid(out_l, out_r) %>%
        nest(km = k:m, s = mc:ue)

    # get estimators for each combination of k and m
    for (i in seq_len(nrow(fig_data))) {
        svMisc::progress(i, nrow(fig_data))

        xy <- get_samples_for_estimators(out_dir, fig_data[i, ], par_name)
        x <- xy$x
        y <- xy$y

        fig_data$s[[i]] <- map2_dfr(
            fig_data$km[[i]]$k,
            fig_data$km[[i]]$m,
            ~unbiased_estimator(x, y, .x, .y, fig_data$tau[i],
                                fig_data$lag[i] / fig_data$sample_interval[i])
        )
    }
    message("estimator terms computed")

    # get mse for each estimator and sd of bias correction
    fig_data <- fig_data %>%
        unnest(c(km, s)) %>%
        select(-tau) %>%
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
        ggplot(aes(x = m, y = mse, colour = as.factor(k))) +
        geom_line(alpha = 0.75) +
        labs(title = sprintf("MSE of Monte Carlo (mc) and unbiased (ue) estimators of %s as k and m vary",
                             par_name),
             subtitle = sprintf("Samples k to m in n = %d coupled chains at interval %.02e",
                                n_distinct(grid_a$c),
                                grid_a$sample_interval[1]),
             x = "m",
             colour = "k") +
         scale_y_continuous(trans = "log1p") +
         facet_wrap(~ L + lambda + lag + estimator, ncol = 3,
                    scales = "free", labeller = "label_both")
     ggsave(sprintf(fig_template, sprintf("%s-mse-trace", par_label)),
            fig1,
            width = 3 * 3 + 2,
            height = 3 * prod(n_distinct(grid_a$L), n_distinct(grid_a$lambda),
                              n_distinct(grid_a$lag)))


    # plot sd of bc
    fig_data$sd_bc <- map_dbl(fig_data$s, ~sd(.$bc))
    fig2 <- fig_data %>%
        ggplot(aes(x = m, y = sd_bc, colour = as.factor(k))) +
        geom_line(alpha = 0.75) +
        labs(title = sprintf("SD of bias correction (bc) in unbiased estimators of %s as k and m vary",
                             par_name),
             subtitle = sprintf("Samples k to m in n = %d coupled chains at interval %.02e",
                                n_distinct(grid_a$c),
                                grid_a$sample_interval[1]),
             x = "m",
             y = "sd(bc)",
             colour = "k") +
         scale_y_continuous(trans = "log1p") +
         facet_wrap(~ L + lambda + lag, ncol = n_distinct(grid_a$lag),
                    scales = "free", labeller = "label_both")
     ggsave(sprintf(fig_template, sprintf("%s-bc-trace", par_label)),
            fig2,
            width = 3 * 3 + 2,
            height = 3 * n_distinct(grid_a$L) * n_distinct(grid_a$lambda))
}

# TraitLab move names
move_names <- function(inds) {
    mns <- c("slide", "exchange (local)", "exchange (wide)", "spr (local)",
             "spr (wide)", "scale (tree)", "scale (subtree)", "mu", NA, NA,
             "leaf time", "scale (top)", "cat (add)", "cat (delete)", "rho",
             "kappa", NA, "cat (move)", "xi (one)", "xi (all)", "beta")
    mn <- mns[inds]
    return(mn)
}
