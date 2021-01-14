# Tools for analysing marginal/coupled TraitLab output
# Each function assumes its relevant files are named <file_stem>_{x,y}.{nex,txt}

library("ape")

# Getting samples for a pair of chains
get_trees <- function(file_stem) {
    file_x <- sprintf("%s_%s.nex", file_stem, "x")
    file_y <- sprintf("%s_%s.nex", file_stem, "y")
    trees_x <- read.tree(file_x, skip = 8, comment.char = "#")
    trees_y <- read.tree(file_y, skip = 7, comment.char = "#")
    return(list(x = trees_x, y = trees_y))
}

get_pars <- function(file_stem) {
    pars_x <- read.table(sprintf("%s_%s.txt", file_stem, "x"), FALSE, skip = 4)
    pars_y <- read.table(sprintf("%s_%s.txt", file_stem, "y"), FALSE, skip = 3)
    pars_x$V13 <- pars_x$V12 / pars_x$V5
    pars_y$V13 <- pars_y$V12 / pars_y$V5
    names(pars_x) <- c("sample", "log prior", "integrated llkd", "root time",
                       "mu", "p", "lambda", "kappa", "rho", "ncat",
                       "log likelihood", "beta", "beta/mu")
    names(pars_y) <- names(pars_x)
    return(list(x = pars_x, y = pars_y))
}

# Distances between chain samples
tree_distances <- function(trees, i_s, i_c) {
    n_s <- length(i_s)
    n_c <- length(i_c)

    d_xy_topo <- matrix(NA, n_s, n_c)
    d_xy_rspr <- matrix(NA, n_s, n_c)
    d_xz_topo <- matrix(NA, n_s, n_c)
    d_xz_rspr <- matrix(NA, n_s, n_c)

    for (j in seq_along(i_c)) {
        t_xj <- trees[[i_c[j]]]$x
        t_yj <- trees[[i_c[j]]]$y
        t_zj <- trees[[i_c[j] + 10]]$x

        # d_xy_topo[, j] <- map_dbl(i_s,
        #                           ~dist.topo(t_xj[[.]], t_yj[[.]], "score"))
        # d_xz_topo[, j] <- map_dbl(i_s,
        #                           ~dist.topo(t_xj[[.]], t_zj[[.]], "score"))

        d_xy_rspr[, j] <- map_dbl(i_s,
                                  ~rspr(t_xj[[.]], t_yj[[.]]))
        d_xz_rspr[, j] <- map_dbl(i_s,
                                  ~rspr(t_xj[[.]], t_zj[[.]]))
    }

    d <- rbind(
        # data.frame(var = "(x,y)", stat = "topo", samp = i_s, d_xy_topo),
        # data.frame(var = "(x,z)", stat = "topo", samp = i_s, d_xz_topo),
        data.frame(var = "(x,y)", stat = "rspr", samp = i_s, d_xy_rspr),
        data.frame(var = "(x,z)", stat = "rspr", samp = i_s, d_xz_rspr)) %>%
        pivot_longer(c(-var, -stat, -samp), "chain", names_prefix = "X") %>%
        mutate(chain = as.factor(chain))
    return(d)
}

par_distances <- function(pars, i_s, i_c, i_p) {
    n_s <- length(i_s)
    n_c <- length(i_c)
    n_p <- length(i_p)

    d_xy_pars <- array(dim = c(n_s, n_c, n_p))
    d_xz_pars <- array(dim = c(n_s, n_c, n_p))

    for (j in seq_along(i_c)) {
        p_xj <- pars[[i_c[j]]]$x
        p_yj <- pars[[i_c[j]]]$y
        p_zj <- pars[[i_c[j] + 10]]$x

        for (k in seq_along(i_p)) {
            d_xy_pars[, j, k] <- abs(p_xj[i_s, i_p[k]] - p_yj[i_s, i_p[k]])
            d_xz_pars[, j, k] <- abs(p_xj[i_s, i_p[k]] - p_zj[i_s, i_p[k]])
        }
    }

    par_names <- names(pars[[1]]$x)[i_p]
    d <- rbind(
        map_dfr(seq_along(i_p),
                ~data.frame(var = "(x,y)", stat = par_names[.], samp = i_s,
                            d_xy_pars[, , .])),
        map_dfr(seq_along(i_p),
                ~data.frame(var = "(x,z)", stat = par_names[.], samp = i_s,
                            d_xz_pars[, , .]))) %>%
        pivot_longer(c(-var, -stat, -samp), "chain", names_prefix = "X") %>%
        mutate(chain = as.factor(chain))
    return(d)
}

# Plotting individual pairs of chains
trace_ind <- function(chain) {
    fig <- d %>%
        filter(chain == !!chain) %>%
        ggplot(aes(x = samp, y = value, colour = var, fill = NULL)) +
        geom_line(size = 0.25) +
        geom_point(aes(y = ifelse(value == 0, 0, NA), fill = "dodgerblue1",
                   colour = NULL, alpha = 0.33)) +
        labs(x = "subsample index", y = "distance") +
        facet_grid(rows = vars(stat), cols = vars(var), scales = "free",
                   switch = "y") +
        scale_y_continuous(trans = "log1p") +
        ggtitle(sprintf("Distances across pair of coupled chains %d", chain)) +
        theme(legend.position = "none")
    return(fig)
}

hist_ind <- function(chain) {
    fig <- d %>%
        filter(chain == !!chain) %>%
        ggplot(aes(x = value, colour = NULL, fill = var)) +
        geom_histogram(position = "dodge") +
        facet_grid(~ stat, scales = "free", switch = "y") +
        scale_x_continuous(trans = "log1p") +
        ggtitle(sprintf(paste0("Distances across pair of coupled chains %d ",
                               "(length %d, burn-in %d)"),
                        chain, n_t, n_b))
    return(fig)
}
