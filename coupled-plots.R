# Trace plot outputs from coupled TraitLabSDLT output
library("ape")
library("tidyverse")
# library("zoo")
# library("latex2exp")

source("~/Workspace/RSPR/R/rspR.r")
source("traitlab-output-tools.R")

fig_stem <- sprintf("figs/%s-L5", lubridate::ymd(Sys.Date()))
dir.create(fig_stem)
fig_template <- sprintf("%s/%%s.pdf", fig_stem)

output_template <- "../output/SIM_B5/SIM_B5-%d"

i_c <- 1:10
n_c <- length(i_c)
output_stems <- map_chr(seq_len(2 * n_c), ~sprintf(output_template, .))

# trees is a list of length n_runs, each containing x and y output trees
trees <- map(output_stems, get_trees)

# Burn-in and total number of samples
n_t <- length(trees[[1]]$x)
n_b <- 0 #floor(n_t / 10)
i_s <- seq(n_b + 1, n_t)
n_s <- length(i_s)

d_tree <- tree_distances(trees, i_s, i_c)

# Getting parameter samples for same chains
pars <- map(output_stems, get_pars)
i_p <- c(3, 4, 5, 12, 13)
n_p <- length(i_p)
d_par <- par_distances(pars, i_s, i_c, i_p)

# Putting all the outputs together
d <- rbind(d_tree, d_par)

# Trace plots
trace_all <- d %>%
    ggplot(aes(x = samp, y = value, colour = var, fill = NULL)) +
    geom_point(alpha = 0.1, size = 0.1) +
    stat_summary(fun = mean, geom = "line", show.legend = TRUE, alpha = 0.8) +
    labs(x = "subsample index", y = "distance") +
    facet_grid(rows = vars(stat), scales = "free_y", switch = "y") +
    scale_y_continuous(trans = "log1p") +
    ggtitle(sprintf("Average distances across %d pairs of coupled chains",
                    n_c)) +
    ggsave(sprintf(fig_template, "trace-all"), width = 10, height = 10)

hist_all<- d %>%
    ggplot(aes(x = value, colour = NULL, fill = var)) +
    geom_histogram(position = "dodge") +
    facet_grid(~ stat, scales = "free", switch = "y") +
    scale_x_continuous(trans = "log1p") +
    ggtitle(sprintf(paste("Distances across %d pairs of coupled chains",
                           "(length %d, burn-in %d)"),
                    n_c, n_t, n_b)) +
    ggsave(sprintf(fig_template, "hist-all"), width = 15, height = 5)

# Same for individual pairs of chains
pdf(sprintf(fig_template, "trace-ind"), width = 10, height = 10)
walk(i_c, ~print(trace_ind(.)))
dev.off()

pdf(sprintf(fig_template, "hist-ind"), width = 15, height = 5)
walk(i_c, ~print(hist_ind(.)))
dev.off()

# Marginal plots
x_pars <- array(dim = c(n_s, n_c, n_p))
y_pars <- array(dim = c(n_s, n_c, n_p))
for (j in seq_along(i_c)) {
    p_xj <- pars[[i_c[j]]]$x
    p_yj <- pars[[i_c[j]]]$y

    for (k in seq_along(i_p)) {
        x_pars[, j, k] <- p_xj[i_s, i_p[k]]
        y_pars[, j, k] <- p_yj[i_s, i_p[k]]
    }
}

par_names <- names(pars[[1]]$x)[i_p]
m <- rbind(
    map_dfr(seq_along(i_p),
            ~data.frame(var = "x", stat = par_names[.], samp = i_s,
                        x_pars[, , .])),
    map_dfr(seq_along(i_p),
            ~data.frame(var = "y", stat = par_names[.], samp = i_s,
                        y_pars[, , .]))) %>%
    pivot_longer(c(-var, -stat, -samp), "chain", names_prefix = "X") %>%
    mutate(chain = factor(chain, levels = i_c))

trace_marg <- m %>%
    ggplot(aes(x = samp, y = value, colour = chain, fill = NULL)) +
    geom_line(size = 0.25, alpha = 0.9) +
    facet_grid(rows = vars(stat), cols = vars(var), scales = "free",
               switch = "y") +
    # scale_y_continuous(trans = "log1p") +
    ggtitle("Trace plot of coupled chains") +
    ggsave(sprintf(fig_template, "trace-marg"), width = 30, height = 15)

hist_marg <- m %>%
    ggplot(aes(x = value, colour = chain, linetype = var)) +
    # geom_histogram(position = "dodge", alpha = 1) +
    stat_ecdf(alpha = 0.66) +
    facet_grid(cols = vars(stat), scales = "free", switch = "y") +
    ggtitle(sprintf(paste("Marginal distributions of %d coupled chains",
                          "(length %d, burn-in %d)"),
                    n_c, i_s[n_s], i_s[1] - 1)) +
    ggsave(sprintf(fig_template, "hist-marg"), width = 10, height = 5)


###### What's going on...
d1 <- filter(d, stat %in% c("beta", "mu"), chain %in% c(1,3))

ggplot(d1, aes(x = samp, y = value, col = var)) +
geom_line(alpha = 0.8) +
facet_grid(rows = vars(stat), cols = vars(chain))
