library("ape")
library("tidyverse")
library("fs")

source("coupling-times-functions.R")

targets <- c("20210113", "20210128", "20210129", "20210130", "20210131",
             "20210132", "20210133")

# Make grids of config and run settings
grid_bs <- map(targets, ~make_grid(file.path("..", ., "config.R"))$grid_b)

# Plotting root times, although they all seem pretty unimodal
for (i in seq_along(targets)) {

    target <- targets[i]
    grid_b <- grid_bs[[i]]

    # Target figure and output templates and directories
    fig_dir <- file.path("..", target, "figs")
    dir_create(fig_dir)
    fig_template <- sprintf("%s/%%s.pdf", fig_dir)

    # Parameters for histograms
    k <- (grid_b$run_length[1] / grid_b$sample_interval[1]) / 2
    m <- grid_b$run_length[1] / grid_b$sample_interval[1]
    inds <- seq.int(k, m)

    out <- grid_b %>%
        expand_grid(inds = inds, root = NA_real_) %>%
        group_by(L, root_time, lambda, mu, beta, run_length, sample_interval)
    out_keys <- group_keys(out)
    out_rows <- group_rows(out)

    out_dir <- file.path("..", target, "output")
    for (j in seq_along(out_rows)) {
        print(j / length(out_rows), digits = 3)
        x <- get_pars_(out_dir, out_keys[j, 1:7], 0, "")$root_time[inds + 1]
        out$root[out_rows[[j]]] <- x
    }

    fig <- out %>%
        ggplot(aes(x = root, fill = as.factor(lambda), colour = NULL)) +
        stat_bin(aes(y = ..density..), bins = 50, position = "dodge") +
        labs(title = "Root time distributions from long marginal chains",
             subtitle = sprintf("Samples %.02e to %.02e, subsample %d",
                                k, m, grid_b$sample_interval[1]),
             fill = NULL) +
        facet_wrap(~ L, ncol = 1, scales = "free", labeller = "label_both") +
        ggsave(sprintf(fig_template, sprintf("root-%s", target)),
               width = 8, height = 10, limitsize = FALSE)

}

################################################################################
# Plotting number of trees with support above threshold
labels <- c("normal", "cat below root", "mix of trees 1", "mix of trees 2",
            "lots of cats", "half data with borrowing", "wrong clade in mcmc")
names(labels) <- targets
tree_out <- map2_dfr(grid_bs, targets,
                     ~expand_grid(.x, config = .y, count = NA_integer_))
for (i in seq_len(nrow(tree_out))) {
    print(i / nrow(tree_out), digits = 3)
    target_file <- make_file_name_(file.path("..", tree_out$config[i],
                                   "output"), tree_out[i, 1:7], 0, "", "supp")
    x <- 0
    try(x <- ape::read.nexus(target_file))
    tree_out$count[i] <- length(x)
}

fig <- tree_out %>%
    ggplot(aes(x = count, colour = NULL)) +
    stat_count() +
    labs(title = "Number of trees with support >10%",
         subtitle = sprintf("Samples %.02e to %.02e, subsample %d",
                            k, m, grid_b$sample_interval[1]),
         fill = NULL) +
    facet_wrap(~ config, ncol = 1, scales = "fixed",
               labeller = as_labeller(labels)) +
    ggsave("../trees-10.pdf", width = 6, height = 10, limitsize = FALSE)

################################################################################
# Plotting posterior distributions on topologies
supp_out <- map2_dfr(grid_bs, targets,
                     ~expand_grid(.x, config = .y, supp = list(NA_real_)))
for (i in seq_len(nrow(supp_out))) {
    print(i / nrow(supp_out), digits = 3)
    target_file <- make_file_name_(file.path("..", supp_out$config[i],
                                   "output"), supp_out[i, 1:7], 0, "", "supps")
    supp_out$supp[[i]] <- scan(target_file)
}

fig_supp_data <- supp_out %>%
    unnest(supp) %>%
    group_by(L, root_time, lambda, mu, beta, run_length, sample_interval,
             config) %>%
    mutate(index = row_number()) %>%
    ungroup()

fig_supp <- fig_supp_data %>%
    ggplot(aes(x = config, y = supp, colour = NULL,
           fill = as.factor(index))) +
    geom_bar(position = "stack", stat = "identity", show.legend = FALSE) +
    labs(title = "Topology posterior support", fill = NULL) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1,
                                     size = 10)) +
    scale_x_discrete(labels = labels) +
    scale_fill_manual(values = rep_len(c("#E69F00", "#56B4E9", "#009E73",
                                       "#F0E442"), 1377)) +
    facet_wrap(~ L + lambda, ncol = 4, labeller = "label_both") +
    ggsave("../figs/topology-support.pdf", width = 10, height = 10)

################################################################################
# Plotting posterior distributions and SPR distances
target <- "20210113"
grid_b <- make_grid(file.path("..", target, "config.R"))$grid_b
supp_spr_out <- expand_grid(grid_b, supp = list(NA_real_),
                            spr = list(NA_integer_))

# read supports and spr distance matrices
for (i in seq_len(nrow(supp_spr_out))) {
    print(i / nrow(supp_spr_out), digits = 3)
    supp_file <- make_file_name_(file.path("..", target, "output"),
                                   supp_spr_out[i, 1:7], 0, "", "supps")
    supp_spr_out$supp[[i]] <- scan(supp_file)
    spr_file <- sub("supps", "supp_dists", supp_file)
    supp_spr_out$spr[[i]] <- as.matrix(read_csv(spr_file, FALSE))
}

for (i in seq_len(nrow(supp_spr_out))) {
    print(i)
    supp <- supp_spr_out$supp[[i]]
    spr <- supp_spr_out$spr[[i]]

    n <- length(supp)
    # if (n == 1) {
    #     next
    # }
    rownames(spr) <- colnames(spr) <- seq_len(n)

    g <- graph_from_adjacency_matrix(spr, "undirected", "dist")
    g <- set_vertex_attr(g, "supp", value = supp)

    h <- as_tbl_graph(g) %>%
        activate(edges) %>%
        mutate(weight = 1 / dist) %>%
        activate(edges) %>%
        mutate(spr = as.factor(dist))

    ggraph(h, "kk") +
    geom_edge_link(aes(edge_alpha = weight, edge_colour = spr),
                   lineend = "round",
                   show.legend = c(edge_alpha = FALSE, edge_colour = TRUE)) +
    geom_node_circle(aes(r = sqrt(supp / (2 * pi)), fill = supp)) +
    labs(title = sprintf("L = %d, lambda = %g", supp_spr_out$L[i],
                         supp_spr_out$lambda[i])) +
    ggsave(sprintf("../20210113/figs/g%02i-L%02d-l%e.pdf", i, supp_spr_out$L[i],
                   supp_spr_out$lambda[i]), width = 5, height = 5)

    ggraph(filter(h, dist == 1), "kk") +
    geom_edge_link(lineend = "round") +
    geom_node_circle(aes(r = sqrt(supp / (2 * pi)), fill = supp)) +
    labs(title = sprintf("L = %d, lambda = %g", supp_spr_out$L[i],
                         supp_spr_out$lambda[i])) +
    ggsave(sprintf("../20210113/figs/h%02i-L%02d-l%e.pdf", i, supp_spr_out$L[i],
                   supp_spr_out$lambda[i]), width = 5, height = 5)
}

# pdfjam --nup 4x4 -o edges-all.pdf --papersize '{15cm,15cm}' g*.pdf
# pdfjam --nup 4x4 -o edges-spr1.pdf --papersize '{15cm,15cm}' h*.pdf
