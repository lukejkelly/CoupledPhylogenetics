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
    labs(title = "Topology posterior support",
         # subtitle = sprintf("Samples %.02e to %.02e, subsample %d",
         #                    k, m, grid_b$sample_interval[1]),
         fill = NULL) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1,
                                     size = 10)) +
    scale_x_discrete(labels = labels) +
    scale_fill_manual(values = rep_len(c("#E69F00", "#56B4E9", "#009E73",
                                       "#F0E442"), 166)) +
    facet_wrap(~ L + lambda, ncol = 4, labeller = "label_both") +
    ggsave("../figs/topology-support.pdf", width = 10, height = 10)


fig_supps_data <- fig_supp_data %>%
    filter(L == 10, lambda == 0.025, config == "20210113")
fig_supp <- fig_supps_data  %>%
    ggplot(aes(x = config, y = supp, fill = as.factor(index))) +
    geom_col() +
    labs(title = "Topology posterior support",
         # subtitle = sprintf("Samples %.02e to %.02e, subsample %d",
         #                    k, m, grid_b$sample_interval[1]),
         fill = NULL) +
    # theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1,
    #                                  size = 10)) +
    # scale_x_discrete(labels = labels) +
    # facet_wrap(~ L + lambda, ncol = 4, labeller = "label_both") +
    ggsave("../figs/topology-support.pdf", width = 10, height = 10)
