output_pattern <- sprintf("L\\d?_r%1$s_l%1$s_m%1$s_b%1$s\\-\\d+_x\\.txt", "\\de[\\+|\\-]\\d{2}")
output_files <- dir_ls(path = output_dir, regexp = output_pattern)

##

L <- as.list(seq(3, 15, 3))
lambda <- rep(list(seq(0.1, 0.5, 0.2)), length(L))
i_c <- rep(list(seq_len(25)), times = length(L))
i_c <- rep(list(seq_len(100)), times = 5)

output <- map_dfr(
    1:length(L),
    ~data.frame(L = rep(L[[.]], each = length(lambda[[.]]) * length(i_c[[.]])),
                lambda = rep(lambda[[.]], each = length(i_c[[.]])),
                i = rep(i_c[[.]], times = length(lambda[[.]])),
                tau = NA))
output_files <- output %>%
    pmap_chr(~sprintf(output_template, ..1, ..2, ..3, ..5, ..6), "x", "txt")
output$tau <- output_files %>%
    map_int(~ifelse(file_exists(.), nrow(read.table(., FALSE, skip = 3)), NA))

output_means <- output %>%
    group_by(L, lambda) %>%
    summarise(tau_bar = mean(tau), .groups = "drop")
