library("ggplot2")

# Plotting bias correction weights and differences
tibble(x = x[-1], i = gl(5, k = floor(length(x) / 5))) %>%
    ggplot(aes(x = x, fill = i)) +
    geom_density(colour = NA, alpha = 0.5) +
    facet_wrap(~ i, ncol = 1)

plot(1:1e4, cummean(x[ind_x(1, 1e4)]), xlim = c(1e3, 1e4))

par(mfrow = c(2, 1))
plot(w, d, "l",
    main = TeX("$ bc = \\Sigma_{i = k+1}^{\\tau - 1} w_i d_i $"),
    xlab = TeX("$ w_i = 1 \\vee \\frac{i - k}{m - k + 1} $"),
    ylab = TeX("$ d_i = x_{i + 1} - y_i $"))
plot(w_inds, cumsum(w * d), "l", main = "cumulative sum", xlab = TeX("$ k + 1 = 101, ..., i, ..., \\tau - 1 = 463 $"))

# Plotting Monte Carlo and unbiased estimators
# I was only plotting one experiment when I wrote this code
out_a %>%
    select(-c(L, root_time, lambda, mu, beta, run_length, sample_interval, c,
              bc)) %>%
    ggplot(aes(x = mc, y = ue, colour = t - 1 > k + 1)) +
    geom_point(alpha = 0.75) +
    facet_wrap(~ k, ncol = 1)

fig_data <- fig_data %>%
    select(-c(L, root_time, lambda, mu, beta, run_length, sample_interval))

a0 <- out_a %>% select(t, k, mc, bc, ue)
a1 <- a0 %>% filter(k == 100)
a2 <- a0 %>% filter(k == 200)
a5 <- a0 %>% filter(k == 500)
all(a5$t - 1 < a5$k + 1)

i1 <- a1$k + 1 < a1$t - 1 # bc1 != 0
all((abs(a2$bc) < abs(a1$bc))[i1]) # |bc| decreases with k
all(!(abs(a2$bc) < abs(a1$bc))[!i1]) # except when its already 0

# so why does mse of ue increase when k goes from 100 to 200
