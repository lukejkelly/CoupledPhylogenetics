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
