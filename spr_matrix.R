#!/usr/bin/Rscript

library("ape")
library("readr")

source("/home/users/kelly/rspr/R/rspR.r")

target <- commandArgs(TRUE)

t <- read.nexus(sprintf("%s.supp", target))
d <- rspr.matrix(t, "full", Inf)

write_csv(data.frame(d), sprintf("%s.supp_dists", target), col_names = FALSE)
