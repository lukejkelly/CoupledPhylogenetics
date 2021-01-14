#!/usr/bin/Rscript

library("ape")
library("purrr")
library("readr")

target <- commandArgs(TRUE)

x <- read.tree(sprintf("%s_%s.nex", target, "x"), skip = 7, comment.char = "#")
y <- read.tree(sprintf("%s_%s.nex", target, "y"), skip = 7, comment.char = "#")

d <- map2_dbl(x[-1], y, dist.topo, "score")

write_lines(d, sprintf("%s.dist", target))
