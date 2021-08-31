#  Analysing output

Steps to analyse output from coupled MCMC experiments. The primary script for  analysing output is `coupling-plots.R`, ignore any lines which are commented out. Many of the others scripts are for checking things and are not needed to plot the output of an experiment

This code uses a modified version of the `RWTY` `R` package to compute ASDSF on disjoint sliding windows samples, steps to install it are below. Alternatively, `RWTY` can be installed from `CRAN` and will produce cumulative ASDSF estimates.

## Description
Start `R` in the directory of a completed experiment; for example, `CoupledPhylogenetics/20210831` and run code in `coupling-plots.R`

* Load some `R` packages, which may need to be installed on your system
* Source various functions (`coupling-functions.R`, ...) to deal with output and create figures
* Read the `config.R` file from an experiment and create a grid where each row corresponds to an experiment
    * `grid_a` is the coupled experiments
    * `grid_b` is the corresponding marginal experiments
* Read the coupling times from files, plot the ECDF and total variation bounds
    * Plotting functions may need to be edited depending on the range of parameters for experiments: currently they `facet_wrap` over `L`, include `mu` if you have more than one value in your experiments, etc.
* Plot histograms or kernel density estimates using (longer) X components of coupled chains
    * This depends on the length of the chain so is difficult to compare if chains are a different random length, can be edited to use samples `k, ..., m` for burn-in `k` and minimum chain length `m`
* Create sliding window ASDSF plot
    * Takes the data from different plot objects and replots it together

---

## Install modified `RWTY`

Download the `disjoint` branch of (https://github.com/lukejkelly/RWTY). The version at time of writing was 1.0.2, so I built it with the following script.
```bash
R CMD build --no-build-vignettes RWTY &&
R CMD check --no-build-vignettes --no-manual rwty_1.0.2.tar.gz &&
R CMD INSTALL --library=~/Workspace/Rlibs/ rwty_1.0.2.tar.gz
```
My `~/Renviron` file contains
```bash
R_LIBS=/Users/kelly/Workspace/Rlibs
```
The modified version of `RWTY` is installed in `~/Workspace/Rlibs` and loaded in `coupling-plots.R` via
```R
library("rwty", lib.loc = "~/Workspace/Rlibs/")
```

---

##  Unit tests
To run all unit tests, execute
```R
source("estimators.R")
source("coupling-functions.R")
testthat::test_dir("tests/testthat")
```
or
```R
source("<file>.R")
testthat::test_file("tests/testthat/test-<file>.R")
```
for individual tests.

Add `reporter = c("summary", "debug")` as an argument to `test_*` to pause execution at failed `expect_*` instances; for example,
```R
testthat::test_dir("tests/testthat", reporter = c("summary", "debug"))
```
