#  Analysing output

Steps to analyse output from coupled MCMC experiments. The primary script for  analysing output is `coupling-plots.R`, ignore any lines which are commented out. Many of the others scripts are for checking things and are not needed to plot the output of an experiment

This code uses a modified version of the `RWTY` `R` package to compute ASDSF on disjoint sliding windows samples, steps to install it are below. (Alternatively, `RWTY` can be installed from `CRAN` and will produce cumulative ASDSF estimates instead --- see instructions below.)

## Description
Start `R` in the directory of a completed experiment --- for example, `CoupledPhylogenetics/20210831` --- and run the code in `coupling-plots.R` which will

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

## Install modified `rwty`

Instead of computing the ASDSF on disjoint or cumulative sliding windows, we take a similar approach to MrBayes and compute it on the most recent `window.lookback` proportion of samples at `window.number` evenly spaced sample indices after discarding `burnin`. In order to do so, we use a modified version of `rwty`.

If `devtools` is available, then execute the following within `R` to install the package.
```R
devtools::install_github("https://github.com/lukejkelly/RWTY", "lagged-asdsf")
```
Otherwise, the code may be installed manually.

This version of `rwty` includes one additional function, `makeplot.asdsf.mb()`, to compute ASDSF on sliding windows which grow in size.
```R
makeplot.asdsf.mb(
    chains,  # chains object, same as makeplot.asdsf
    burnin = 0,  # burnin to discard, same as makeplot.asdsf
    window.lookback = 0.75,  # percentage of most recent samples to use for each window
    window.number = 10,  # number of evenly spaced window endpoints
    min.freq = 0.0,  # only consider split frequencies above this threshold, same as makeplot.asdsf
    log.y = TRUE  # plot ASDSF on log axis, same as makeplot.asdsf
)
```
For example, if each of our chains have 111 samples `X_0, ..., X_110`, or `X[1], ..., X[111]` in R's indexing, and we set

 - `burnin = 11`
 - `window.lookback = 0.75`
 - `window.number = 4`

then we have 100 samples after discard burnin and will calculate the ASDSFs at evenly spaced endpoints using samples `X[f(l):f(u)]`

| `u` | `l` | `(u - l + 1) / u` |
|---- | -- | ----------------- |
| 25  | 7  | 0.76  
| 50  | 13 | 0.76
| 75  | 20 | 0.747
| 100 | 26 | 0.75

and `f(i) = i + b`.

The output agrees with `makeplot.asdsf` when `window.lookback = 1` provided the endpoint indices are aligned.

---

##  Unit tests
To run all unit tests, execute the following from the `analyses/` directory.
```R
source("estimators.R")
source("coupling-functions.R")
testthat::test_dir("tests/testthat")
```
Add `reporter = c("summary", "debug")` as an argument to `test_*` to pause execution at failed `expect_*` instances; for example,
```R
testthat::test_dir("tests/testthat", reporter = c("summary", "debug"))
```
