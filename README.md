# CoupledPhylogenetics

## Running experiments with marginal/coupled chains in `TraitLabSDLT-coupled`
The code uses  relative path names so `TraitLabSDLT-coupled` and `CoupledPhylogenetics` must be in the same directory.

### Set up experiment
Matlab code to generate synthetic data sets, parameter files and CEREMADE cluster submission scripts is in `simulate/`.

Edit `paramConfig.m` to configure the synthetic data generation and parameters for experiments

Each experiment is a folder named `<yyyymmdd>` (which is `.gitignore`d)
    * `data/` holds the `.nex` files and `pars/` the corresponding `.par` files
    * `output/` contains the output of running `TraitLabSDLT-coupled` and `figs/` any figures created by the scripts in `analyses/`
    * `config.R` is `make-experiments/paramConfig.m` in a format that `R` can read
    * change the PBS directives in `submit.sh`, if necessary
    * `job-a.pbs` is for the coupled chains and `job-b.pbs` the long, single chains to get ground truth estimates.

### Make additional edits, if necessary
If making multiple copies of a folder, one can use `sed` to update the references; for example,
```bash
cp -r 20210615 20210615a
cd 20210615a
grep -r -l 20210615 * | xargs sed -i "" 's/\(20210615\)/\1a/g'
```
on Mac, there is no argument to `-i` on Linux.

### Run experiment
From the directory of the experiment, execute `bash submit.sh` to submit the jobs to the queue.

If you are rerunning jobs in the same folder then make sure to `rm output/*` as `batchTraitLab` will not overwrite existing log files.

### Analyse output
The `analyses/` directory contains scripts to analyse output and create figures in `<experiment>/figs` for a single experiment, or `figs/` for multiple experiments.

## Documents
 The `docs/` folder includes various documents we've created to describe our method, such as the kernel which includes snippets of code from `TraitLabSDLT-coupled`. It also contains code for making Tikz figures to place on Overleaf (which does not seem to compile and use them online).
