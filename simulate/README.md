# Create experiment files with coupled and marginal chains

Sample synthetic data and create slurm submission scripts. The output can also be used as a template for real data, just change the paths to the data set, and so on.

## Setup
The code uses relative path names so `TraitLab` and `CoupledPhylogenetics` must be in the same directory. The email address and other parameters of the slurm submission scripts need to be updated with your details, either by editing `makeJobFile` before executing `generateDataSets` or editing the run files after they have been generated.

## Create synthetic data experiment
Edit  `paramConfig` to set the parameters of the experiments. A data set will be generated for every combination of parameters. In the following example, for each tree a single data set will be generated and slurm files to run 100 coupled chains at each of lags 1e5, 3e5 and 5e5, for a minimum of 1e4 iterations, and 1 marginal chain for 1e6 iterations.
```matlab
list_L = [8, 12, 16]; % Number of taxa
list_root_time = 1e3; % Root time of synthetic tree
list_lambda = 1e-1; % Birth rate
list_mu = 2.5e-4;  % Death rate
list_beta = 0; % Lateral transfer rate
list_run_length = struct('coupled', 1e4, 'marginal', 1e6);
list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
list_lag = [1e5, 3e5, 5e5];
n_chains = 100;
extras = struct(... %  Set extra components to 0 to disable
    'missing', 0, ...
    'clades', 1, ... % Must be non-zero if catastrophes are included
    'ncats', 1, ... % Placed on branches constrained by clades
    'kappa', 1 / 3 ...
);
```
Note: The `lag` must be an integer multiple of the `sample_interval`.

In this example, we do include lateral transfer so it is disabled in the corresponding `.par` files for TraitLab. If you do not want run a marginal chain then set its `run_length` to 0.

To create all the files, start Matlab in this folder, `CoupledPhylogenetics/simulate` and execute `generateDataSets`.

1. A `startup` scripts adds TraitLab to the path.
2. This will load the contents of `paramConfig`.
3. Create a folder with today's date in the main `CoupledPhylogenetics` directory
4. Generate grids of parameters from the various lists
5. For each row in the (taxa number)-(root time) grid
    * Simulate a tree, catastrophes and and clades, etc.
        * Random catastrophe locations on branches which have clade constraints
    * For each row in the SD model parameter grid (lambda, mu, beta)
        * Simulate synthetic data
        * Make marginal and coupled `.par` files
6. Write the experiment parameters to a file for R to read when analysing output
7. Make slurm job scripts and a bash script to submit them
    * Coupled analyses are run as a grid array
    * The scripts specify a different numbers of cores for each job when fitting models with (2 cores, submission script uses `job-a2.sh` for coupled runs and `job-b2.sh` for marginal runs) and without (1 core, `job-a1.sh` and `job-b1.sh`) lateral transfer

Some of the TraitLab settings are hard-coded into the template which is populated to create the `.par` files; fixing `mu` for example. As well as manually changing files, the template can be changed before calling `generateDataSets`.

## Run experiment

### Single local experiment
Execute `batchTraitLab` in Matlab, with first argument the relevant `.par` file.

### Run all experiments on cluster

My directory is structured as
```
TraitLab/
CoupledPhylogenetics/
    20210830/
    20210831/
    ...
```
If I want to run the experiments in `20210830`, then I `cd` to this folder on the cluster and execute `bash submit.sh` to submit the jobs to the queue.

If you are rerunning jobs in the same folder then make sure to `rm output/* log/*` first.

## Analyse output
See the README in the `analyses` directory for instructions and examples.
