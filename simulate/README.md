# Create experiment files with coupled and marginal chains

Sample synthetic data and create slurm submission scripts. The output can also be used as a template for real data, just change the paths to the data set, etc.

## Setup
The code uses relative path names so `TraitLab` and `CoupledPhylogenetics` must be in the same directory. The email address and other parameters of the slurm submission scripts need to be updated with your details.

## Create synthetic data experiment
Edit  `paramConfig` to set the parameters of the experiments. A data set will be generated for every combination of parameters. In the following example, for each tree a single data set will be generated and slurm files to run 100 coupled chains at each of lag 1e5, 3e5 and 5e5 for a minimum of 1e4 iterations, and 1 marginal chain for 1e6 iterations.
```matlab
list_L = [8, 12, 16]; % Number of taxa
list_root_time = 1e3; % Root time of synthetic tree
list_lambda = 1e-1; % Birth rate
list_mu = 2.5e-04; % Death rate
list_beta = 0; % Lateral transfer rate
list_run_length = struct('coupled', 1e4, 'marginal', 1e6);
list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
list_lag = [1e5, 3e5, 5e5]; % Lags for coupled experiments
n_chains = 100;
extras = struct('missing', 1, ... % Set these to 0 to disable
                'clades', 1, ...
                'ncats', 2, ...
                'kappa', 1 / 3);
```
Note: The `lag` must be an integer multiple of the `sample_interval`.

In this example, we do include lateral transfer so it is disabled in the corresponding `.par` files for TraitLab. If you do not want run a marginal chain then set its run_length to 0.

To create  all the files, start Matlab in this folder and execute `generateDataSets`.

1. This will load the contents of `paramConfig`.
2. Create a folder with today's date in the main `CoupledPhylogenetics` directory
3. Generate grids of parameters from the various lists
4. For each row in the (taxa number)-(root time) grid
    * Simulate a tree, catastrophes and and clades, etc.
        * Random catastrophe locations on the tree may cause identifiability issues, better pause execution in the debugger and place them on branches which are constrained by a clade
    * For each row in the SD model parameter grid (lambda, mu, beta)
        * Simulate synthetic data
        * Make marginal and coupled `.par` files
5. Write the experiment parameters to a file for `R` to read when analysing output
6. Make slurm job scripts and a bash script to submit them
    * Coupling is via a grid array
    * I use different numbers of cores when fitting models with and without lateral transfer

Some of the TraitLab settings are hard-coded into the template which is populated to create the `.par` files; fixing `mu` for example. As well as manually changing files, the template can be changed before calling `generateDataSets` or a tool such as `sed` can be used to edit them all at once.

## Run experiment

### Singe local experiment
Execute `batchTraitLab` in Matlab, with first argument the relevant `.par` file.

## Run all experiments on cluster

I use the `send_folder.sh` script to `rsync` a folder of submission scripts to the `CoupledPhylogenetics` folder on the compute cluster. My directory is structured as
```
TraitLab/
CoupledPhylogenetics/
    20210830/
    20210831/
    ...
```
If I want to run the experiments in `20210830`, then I `cd` to this folder and execute `bash submit.sh` to submit the jobs to the queue.

If you are rerunning jobs in the same folder then make sure to `rm output/* log/*`.

## Analyse output
See the README in the `analyses` directory for instructions and examples.
