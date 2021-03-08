# CoupledPhylogenetics

Running experiments with marginal/coupled chains in `TraitLabSDLT-coupled`:
 * TODO: data/ contains Matlab code to generate synthetic data sets, parameter files and CEREMADE cluster submission scripts
   * edit +generateDataSets/paramConfig.m to configure the synthetic data generation
 * each experiment is a folder named <yyyymmdd>
    * data/ holds the .nex files and pars/ the corresponding .par files
    * output/ contains the output of running TraitLabSDLT-coupled and figs/ any figures
    * config.R is +generateDataSets/paramConfig.m in a format that R can read
    * execute submit.sh on the CEREMADE cluster to add experiments to the queue, job-a.pbs is for the coupled chains and job-b.pbs the long, single chains to get ground truth estimates
 * analyses/ is a mixture of scripts to analyse output and create figures in <experiment>/figs for a single experiment or figs/ for multiple experiments.

 The docs/ folder includes various documents we've created to describe our method. Some are now on Overleaf so I will delete them shortly.
