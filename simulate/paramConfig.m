% 20210805
list_L = 10;
list_root_time = 1e3;
list_lambda = 1e-1;
list_mu = 2.5e-04;
list_beta = 0; %, 2.5e-4];
list_run_length = struct('coupled', 1e6, 'marginal', 0);
list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
list_lag = [1e5, 3e5, 5e5];
n_chains = 100;
extras = struct('missing', 0, ...
                'clades', 1, ...
                'ncats', 0, ...
                'kappa', 1 / 3);
