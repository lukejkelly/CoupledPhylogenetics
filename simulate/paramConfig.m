% 20210730
list_L = 10;
list_root_time = 1e3;
list_lambda = [1e-1, 2e-1];
list_mu = 5e-04;
list_beta = 0;
list_run_length = struct('coupled', 1e6, 'marginal', 3e6);
list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
list_lag = [1e5, 3e5, 5e5];
n_chains = 25;
extras = struct('missing', 1, ...
                'clades', 1, ...
                'ncats', 1, ...
                'kappa', 0.5);
