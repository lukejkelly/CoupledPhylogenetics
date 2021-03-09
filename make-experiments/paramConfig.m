% 20210218
list_L = [8, 10];
list_root_time = 1e3;
list_lambda = [2.5e-2, 5e-2];
list_mu = 5e-04;
list_beta = 0;
grid_run_length = struct('coupled', 5e5, 'marginal', 5e6);
grid_sample_interval = struct('coupled', 5e2, 'marginal', 5e2);
grid_lag = [5e2, 5e3, 5e4, 5e5];
