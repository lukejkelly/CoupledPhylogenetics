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
