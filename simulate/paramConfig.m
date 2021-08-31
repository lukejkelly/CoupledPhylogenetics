% % 20210815
% list_L = [8, 12, 16];
% list_root_time = 1e3;
% list_lambda = 1e-1;
% list_mu = 2.5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 0, 'marginal', 0);
% list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
% list_lag = [1e5, 3e5, 5e5];
% n_chains = 100;
% extras = struct('missing', 0, ...
%                 'clades', 0, ...
%                 'ncats', 0, ...
%                 'kappa', 1 / 10);

% % 20210818
% list_L = [8, 12, 16];
% list_root_time = 1e3;
% list_lambda = 1e-1;
% list_mu = 2.5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 0, 'marginal', 0);
% list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
% list_lag = [1e5, 3e5, 5e5];
% n_chains = 100;
% extras = struct('missing', 1, ...
%                 'clades', 1, ...
%                 'ncats', 3, ...
%                 'kappa', 1 / 20);

% 20210824: same as above but stronger kappa and more stringent cat prior
list_L = [8, 12, 16];
list_root_time = 1e3;
list_lambda = 1e-1;
list_mu = 2.5e-04;
list_beta = 0;
list_run_length = struct('coupled', 0, 'marginal', 0);
list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
list_lag = [1e5, 3e5, 5e5];
n_chains = 100;
extras = struct('missing', 1, ...
                'clades', 1, ...
                'ncats', 2, ...
                'kappa', 1 / 3);
