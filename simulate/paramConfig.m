% % 20210218
% list_L = [8, 10];
% list_root_time = 1e3;
% list_lambda = [2.5e-2, 5e-2];
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 5e5, 'marginal', 5e6);
% list_sample_interval = struct('coupled', 5e2, 'marginal', 5e2);
% list_lag = [5e2, 5e3, 5e4, 5e5];
% n_chains = 1e2;

% % 20210311
% list_L = [8];
% list_root_time = 1e3;
% list_lambda = [2.5e-2];
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 2.5e5, 'marginal', 2.5e6);
% list_sample_interval = struct('coupled', 2.5e2, 'marginal', 2.5e2);
% list_lag = [2.5e2, 2.5e3, 2.5e4];
% n_chains = 1e2;

% % 20210319
% list_L = [6];
% list_root_time = 1e3;
% list_lambda = [2.5e-2];
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 5e4, 'marginal', 5e6);
% list_sample_interval = struct('coupled', 5e1, 'marginal', 5e2);
% list_lag = [5e1, 5e2, 5e3];
% n_chains = 1e2;

% % 20210406
% list_L = 8;
% list_root_time = 1e3;
% list_lambda = 2.5e-2;
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 0, 'marginal', 2.5e6);
% list_sample_interval = struct('coupled', 1, 'marginal', 2.5e2);
% list_lag = [1, 1e2, 1e4];
% n_chains = 1e2;

% % 20210408
% list_L = [8, 10];
% list_root_time = 1e3;
% list_lambda = 2.5e-2;
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 5e5, 'marginal', 5e6);
% list_sample_interval = struct('coupled', 5e2, 'marginal', 5e2);
% list_lag = [5e4, 1e5];
% n_chains = 1e2;

% % 20210428{a,b}{1,2}
% list_L = [10];
% list_root_time = 1e3;
% list_lambda = [1, 2, 3, 4];
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 1e6, 'marginal', 1e6);
% list_sample_interval = struct('coupled', 1e3, 'marginal', 1e3);
% list_lag = 1e5;
% n_chains = 10;

% % 20210601
% list_L = [8, 10];
% list_root_time = 1e3;
% list_lambda = [2.5e-2, 5e-2];
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 1e6, 'marginal', 1e6);
% list_sample_interval = struct('coupled', 1e3, 'marginal', 1e3);
% list_lag = 1e5;
% n_chains = 10;
% extras = struct('missing', 1, ...
%                 'clades', 1, ...
%                 'ncats', 2, ...
%                 'kappa', 0.25 + rand * 0.75);

% % 20210603
% list_L = 8;
% list_root_time = 1e3;
% list_lambda = 2.5e-2;
% list_mu = 5e-04;
% list_beta = list_mu;
% list_run_length = struct('coupled', 1e6, 'marginal', 1e6);
% list_sample_interval = struct('coupled', 1e3, 'marginal', 1e3);
% list_lag = 1e5;
% n_chains = 10;
% extras = struct('missing', 1, ...
%                 'clades', 1, ...
%                 'ncats', 0, ...
%                 'kappa', 0.25 + rand * 0.75);

% % 20210604
% list_L = 8;
% list_root_time = 1e3;
% list_lambda = 2.5e-2;
% list_mu = 5e-04;
% list_beta = 0;
% list_run_length = struct('coupled', 0, 'marginal', 1e6);
% list_sample_interval = struct('coupled', 1e1, 'marginal', 1e3);
% list_lag = 1e5;
% n_chains = 20;
% extras = struct('missing', 1, ...
%                 'clades', 1, ...
%                 'ncats', 5, ...
%                 'kappa', 0.25 + rand * 0.75);

% 20210619
list_L = 10;
list_root_time = 1e3;
list_lambda = 2e-1;
list_mu = 5e-04;
list_beta = 0;
list_run_length = struct('coupled', 1e6, 'marginal', 1e6);
list_sample_interval = struct('coupled', 1e2, 'marginal', 1e2);
list_lag = 1e5;
n_chains = 20;
extras = struct('missing', 1, ...
                'clades', 1, ...
                'ncats', 2, ...
                'kappa', 0.221199);
