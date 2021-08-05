rng(yyyymmdd(datetime));
% Global variables from TraitLab
GlobalSwitches
GlobalValues

% Run script to get parameters
paramConfig;
if any(list_lag < 0 | mod(list_lag / list_sample_interval.coupled, 1) ~= 0)
    error('Coupling lags must be an integer multiple of the sample interval');
end

dest_dir = fileDest();
cellfun(@(x) mkdir(dest_dir, x), {'data', 'output', 'pars', 'log'});

% Tree parameter grids
[grid_L, grid_root_time] = ndgrid(list_L, list_root_time);
% Trait process parameter grids
[grid_lambda, grid_mu, grid_beta] = ndgrid(list_lambda, list_mu, list_beta);

% Generate a tree for each grid element and data for each trait grid
for i = 1:numel(grid_L)
    % We'll use the same tree for each parameter combination
    treePars = {grid_L(i), grid_root_time(i)};
    fprintf('Sampling tree with %i leaves and root at %e\n', treePars{:});
    [s, clade] = sampleSyntheticTree(treePars{:}, extras);
    for j = 1:numel(grid_lambda)
        % Sample and write data
        traitPars = {grid_lambda(j), grid_mu(j), grid_beta(j)};
        fprintf('Data when lambda = %e, mu = %e, beta = %e\n', traitPars{:});
        sampleSyntheticData(s, traitPars{:}, extras, clade);
        % Write coupled and marginal run files
        arrayfun(@(lag) makeParFile(treePars{:}, traitPars{:}, ...
                                    list_run_length.coupled, ...
                                    list_sample_interval.coupled, lag, ...
                                    extras), ...
                 list_lag);
        makeParFile(treePars{:}, traitPars{:}, list_run_length.marginal, ...
                    list_sample_interval.marginal, [], extras);
    end
end
% Write file with settings in R
writeConfig(list_L, list_root_time, list_lambda, list_mu, list_beta, ...
            list_run_length, list_sample_interval, list_lag, n_chains, extras);
% Bash script to submit to cluster
makeSubmitFile(list_L, list_root_time, list_lambda, list_mu, list_beta, ...
               list_lag, n_chains);
% PBS files for coupled (a) and marginal (b) experiments
makeJobFile('a1');
makeJobFile('a2');
makeJobFile('b1');
makeJobFile('b2');
