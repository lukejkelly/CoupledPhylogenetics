% Run script to get parameters
paramConfig;
if any(grid_lag < 0 | mod(grid_lag / grid_sample_interval.coupled, 1) ~= 0)
    error('Coupling lags must be an integer multiple of the sample interval');
end

% Tree parameter grids
[grid_L, grid_root_time] = ndgrid(list_L, list_root_time);
% Trait process parameter grids
[grid_lambda, grid_mu, grid_beta] = ndgrid(list_lambda, list_mu, list_beta);

dest_dir = fileDest();
cellfun(@(x) mkdir(dest_dir, x), {'data', 'output', 'pars'});

for i = 1:numel(grid_L)
    [L, root_time] = deal(grid_L(i), grid_root_time(i));
    % We'll use the same tree for each parameter combination
    s = sampleSyntheticTree(L, root_time);
    for j = 1:numel(grid_lambda)
        [lambda, mu, beta, lag] = deal(grid_lambda(j), grid_mu(j),
                                       grid_beta(j), grid_lag(j));
        sampleSyntheticData(s, lambda, mu, beta);
        for k = 1:numel(grid_run_length)
            [run_length, sample_interval] = deal(grid_run_length(k), ...
                                                 grid_sample_interval(k));
            if k == 1
                for lag = list_lag
                    makeParFile(L, root_time, lambda, mu, beta, run_length, ...
                                sample_interval, lag);
                end
            else
                makeParFile(L, root_time, lambda, mu, beta, run_length,
                            sample_interval);
            end
        end
    end
end

writeConfig(dest_dir, list_L, list_root_time, list_lambda, ...
    list_mu, list_beta, grid_run_length, grid_sample_interval, list_lag);
makeSubmitFile(dest_dir, list_L, list_root_time, ...
    list_lambda, list_mu, list_beta, grid_run_length, grid_sample_interval, list_lag);
makeJobFile(dest_dir, 'a');
makeJobFile(dest_dir, 'b');
