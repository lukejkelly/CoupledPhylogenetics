addpath tools
GlobalSwitches; GlobalValues;

% Run script to get parameters
generateDataSets.paramConfig;

[grid_L, grid_root_time] = ndgrid(list_L, list_root_time);
[grid_lambda, grid_mu, grid_beta] = ndgrid(list_lambda, list_mu, list_beta);

if numel(grid_run_length) ~= 2 || numel(grid_sample_interval) ~= 2
    error('Edit makeSubmitFile to account for this');
end

dest_dir = generateDataSets.fileDest();
cellfun(@(x) mkdir(dest_dir, x), {'data', 'output', 'pars'});

for i = 1:numel(grid_L)
    [L, root_time] = deal(grid_L(i), grid_root_time(i));
    % We'll use the same tree for each parameter combination
    s = generateDataSets.sampleSyntheticTree(L, root_time);
    for j = 1:numel(grid_lambda)
        [lambda, mu, beta] = deal(grid_lambda(j), grid_mu(j), grid_beta(j));
        generateDataSets.sampleSyntheticData(s, lambda, mu, beta);
        for k = 1:numel(grid_run_length)
            [run_length, sample_interval] = deal(grid_run_length(k), ...
                                                 grid_sample_interval(k));
            generateDataSets.makeParFile(L, root_time, lambda, mu, beta, ...
                                         run_length, sample_interval);
        end
    end
end

generateDataSets.writeConfig(dest_dir, list_L, list_root_time, list_lambda, ...
    list_mu, list_beta, grid_run_length, grid_sample_interval);
generateDataSets.makeSubmitFile(dest_dir, list_L, list_root_time, ...
    list_lambda, list_mu, list_beta, grid_run_length, grid_sample_interval);
generateDataSets.makeJobFile(dest_dir, 'a');
generateDataSets.makeJobFile(dest_dir, 'b');
