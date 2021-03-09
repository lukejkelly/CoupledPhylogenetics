function makeParFile(L, root_time, lambda, mu, beta, run_length, ...
        sample_interval, lag)

    coupling = exist('lag', 'var');
    if coupling
        [file_dir, nex_stem, par_stem] = fileDest(...
            L, root_time, lambda, mu, beta, run_length, sample_interval, lag);
    else
        lag = nan;
        [file_dir, nex_stem, par_stem] = fileDest(...
            L, root_time, lambda, mu, beta, run_length, sample_interval);
    end

    file_nex = sprintf('%s.nex', nex_stem);
    file_par = sprintf('%s.par', par_stem);

    fid = fopen(fullfile(file_dir, 'pars', file_par), 'w');

    fprintf(fid, '%% FULL PATH OF DATA FILE INCLUDE .NEX EXTENSION\n');
    fprintf(fid, 'Data_file_name = %s/data/%s\n', file_dir, file_nex);
    fprintf(fid, '\n');
    fprintf(fid, '%% ONE OF THE FOLLOWING THREE OPTIONS MUST BE ONE, THE OTHERS 0\n');
    fprintf(fid, 'Start_from_rand_tree = 1\n');
    fprintf(fid, 'Start_from_tree_in_output = 0\n');
    fprintf(fid, 'Start_from_true_tree = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% VALUE OF THETA IGNORED UNLESS Start_from_rand_tree == 1\n');
    fprintf(fid, 'Theta = 0.010000\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% NEXT TWO FIELDS IGNORED UNLESS Start_from_tree_in_output == 1\n');
    fprintf(fid, '%% FULL PATH OF OLD OUTPUT FILE INCLUDE .NEX EXTENSION\n');
    fprintf(fid, 'Tree_file_name =\n');
    fprintf(fid, 'Use_tree = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Omit_taxa = 0\n');
    fprintf(fid, '%% LIST IS IGNORED UNLESS Omit_taxa == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fprintf(fid, 'Omit_taxa_list =\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Omit_traits = 0\n');
    fprintf(fid, '%% LIST IS IGNORED UNLESS Omit_traits == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fprintf(fid, 'Omit_trait_list =\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% ONE OF THE FOLLOWING TWO OPTIONS MUST BE 1, THE OTHER 0\n');
    fprintf(fid, 'Yule_prior_on_tree = 0\n');
    fprintf(fid, 'Flat_prior_on_tree = 1\n');
    fprintf(fid, '%% ONE OF THE FOLLOWING TWO OPTIONS MUST BE 1, THE OTHER 0\n');
    fprintf(fid, 'Uniform_prior_on_tree_topologies = 1\n');
    fprintf(fid, 'Uniform_prior_on_labelled_histories = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED UNLESS Flat_prior_on_tree == 1\n');
    fprintf(fid, 'Max_root_age = 10000\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Vary_topology = 1\n');
    fprintf(fid, 'Account_rare_traits = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Impose_clades = 0\n');
    fprintf(fid, '%% LIST IS IGNORED UNLESS Impose_clades == 1 CAN USE MATLAB VECTOR NOTATION\n');
    fprintf(fid, 'Omit_clade_list =\n');
    fprintf(fid, 'Omit_clade_ages_list =\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Vary_loss_rate = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Random_initial_loss_rate == 1\n');
    fprintf(fid, 'Initial_loss_rate = %g\n', LossRate(mu));
    fprintf(fid, 'Random_initial_loss_rate = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Account_for_lateral_transfer = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Account_for_lateral_transfer == 0\n');
    fprintf(fid, 'Vary_borrowing_rate = 1\n');
    fprintf(fid, 'Random_initial_borrowing_rate = 1\n');
    fprintf(fid, '%% NEXT LINE IS IGNORED WHEN Random_initial_borrowing_rate == 1\n');
    fprintf(fid, 'Initial_borrowing_rate = %g\n', beta);
    fprintf(fid, '\n');
    fprintf(fid, 'Include_catastrophes = 0\n');
    fprintf(fid, '%% NEXT 6 LINES ARE IGNORED WHEN Include_catastrophes = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Random_initial_cat_death_prob = 1\n');
    fprintf(fid, 'Initial_cat_death_prob = 0.5\n');
    fprintf(fid, 'Random_initial_cat_death_prob = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Random_initial_cat_rate = 1\n');
    fprintf(fid, 'Initial_cat_rate = 0.00015\n');
    fprintf(fid, 'Random_initial_cat_rate = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Model_missing = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, 'Run_length = %g\n', run_length);
    fprintf(fid, 'Sample_interval = %g\n', sample_interval);
    fprintf(fid, '\n');
    fprintf(fid, 'Seed_random_numbers = 0\n');
    fprintf(fid, '%% FOLLOWING IS IGNORED UNLESS Seed_random_numbers == 1\n');
    fprintf(fid, 'With_seed = 0\n');
    fprintf(fid, '\n');
    fprintf(fid, '%% OUTPUT FILE NAME OMITTING PATH AND ANY EXTENSIONS\n');
    fprintf(fid, 'Output_file_name = %s\n', par_stem);
    fprintf(fid, '%% FULL PATH FOR DIRECTORY TO OUTPUT FILES\n');
    fprintf(fid, 'Output_path_name = %s/output/\n', file_dir);
    fprintf(fid, '\n');
    fprintf(fid, '%% Gaps are treated as missing data. To change this, edit GlobalValues.m.\n\n');

    % Coupling and lags LJK 19/2/21
    fprintf(fid, 'Coupled_markov_chains = %g\n', coupling);
    fprintf(fid, '%% FOLLOWING IS IGNORED WHEN Coupled_markov_chains == 0\n');
    fprintf(fid, 'Coupling_lag = %g\n', lag);

    fclose(fid);

end
