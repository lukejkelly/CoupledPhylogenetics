function writeConfig(list_L, list_root_time, list_lambda, list_mu, ...
        list_beta, list_run_length, list_sample_interval, list_lag, n_chains)

    fid = fopen(fullfile(fileDest(), 'config.R'), 'w');

    fprintf(fid, 'list_L <- c(%s)\n', formatList(list_L, '%d'));
    fprintf(fid, 'list_root_time <- c(%s)\n', formatList(list_root_time, '%e'));
    fprintf(fid, 'list_lambda <- c(%s)\n', formatList(list_lambda, '%e'));
    fprintf(fid, 'list_mu <- c(%s)\n', formatList(list_mu, '%e'));
    fprintf(fid, 'list_beta <- c(%s)\n', formatList(list_beta, '%e'));
    fprintf(fid, ...
            'list_run_length <- list(marginal = %e, coupled = %e)\n', ...
            list_run_length.marginal, list_run_length.coupled);
    fprintf(fid, ...
            'list_sample_interval <- list(marginal = %e, coupled = %e)\n', ...
            list_sample_interval.marginal, list_sample_interval.coupled);
    fprintf(fid, 'list_lag <- c(%s)\n', formatList(list_lag, '%e'));
    fprintf(fid, 'list_c <- seq_len(%e)\n', n_chains);

    fclose(fid);

end

function fl = formatList(l, f)
    cl = arrayfun(@(x) sprintf(f, x), l, 'UniformOutput', false);
    fl = strjoin(cl, ', ');
end
