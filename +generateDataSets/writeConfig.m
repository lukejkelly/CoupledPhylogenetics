function writeConfig(dest_dir, list_L, list_root_time, list_lambda, ...
    list_mu, list_beta, grid_run_length, grid_sample_interval)


    fid = fopen(fullfile(dest_dir, 'config.R'), 'w');

    fprintf(fid, 'list_L <- c(%s)\n', format_list(list_L, '%d'));
    fprintf(fid, 'list_root_time <- c(%s)\n', format_list(list_root_time, '%e'));
    fprintf(fid, 'list_lambda <- c(%s)\n', format_list(list_lambda, '%e'));
    fprintf(fid, 'list_mu <- c(%s)\n', format_list(list_mu, '%e'));
    fprintf(fid, 'list_beta <- c(%s)\n', format_list(list_beta, '%e'));
    fprintf(fid, 'grid_run_length <- c(%s)\n', format_list(grid_run_length, '%e'));
    fprintf(fid, 'grid_sample_interval <- c(%s)\n', format_list(grid_sample_interval, '%e'));
    fprintf(fid, 'grid_c <- seq_len(100)\n');

    fclose(fid);

end

function fl = format_list(l, f)

    cl = arrayfun(@(x) sprintf(f, x), l, 'UniformOutput', false);
    fl = strjoin(cl, ', ');

end
