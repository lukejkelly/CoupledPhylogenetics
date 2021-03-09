function makeSubmitFile(dest_dir, list_L, list_root_time, list_lambda, ...
        list_mu, list_beta, grid_run_length, grid_sample_interval, grid_lag)

    job_name = '-N ${L}-${ROOT_TIME}-${LAMBDA}-${MU}-${BETA}';
    var_list = '-v L=${L},ROOT_TIME=${ROOT_TIME},LAMBDA=${LAMBDA},MU=${MU},BETA=${BETA}';

    pars_a = '-t 1-100 -l walltime=24:00:00';
    pars_b = '-t 0 -l walltime=240:00:00';

    fid = fopen(fullfile(dest_dir, 'submit.sh'), 'w');

    fprintf(fid, '#!/bin/bash\n');

    fprintf(fid, '\n');
    fprintf(fid, 'for L in%s; do\n', sprintf(' %d', list_L));

    fprintf(fid, '%s', repmat(' ', 1, 4));
    fprintf(fid, 'for ROOT_TIME in%s; do\n', sprintf(' %e', list_root_time));

    fprintf(fid, '%s', repmat(' ', 1, 8));
    fprintf(fid, 'for LAMBDA in%s; do\n', sprintf(' %e', list_lambda));

    fprintf(fid, '%s', repmat(' ', 1, 12));
    fprintf(fid, 'for MU in%s; do\n', sprintf(' %e', list_mu));

    fprintf(fid, '%s', repmat(' ', 1, 16));
    fprintf(fid, 'for BETA in%s; do\n', sprintf(' %e', list_beta));

    fprintf(fid, '%s', repmat(' ', 1, 20));
    for lag = grid_lag
        fprintf(fid, ...
                'qsub %s %s %s,RUN_LENGTH=%e,SAMPLE_INTERVAL=%e,LAG=%e job-a.pbs\n', ...
                job_name, pars_a, var_list, grid_run_length(1), ...
                grid_sample_interval(1), lag);
    end
    fprintf(fid, '%s', repmat(' ', 1, 20));
    fprintf(fid, ...
            'qsub %s %s %s,RUN_LENGTH=%e,SAMPLE_INTERVAL=%e job-b.pbs\n', ...
            job_name, pars_b, var_list, grid_run_length(2), ...
            grid_sample_interval(2));

    fprintf(fid, '%sdone\n', repmat(' ', 1, 16));
    fprintf(fid, '%sdone\n', repmat(' ', 1, 12));
    fprintf(fid, '%sdone\n', repmat(' ', 1, 8));
    fprintf(fid, '%sdone\n', repmat(' ', 1, 4));
    fprintf(fid, 'done\n');

    fclose(fid);

end
