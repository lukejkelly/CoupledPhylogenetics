function makeSubmitFile(list_L, list_root_time, list_lambda, list_mu, ...
        list_beta, list_lag, n_chains)

    job_name = '-N ${L}-${ROOT_TIME}-${LAMBDA}-${MU}-${BETA}';
    var_list = '-v L=${L},ROOT_TIME=${ROOT_TIME},LAMBDA=${LAMBDA},MU=${MU},BETA=${BETA}';

    pars_a = sprintf('-t 1-%d', n_chains);
    pars_b = '';

    fid = fopen(fullfile(fileDest(), 'submit.sh'), 'w');

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
    fprintf(fid, 'for LAG in%s; do\n', sprintf(' %e', list_lag));

    % Coupled chains get lag argument
    fprintf(fid, '%s', repmat(' ', 1, 24));
    fprintf(fid, ...
            'qsub %s-${LAG} %s %s,LAG=${LAG} job-a.pbs\n', ...
            job_name, pars_a, var_list);
    fprintf(fid, '%sdone\n', repmat(' ', 1, 20));

    % Ground truth longer chain
    fprintf(fid, '%s', repmat(' ', 1, 20));
    fprintf(fid, 'qsub %s %s %s job-b.pbs\n', job_name, pars_b, var_list);

    fprintf(fid, '%sdone\n', repmat(' ', 1, 16));
    fprintf(fid, '%sdone\n', repmat(' ', 1, 12));
    fprintf(fid, '%sdone\n', repmat(' ', 1, 8));
    fprintf(fid, '%sdone\n', repmat(' ', 1, 4));
    fprintf(fid, 'done\n');

    fclose(fid);

end
