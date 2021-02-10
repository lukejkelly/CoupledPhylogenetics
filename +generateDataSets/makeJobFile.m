function makeJobFile(dest_dir, job_type)

    fid = fopen(fullfile(dest_dir, sprintf('job-%s.pbs', job_type)), 'w');

    switch job_type
    case 'a'
        batchCommand = 'batchTraitLabCoupled';
    case 'b'
        batchCommand = 'batchTraitLab';
    otherwise
        error('Job type must be a or b');
    end

    fprintf(fid, '#!/bin/bash\n\n');

    fprintf(fid, '#PBS -M kelly@ceremade.dauphine.fr\n');
    fprintf(fid, '#PBS -m e\n');
    fprintf(fid, '#PBS -l nodes=1:ppn=1\n');
    fprintf(fid, '#PBS -o log/kelly.o$PBS_JOBID\n');
    fprintf(fid, '#PBS -j oe\n\n');

    fprintf(fid, 'sleep $[ ($RANDOM %% 60) + 1 ]s\n\n');

    fprintf(fid, 'PAR_FILE=$(LC_NUMERIC="en_GB.UTF-8" \\\n');
    fprintf(fid, '           printf "L%%d_r%%e_l%%e_m%%e_b%%e_n%%e_s%%e" \\\n');
    fprintf(fid, '                  "$L" "$ROOT_TIME" "$LAMBDA" "$MU" "$BETA" "$RUN_LENGTH" \\\n');
    fprintf(fid, '                  "$SAMPLE_INTERVAL")\n\n');

    fprintf(fid, 'cd ~/TraitLabSDLT-coupled\n');
    fprintf(fid, '/usr/local/bin/matlab -nodesktop -nodisplay \\\n');
    fprintf(fid, '    -r "%s(''%s/${PAR_FILE}.par'', ${PBS_ARRAYID}); exit"\n', ...
            batchCommand, fullfile(dest_dir, 'pars'));

end
