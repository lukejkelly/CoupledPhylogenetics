function makeJobFile(job_type)

    [dest_dir_local, dest_dir_traitlab] = fileDest();
    fid = fopen(fullfile(dest_dir_local, sprintf('job-%s.sh', job_type)), 'w');

    fprintf(fid, '#!/bin/bash\n\n');

    fprintf(fid, '#SBATCH --mail-user=kelly@ceremade.dauphine.fr\n');
    fprintf(fid, '#SBATCH --mail-type=END\n');
    fprintf(fid, '#SBATCH -N 1\n');
    fprintf(fid, '#SBATCH -n 1\n');

    switch job_type
    case 'a'
        fprintf(fid, '#SBATCH -o log/slurm-%%A_%%a.out\n');
        fprintf(fid, '#SBATCH -t 48:00:00\n\n');
    case 'b'
        fprintf(fid, '#SBATCH -o log/slurm-%%j.out\n');
        fprintf(fid, '#SBATCH -t 240:00:00\n\n');
    end

    fprintf(fid, 'sleep $[ ($RANDOM %% 60) + 1 ]s\n\n');

    % fprintf(fid, 'PAR_FILE=$(LC_NUMERIC="en_GB.UTF-8" \\\n');
    fprintf(fid, 'PAR_FILE=$( \\\n');
    switch job_type
    case 'a'
        fprintf(fid, '           printf "L%%d_r%%e_l%%e_m%%e_b%%e_l%%e" \\\n');
        fprintf(fid, '                  "$L" "$ROOT_TIME" "$LAMBDA" "$MU" "$BETA" "$LAG")\n\n');
    case 'b'
        fprintf(fid, '           printf "L%%d_r%%e_l%%e_m%%e_b%%e" \\\n');
        fprintf(fid, '                  "$L" "$ROOT_TIME" "$LAMBDA" "$MU" "$BETA")\n\n');
    end

    % Cluster job starts from home directory
    fprintf(fid, 'cd "$HOME"/TraitLabSDLT-coupled\n');
    fprintf(fid, '/usr/local/bin/matlab -nodesktop -nodisplay \\\n');
    switch job_type
    case 'a'
        fprintf(fid, ...
                '    -r "batchTraitLab(''%s/${PAR_FILE}.par'', $SLURM_ARRAY_TASK_ID); exit"\n', ...
                fullfile(dest_dir_traitlab, 'pars'));
    case 'b'
        fprintf(fid, ...
                '    -r "batchTraitLab(''%s/${PAR_FILE}.par''); exit"\n', ...
                fullfile(dest_dir_traitlab, 'pars'));
    end

end
