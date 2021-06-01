function sampleSyntheticData(s, lambda, mu, beta, extras)
    % Sample synthetic data and write it to file
    global LEAF ROOT

    % Parameters of SDLT process
    L = length(s) / 2;
    root_node = [s.type] == ROOT;
    root_time = s(root_node).time;
    tr = [lambda; mu; beta; extras.kappa];

    % Read events and generate data
    [tEvents, rl] = stype2Events(s);
    D = simBorCatDeath(tEvents, tr);

    % Removing empty site-patterns
    D = D(sum(D, 2) > 0, :);

    % Masking matrix to incorporate missing data
    if extras.missing > 0
        disp('Generating missingness masks');
        xi = [s(rl).xi];
        M = (rand(size(D)) > repmat(xi, size(D, 1), 1));
        D(M) = 2;
    end

    % Adding data to leaves
    for l = 1:L
      s(rl(l)).dat = D(:, L + 1 - l)';
    end

    % Make some clades
    if extras.clades > 0
        fprintf('Synthesising %i clades\n', extras.clades);
        clade = synthclades(s, extras.clades, 2, 1 - rand^3);
    else
        clade = [];
    end

    [file_dir, ~, nex_stem] = fileDest(L, root_time, lambda, mu, beta);

    sFile = stype2nexus(s, '', 'BOTH', '', clade);
    fid = fopen(fullfile(file_dir, 'data', sprintf('%s.nex', nex_stem)), 'w');
    fprintf(fid, sFile);
    fclose(fid);

    if extras.ncats > 0
        cats = cellfun('length', {s.catloc});
        catTree = wnexcattree(s, root_node, cats(:));
        fid = fopen(fullfile(file_dir, 'data', sprintf('%s.cat', nex_stem)), 'w');
        fprintf(fid, catTree);
        fclose(fid);
    end

    if extras.missing > 0
        fid = fopen(fullfile(file_dir, 'data', sprintf('%s.xi', nex_stem)), 'w');
        fprintf(fid, 'name,xi\n');
        for j = find([s.type] == LEAF)
            fprintf(fid, '%s,%e\n', s(j).Name, s(j).xi);
        end
        fclose(fid);
    end
end
