function sampleSyntheticData(s, lambda, mu, beta)
    % Sample synthetic data and write it to file
    global ROOT

    % Parameters of SDLT process
    L = length(s) / 2;
    root_time = s([s.type] == ROOT).time;
    kappa = NaN;
    tr = [lambda; mu; beta; kappa];

    % Read events and generate data
    [tEvents, rl] = stype2Events(s);
    D = simBorCatDeath(tEvents, tr);

    % Removing empty site-patterns
    D = D(sum(D, 2) > 0, :);

    % % Masking matrix to incorporate missing data
    % xi = [s(rl).xi];
    % M = (rand(size(D)) > repmat(xi, size(D, 1), 1));
    % D(M) = 2;

    % Adding data to leaves
    for l = 1:L
      s(rl(l)).dat = D(:, L + 1 - l)';
    end

    % % Make some clades
    % global LEAF
    % i = randsample(find([s.type] == LEAF), 1);
    % j = randsample(setdiff(find([s.type] == LEAF), [i, s(i).sibling]), 1);
    % clades = {struct('name', 'dummy_clade', ...
    %                  'rootrange', [], ...
    %                  'adamrange', [], ...
    %                  'language', {{num2str(i), num2str(j)}})};
    % clades = {struct('name', 'root_clade', ...
    %                  'rootrange', [750, 1250], ...
    %                  'adamrange', [], ...
    %                  'language', {{'1', '2', '3', '4', '5'}})};

    % Write to Nexus file using this command
    % sFile = stype2nexus(s, '', 'BOTH', '', clades);

    [file_dir, ~, nex_stem] = fileDest(L, root_time, lambda, mu, beta);

    fid = fopen(fullfile(file_dir, 'data', sprintf('%s.nex', nex_stem)), 'w');
    fprintf(fid, stype2nexus(s, '', 'BOTH', ''));
    fclose(fid);
end
