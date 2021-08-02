function [s, clade] = sampleSyntheticTree(L, root_time, extras)
    % Node times scaled so that root age is root_time
    global ROOT ANST

    % Sample tree
    % theta = -log(rand);
    % s = ExpTree(L, theta);
    [s, ~, ~, clade] = nexus2stype('../../TraitLabSDLT/data/SIM_B.nex');

    % Rescale node times
    old_root_time = s([s.type] == ROOT).time;
    for j = find([s.type] == ANST | [s.type] == ROOT)
       s(j).time = s(j).time * root_time / old_root_time;
    end

    % % Make some clades
    if extras.clades > 0
    %     fprintf('Synthesising %i clades\n', extras.clades);
    %     clade = synthclades(s, extras.clades, 2, 1 - rand^3);
        rootmax = 2e3; % rootmax = (1 + rand) * s([s.type] == ROOT).time;
        prior = unitTests.clade2prior(clade, rootmax);
        s = treeclades(s, prior.clade);
    else
        clade = [];
    end

    % Sample catastrophes
    if extras.ncats > 0
        fprintf('Sampling %i catastrophe locations\n', extras.ncats);
        bInds = find(~cellfun(@isempty, {s.clade}));
        % bLengths = arrayfun(@(i) s(s(i).parent).time - s(i).time, bInds);
        % cInds = bInds(randsample(length(bInds), extras.ncats, true, bLengths));
        % for j = cInds
        %     s(j).cat = s(j).cat + 1;
        %     s(j).catloc = [s(j).catloc, rand];
        % end
        j7 = bInds(arrayfun(@(j) strcmp(s(j).Name, '7'), bInds));
        p67 = s(j7).parent;
        s(p67).cat = 1;
        s(p67).catloc = 2 / 3;

        j9 = bInds(arrayfun(@(j) strcmp(s(j).Name, '9'), bInds));
        s(j9).cat = 1;
        s(j9).catloc = 1 / 3;
    end
end
