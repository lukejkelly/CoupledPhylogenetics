function s = sampleSyntheticTree(L, root_time)
    % Node times scaled so that root age is root_time
    global ROOT ANST

    THETA = -log(rand);
    s = ExpTree(L, THETA);

    % % Add a bunch of catastrophes
    % for i = 1:10
    %     j = randsample(find([s.type] <= ANST), 1);
    %     s(j).cat = s(j).cat + 1;
    %     s(j).catloc = [s(j).catloc, rand];
    % end

    % Rescale node times
    old_root_time = s([s.type] == ROOT).time;
    for j = find([s.type] == ANST | [s.type] == ROOT)
       s(j).time = s(j).time * root_time / old_root_time;
    end
end
