function [d, e, f, g] = fileDest(L, root_time, lambda, mu, beta, lag)

    % Destination from current directory
    c = num2str(yyyymmdd(datetime));
    d = fullfile('..', c);

    % Destination from TraitLabSDLT-coupled
    e = fullfile('..', 'CoupledPhylogenetics', c);

    % For data sets
    if nargin >= 5
        f = sprintf('L%d_r%e_l%e_m%e_b%e', L, root_time, lambda, mu, beta);
    else
        f = nan;
    end

    % For run files
    if nargin == 5
        g = f;
    elseif nargin == 6
        g = sprintf('%s_l%e', f, lag);
    else
        g = nan;
    end
end
