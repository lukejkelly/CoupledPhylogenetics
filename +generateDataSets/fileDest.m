function [d, f, g] = fileDest(L, root_time, lambda, mu, beta, run_length, ...
        sample_interval)

    % Destination folder
    d = fullfile('..', 'CoupledPhylogeneticAnalyses', ...
                 num2str(yyyymmdd(datetime)));

    % For data sets
    if nargin >= 5
        f = sprintf('L%d_r%e_l%e_m%e_b%e', ...
                    L, root_time, lambda, mu, beta);
    else
        f = nan;
    end

    % For run files
    if nargin == 7
        g = sprintf('%s_n%e_s%e', f, run_length, sample_interval);
    else
        g = nan;
    end
end
