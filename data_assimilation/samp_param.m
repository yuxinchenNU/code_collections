%%Sample from a prior: provide minimum and maximum guesses in column vector
%%form, and a total number of particles, N, to sample.

%uniform prior
function [theta] = samp_param(min,max,N)
    dim = length(min);

    diff = max-min;

    min_N = repmat(min,1,N);
    diff_N = repmat(diff,1,N);

    theta = min_N + diff_N.*lhsdesign(dim,N);
