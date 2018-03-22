function [baseline] = baseline(SPIKES, baseline_thresh)
%BASELINE evaluate basline of fluorescence traces
%
% SPIKES:           (N,M), N datapoints for M traces
% baseline thresh:  scalar, threshold for baseline calculation
% baseline:         (1,M), signal baseline
% 
%   giobarbera@neuralmapper.com


M = size(SPIKES,2);
N = size(SPIKES,1);

baseline = zeros(1,M);

zero_idx = zeros(N,M);

% Calculate baseline
for m = 1:M
    for n = 15:N-15
        if SPIKES(n,m) > baseline_thresh
            zero_idx(n-14:n+15,m) = ones(30,1);
        end
    end
end