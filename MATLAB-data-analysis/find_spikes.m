function [OUTPUT] = find_spikes(INPUT, pos_thresh, neg_thresh, delay, use_neg_delay)
%FIND_SPIKES Calculate onset of spikes (above or below threshold) and
%returns a vector with positive and negative Kronecker deltas.
%
% INPUT:            (N,M), input, N points, M traces
% pos_thresh:       (1), positive threshold for spike detection
% neg_thresh:       (1), negative threshold for spike detection
% delay:            (1), average time instants of spike onset before threshold is reached
% use_neg_delay:    (1), set to 1 to use negative delay with negative spikes
%
% OUTPUT:           (N,M), output spike onsets
% 
%   giobarbera@neuralmapper.com


M = size(INPUT,2);
N = size(INPUT,1);

OUTPUT = zeros(N,M);

spike_idx_p = false(N,M);
spike_idx_n = false(N,M);


% Calculate spike events 
for m = 1:M
    for n = delay:N-2*delay
        if INPUT(n,m) > pos_thresh
            spike_idx_p(n-delay+1:n+delay,m) = true(2*delay,1);
        end
        if INPUT(n,m) < neg_thresh
            spike_idx_n(n-delay+1:n+delay,m) = true(2*delay,1);
        end
    end
end

% figure; plot(spike_idx_p); hold on; plot(-spike_idx_n,'r');


% Calculate spike onsets
POS_EDGES = false(N,M);
NEG_EDGES = false(N,M);
for m = 1:M
    tmp_p = [spike_idx_p(2:end,m); NaN];
    tmp_p2 = find(spike_idx_p(:,m) < 0.5 & tmp_p > 0.5);
    % tmp_p2 = tmp_p2+((delay).*ones(length(tmp_p2),1));
    if ~isempty(tmp_p2)
        if tmp_p2(1) == 0 % could happen
            tmp_p2(1) = 1;
        end
    end
    
    % N_spikes(m) = size(tmp2,1);
    %if ~isempty(tmp_p2)
        POS_EDGES(tmp_p2,m)=true;
    %end
    
    tmp_n = [spike_idx_n(2:end,m); NaN];
    tmp_n2 = find(spike_idx_n(:,m) < 0.5 & tmp_n > 0.5);
    
    %Use negative delay
    if use_neg_delay == 1
        tmp_n2 = tmp_n2+((2*delay).*ones(size(tmp_n2)));
    end
    
    if ~isempty(tmp_n2)
        if tmp_n2(1)==0
            tmp_n2(1) = 1;
        end
    end
    % N_spikes(m) = size(tmp2,1);
    %if ~isempty(tmp_n2)
        NEG_EDGES(tmp_n2,m)=true;
    %end
   
end

OUTPUT = POS_EDGES-NEG_EDGES;