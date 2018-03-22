function [N_spikes, activity, spike_amp, spike_freq, baseline, SPIKES2, SYNC] = neural_stat(SPIKES, bin_size, baseline_thresh, event_thresh)
%NEURAL_STAT evaluate statistics of multiple calcium traces
%
% SPIKES:           (N,M), N datapoints for M traces
%
% N_spikes:         (1,M), total number of spikes per trace
% bin_size:         scalar, window bin for frequency calculation
% baseline thresh:  scalar, threshold for baseline calculation
% event_thresh:     scalar, threshold for transient identification
% activity:         (1,M), integral sum of each calcium trace
% spike_amp:        (M,P), average bin transient amplitude
% spike_freq:       (M,P), average bin transient frequency
% baseline:         (1,M), signal baseline
% SPIKES2:          (N,M), baseline removed calcium traces
% SYNC:             (1,N), neurons active simultaneously at each time
% 
%   giobarbera@neuralmapper.com


M = size(SPIKES,2);
N = size(SPIKES,1);
P = floor(N/bin_size);

N_spikes = zeros(1,M);
activity = zeros(1,M);
spike_amp = zeros(M,P);
spike_freq = zeros(M,P);
baseline = zeros(1,M);
SPIKES2 = zeros(N,M);
SYNC = zeros(1,N);


zero_idx = zeros(N,M);
spike_idx = zeros(N,M);

% Calculate spike events and baseline

% % First calculate maximum spike amplitude
% mn = zeros(1,M);
% tmp_prev = 0;
% tmp_val = 0;
% tmp_max = 0;
% for m = 1:M  
%     for n = 1:N
%         if SPIKES(n,m) >= tmp_prev
%             tmp_val = tmp_val + (SPIKES(n,m)-tmp_prev);
%         else
%             if tmp_val > mn(m)
%                 mn(m) = max(tmp_val, mn(m));
%                 tmp_max = SPIKES(n,m);
%             end
%             
%             tmp_val = 0;
%         end
%         tmp_prev = SPIKES(n,m);
%     end
%     mn(m) = min(SPIKES(:,m)) + ((tmp_max-mn(m)) - min(SPIKES(:,m)))/2  ;
%     SPIKES2(:,m) = SPIKES(:,m) - mn(m);
% end
% 
% % Then subtract it to the spike values and calculate spikes and baseline
% % vectors
% for m = 1:M
%     for n = 30:N-30
%         if abs(SPIKES2(n,m)) > baseline_thresh
%             zero_idx(n-29:n+30,m) = ones(60,1);
%         end
%         if (SPIKES2(n,m)) > event_thresh
%             spike_idx(n-4:n+5,m) = ones(10,1);
%         end
%     end
% end




% First calculate baseline when signal is steady (within +-5)
% for 3 seconds

mn = zeros(1,M);
tt = 50;
thresh = 3;
tmp0 = 0;
tmp_prev = 0;

for m = 1:M  
    for n = tt:N
        if abs(SPIKES(n,m) - tmp_prev) < thresh
            tmp0 = tmp0 + 1;
        else
            tmp0 = 0;
        end
        tmp_prev = SPIKES(n,m);
        
        mn(m) = mean(SPIKES(:,m));
        
        if tmp0 > tt
            if (max(SPIKES(n-tt+1:n,m)-SPIKES(n-round(tt/2),m)) < thresh) & (min(SPIKES(n-tt+1:n,m)-SPIKES(n-round(tt/2),m)) > -thresh)
              %  if mn(m) == 0
                   mn(m) = mean(SPIKES(n-tt+1:n,m));
                   break;
              %  else
              %     mn(m) = mean([mn(m),mean(SPIKES(n-tt+1:n,m))]);
              %  end
            end
        end
        
    end
    SPIKES2(:,m) = SPIKES(:,m) - mn(m);
end






% Then subtract it to the spike values and calculate spikes and baseline
% vectors
for m = 1:M
    for n = 30:N-30
        if abs(SPIKES2(n,m)) > baseline_thresh
            zero_idx(n-29:n+30,m) = ones(60,1);
        end
        if (SPIKES2(n,m)) > event_thresh
            spike_idx(n-4:n+5,m) = ones(10,1);
        end
    end
end





% Remove baseline signal
for m = 1:M
    if length(find(zero_idx(:,m)==0))>0
        baseline(m) = mean(SPIKES2(find(zero_idx(:,m)==0),m));
        SPIKES2(:,m) = SPIKES2(:,m)-baseline(m);
        % SPIKES2(find(zero_idx(:,m)==0),m) = 0;
    end
        
end
% Remove negative spikes
% SPIKES2(find(SPIKES2<0))=0;

% figure; plot(SPIKES2(:,60)); hold on; plot(50.*zero_idx(:,60));  plot(30.*spike_idx(:,60));





% Caclulate cell activity
for m = 1:M
    activity(m) = sum(SPIKES2(:,m));
end


% Calculate number of calcium transients

POS_EDGES = zeros(N,M);
for m = 1:M
    tmp = [spike_idx(2:end,m); NaN];
    tmp2 = find(spike_idx(:,m) < 0.5 & tmp > 0.5);
    N_spikes(m) = size(tmp2,1);
    POS_EDGES(tmp2,m)=1;
end




% Calculate frequency
last_spike = 0;
for m = 1:M
    last_spike = 0;
    for p = 1:P
        % spike_freq(m,p) = sum(POS_EDGES((p-1)*bin_size+1:p*bin_size,m));
  
        spikes = find(POS_EDGES((p-1)*bin_size+1:p*bin_size,m)==1)+(p-1)*bin_size;

        for h = 1:length(spikes)
            if h == 1
                spike_freq(m,p) = spike_freq(m,p) + (spikes(h)-last_spike);
            else
                spike_freq(m,p) = spike_freq(m,p) + (spikes(h)- spikes(h-1));
            end
            last_spike = spikes(h); 
        end
        if isempty(spikes)==0
            spike_freq(m,p) = 60*10./(spike_freq(m,p)./length(spikes));
        end
       
    end
end

% MM = 8;
% figure; plot(SPIKES2(:,MM)); hold on; plot(SPIKES(:,MM),'r'); plot(40.*POS_EDGES(:,MM),'g--'); % plot(40.*spike_idx(:,MM),'r--');
% figure; plot(spike_freq(MM,:));
% spike_freq(MM,:)

% Caclulate calcium events amplitudes
for m = 1:M
    for p = 1:P
        spike_amp(m,p) = max(SPIKES2((p-1)*bin_size+1:p*bin_size,m));
    end
end

for n = 1:N
    SYNC(n) = sum(spike_idx(n,:));
end

% figure; plot(SPIKES2(:,13)); hold on; plot(30.*zero_idx(:,13),'r'); plot(50.*spike_idx(:,13),'g');
% figure; plot(N_spikes);
% figure; plot(spike_freq(13,:));
% figure; plot(spike_amp(13,:));

