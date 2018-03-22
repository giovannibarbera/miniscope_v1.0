function [N_spikes, intertimes] = count_spikes(bin_spikes)
%COUNT_SPIKES Returns the number of positive edges in a binary array
%representing single cell firing
%
% bin_spikes:       (N_frames), binary array of spikes
%
% N_SPIKES:         scalar, total number of spikes  
% intertimes        (N_spikes), delay between spikes
% 
%   giobarbera@neuralmapper.com

N_spikes = 0;
intertimes = 0;
size_spk = size(bin_spikes);
prev = 0;
cnt = 0;

for k = 1:size_spk
    if (bin_spikes(k) == 1) && (prev == 0)
        N_spikes = N_spikes + 1;
        if N_spikes>1 % skip first spike
            intertimes = [intertimes; cnt];
        end
        cnt = 0;
    end
    
    cnt = cnt + 1;
    prev = bin_spikes(k);
    
   
end

if length(intertimes) >1
    intertimes = intertimes(2:end);
end



