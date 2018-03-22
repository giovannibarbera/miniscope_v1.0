
% close all; 
% clear all; 
% clc


%% Load single mouse data


BASEDIR             = '..\data'; 


cell_type = 1;
mouse_numb = 3;



% workspace_name = 'D736_042716';  % ***
% workspace_name = 'D736';

eval(['MOUSE = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL;'])
% eval(['MOUSE = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),';'])

[FRAMES, N_days, N_ss, T, width, height, days, sessions, cell_mask_size, cell_size, min_cell_distance, rms_fact, area_min, T_cell, N_det, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, X_SHIFT, Y_SHIFT, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CENT_X_DAY, CENT_Y_DAY, CENT_X_GLOBAL, CENT_Y_GLOBAL, IDX_MAP, CA_TRACES_RAW, CA_TRACES_FILT, CA_TRACES_BIN, CA_TRACES_ONSET, CA_TRACES_RAW_G, CA_TRACES_FILT_G, CA_TRACES_BIN_G, CA_TRACES_ONSET_G, clust_idx, clust_idx_s, sort_idx,  CLUST_CONN_FULL, clust_conn_thresh, clust_iterations, clust_bin_size, clust_curr_days, clust_curr_sessions, clust_metric, clust_tot_kmeans_iter, clust_max_iter, global_clusters, eartag, ref_day, ref_angle] = load_mouse_full(MOUSE);
clear MOUSE
colors;


fprintf('Loaded! Ready to go.\n');





%% Set parameters


 

% Butterworth filter parameters
[butter_b,butter_a] = butter(16,0.33,'low'); 

bin_thresh = 3;
onset_thresh = 6;
noise_thresh = 2;

  
max_delta_F = 60;

bin_size = 50;
baseline_thresh = 3;
event_thresh = 5;       % 5 percent deltaF/F
mov_avg = 500;               % moving average 





% % For 1000um GRIN lens
% r_roi = 4;   
% R1 = 7;        
% R2 = 11;


% For 1000um GRIN lens
r_roi = 4;   
R1 = 5;        
R2 = 9;

% % % For 1000um GRIN lens
%  r_roi = 5;   
%  R1 = 7;        
%  R2 = 13;




ROI = false(2*r_roi,2*r_roi);
for m = -r_roi:r_roi
    for n = -r_roi:r_roi
         if  m^2 + n^2 <= r_roi^2
            ROI(m+r_roi+1,n+r_roi+1) = true;
         end
    end
end




% % For 500um GRIN lens
% r_roi = 7;   
% R1 = 8;        
% R2 = 18;
% 
% ROI = false(2*r_roi,2*r_roi);
% for m = -r_roi:r_roi
%     for n = -r_roi:r_roi
%          if  m^2 + n^2 <= r_roi^2
%             ROI(m+r_roi+1,n+r_roi+1) = true;
%          end
%     end
% end



% % For 500um GRIN lens
% r_roi = 8;   
% R1 = 10;        
% R2 = 16;
% 
% ROI = false(2*r_roi,2*r_roi);
% for m = -r_roi:r_roi
%     for n = -r_roi:r_roi
%          if  m^2 + n^2 <= r_roi^2
%             ROI(m+r_roi+1,n+r_roi+1) = true;
%          end
%     end
% end






CIRCLE = false(2*R2,2*R2);
for m = -R2:R2
    for n = -R2:R2
         if (m^2 + n^2 >= R1^2) && (m^2 + n^2 <= R2^2)
            CIRCLE(m+R2+1,n+R2+1) = true;
         end
    end
end


figure;
hold on;
imagesc([-R2:R2], [-R2:R2],0.2.*double(CIRCLE));
imagesc([-r_roi:r_roi], [-r_roi:r_roi],0.5.*double(ROI));
 
fprintf('Finished with Set Parameters.\n');

%% Calculate neural activity (daily cells)
  


CA_TRACES_RAW = cell(N_days,N_ss);
CA_TRACES_FILT = cell(N_days,N_ss);
CA_TRACES_BIN = cell(N_days,N_ss);
CA_TRACES_ONSET = cell(N_days,N_ss);
CA_F = cell(N_days,N_ss);

CA_TRACES = cell(5,1);

wbar = waitbar(0,'Calculating cell activity for each session');
for curr_day = days
    curr_day_ind = find(days==curr_day);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss); 
        waitbar(((N_ss*(curr_day_ind-1))+curr_ss_ind)/(N_days*N_ss),wbar,'Calculating cell activity for each session')
        
        if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            
            % Load frames for each experiment
            load(strcat(BASEDIR,slsh,'Frames',slsh,'Frames_r_D',num2str(cell_type),'M',num2str(mouse_numb),'_D',num2str(curr_day),'S',num2str(curr_ss),'.mat'));
            % load(strcat(BASEDIR,slsh,'Frames',slsh,'Frames_D',num2str(cell_type),'M',num2str(mouse_numb),'_D',num2str(curr_day),'S',num2str(curr_ss),'.mat'));
            
    
            % Day map
            CURR_CENT_X = cell2mat(CENT_X_DAY(curr_day_ind));
            CURR_CENT_Y = cell2mat(CENT_Y_DAY(curr_day_ind));

    %         % Day 1 map
    %         CURR_CENT_X = cell2mat(CENT_X_DAY(2));
    %         CURR_CENT_Y = cell2mat(CENT_Y_DAY(2));

% %             Single session map
%             CURR_CENT_X = CENT_X{curr_day_ind,curr_ss_ind};
%             CURR_CENT_Y = CENT_Y{curr_day_ind,curr_ss_ind};


            curr_identified_cells = length(CURR_CENT_X);
            
            
             min_frame = MIN_FRAME{curr_day_ind,curr_ss_ind};
             
             
             
%             p = 6;
%             aaa = 5; bbb = 5;
%             F00 = min_frame(max(round(CURR_CENT_Y(p)-floor(aaa/2)),1):min(round(CURR_CENT_Y(p)+ceil(aaa/2)-1),size(frames,1)), max(round(CURR_CENT_X(p)-floor(bbb/2)),1):min(round(CURR_CENT_X(p)+ceil(bbb/2)-1),size(frames,2)));
%             F00 = double(F00);
%             F00(F00<1) = nan;
%             F0(p) = nanmean(nanmean(F00));
%             F0(p)



            % Calculate calcium transients
           
            
            [SPIKES, FF] = cell_activity3(curr_identified_cells, CURR_CENT_X, CURR_CENT_Y, ROI, CIRCLE, frames, min_frame);
            
       

            SPIKES_FILT = SPIKES';
            SPIKES_RAW = SPIKES';
            
            CURR_TRACES_BIN = false(curr_identified_cells,size(SPIKES_RAW,2));
            CURR_TRACES_DER_BIN = false(curr_identified_cells,size(SPIKES_RAW,2));
            CURR_TRACES_ONSET = false(curr_identified_cells,size(SPIKES_RAW,2));
            CURR_TRACES_DER = zeros(size(SPIKES_RAW));
            
                    

    %         % Remove moving average
    %         for k = 1:size(SPIKES_FILT,1)
    %             curr_trace = SPIKES_FILT(k,:);
    %             curr_avg = tsmovavg(curr_trace','s',mov_avg,1);
    %             curr_avg(isnan(curr_avg)) = curr_avg(mov_avg+1);
    %             SPIKES_FILT(k,:) = curr_trace - curr_avg';
    %         end

    
    
            % Remove trace average
            for k = 1:size(SPIKES_FILT,1)
                curr_trace = SPIKES_FILT(k,:);
                
%                 % Trace average:
%                 curr_avg = mean(curr_trace);
                
                % Mode:
                [tmp1 tmp2] = hist(curr_trace,100);
                [tmp3 tmp4] = max(tmp1);
                curr_avg = tmp2(tmp4);

                % SPIKES_FILT(k,:) = curr_trace;
                SPIKES_FILT(k,:) = curr_trace - curr_avg;
            end
            
            % Remove trace average
            for k = 1:size(SPIKES_RAW,1)
                curr_trace = SPIKES_RAW(k,:);
                
%                 % Trace average:
%                 curr_avg = mean(curr_trace);
                
                % Mode:
                [tmp1 tmp2] = hist(curr_trace,100);
                [tmp3 tmp4] = max(tmp1);
                curr_avg = tmp2(tmp4);

                % SPIKES_FILT(k,:) = curr_trace;
                SPIKES_RAW(k,:) = curr_trace - curr_avg;
            end
            

            
            
            % Apply Butterworth filter
            for k = 1:size(SPIKES_FILT,1)
                curr_trace = SPIKES_FILT(k,:);
                SPIKES_FILT(k,:) = filtfilt(butter_b,butter_a,curr_trace); 
            end
   

            % Calculate calcium transient statistics
            [N_spikes, activity, spike_amp, spike_freq, baseline, SPIKES2, SYNC] = neural_stat(SPIKES, bin_size, baseline_thresh, event_thresh);

            % SPIKES2 = SPIKES;


            % Remove time average
            % SPIKES_FILT = SPIKES2';
            
            %         filt_coeff = ones(1, avg)/avg;
            %         for k = 1:curr_identified_cells
            %             curr_trace = SPIKES2(:,k);
            %             curr_trace2 = [zeros(avg,1); curr_trace];
            %             curr_trace_filt2 = filtfilt(filt_coeff, 1, curr_trace2);
            %             curr_trace_filt = curr_trace_filt2(avg+1:end);
            %
            %             SPIKES_FILT(k,:) = curr_trace - curr_trace_filt;
            %         end
            
            % figure; plot(curr_trace); hold on; plot(curr_trace_filt)
            
            
            
            
            
            %         SPIKES_FILT = SPIKES2';
            %         for k = 1:curr_identified_cells
            %             curr_trace = SPIKES2(:,k);
            %             curr_thresh = 3.*rms(curr_trace);
            %             % curr_thresh = 5;
            %             SPIKES_FILT(find(curr_trace<curr_thresh)) = 0;
            %         end


    
    
    
            
            
            
            
            
            
            
            
            
            % %  Calculate spike onsets
            for k = 1:size(SPIKES_FILT,1)
                
                curr_trace = SPIKES_RAW(k,:);
                curr_trace(curr_trace<noise_thresh) = 0;
                SPIKES_FILT2 = SPIKES_FILT;
                SPIKES_FILT2(SPIKES_FILT2<noise_thresh) = 0;
                
                CURR_TRACES_DER(k,2:end) = diff(SPIKES_FILT2(k,:));
                CURR_TRACES_DER(k,1) = CURR_TRACES_DER(k,2);
                
                % Binary derivative traces
                n_outliers = find(abs(curr_trace) < 5);
                curr_baseline = curr_trace(n_outliers);
                curr_rms = rms(curr_baseline);
                % curr_thresh = 2.*rms(curr_trace);
                curr_thresh = bin_thresh*curr_rms; % *** binary threshold
                CURR_TRACES_BIN(k,find(curr_trace>curr_thresh)) = true;
                % SPIKES_FILT(find(curr_trace<curr_thresh)) = 0;
                
                % Trace onset
                curr_rms = rms(CURR_TRACES_DER(k,n_outliers));
                curr_thresh = onset_thresh*curr_rms; % *** binary threshold
                CURR_TRACES_DER_BIN(k,find(CURR_TRACES_DER(k,:)>curr_thresh)) = true;
                % Calculate spike onsets
                
                tmp_p = [CURR_TRACES_DER_BIN(k,2:end) NaN];
                tmp_p2 = find(CURR_TRACES_DER_BIN(k,:) < 0.5 & tmp_p > 0.5);
                % tmp_p2 = tmp_p2+((delay).*ones(length(tmp_p2),1));
                if ~isempty(tmp_p2)
                    if tmp_p2(1) == 0 % could happen
                        tmp_p2(1) = 1;
                    end
                end
                
                % N_spikes(m) = size(tmp2,1);
                %if ~isempty(tmp_p2)
                CURR_TRACES_ONSET(k,tmp_p2)=true;
                %end
            end
            
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%             % Calculate binary data
%             SPIKES_BIN = false(curr_identified_cells,size(SPIKES,1));
%             for k = 1:curr_identified_cells
%                 curr_trace = SPIKES_FILT(k,:);
% 
%                 % curr_thresh = 2.*rms(curr_trace);
%                 curr_thresh = 5; % *** binary threshold
% 
%                 SPIKES_BIN(k,find(curr_trace>curr_thresh)) = true; 
%                 % SPIKES_FILT(find(curr_trace<curr_thresh)) = 0;
%             end
            
            
            for kk = 1:curr_identified_cells
                curr_trace = SPIKES_RAW(kk,:);
                curr_trace(curr_trace<noise_thresh) = 0;
                n_outliers = find(abs(curr_trace) < 5);
                curr_baseline = curr_trace(n_outliers);
                curr_rms = rms(curr_baseline);
                % curr_thresh = 2.*rms(curr_trace);
                curr_thresh = bin_thresh*curr_rms; % *** binary threshold
                CURR_TRACES_BIN(kk,find(curr_trace>curr_thresh)) = true;
                % SPIKES_FILT(find(curr_trace<curr_thresh)) = 0;
            end



            % Filter binary data by minimum duration




%             % Get onset of spikes
%             SPIKES_ONSET = false(size(SPIKES_BIN));
%             for n = 1:curr_identified_cells
%                 tmp = [SPIKES_BIN(n,2:end) NaN];
%                 tmp2 = find(SPIKES_BIN(n,:) < 0.5 & tmp > 0.5);
%                 N_events = size(tmp2,1);
%                 SPIKES_ONSET(n,tmp2)=1;
%             end

            SPIKES_ONSET = CURR_TRACES_ONSET;
            SPIKES_BIN = CURR_TRACES_BIN;



            % Save data

            CA_TRACES_RAW{curr_day_ind,curr_ss_ind} = SPIKES_RAW;
            CA_TRACES_FILT{curr_day_ind,curr_ss_ind} = SPIKES_FILT;
            CA_TRACES_BIN{curr_day_ind,curr_ss_ind} = SPIKES_BIN;
            CA_TRACES_ONSET{curr_day_ind,curr_ss_ind} = SPIKES_ONSET;
            CA_F{curr_day_ind,curr_ss_ind} = FF;

            clear FF
        end
    end
end
close(wbar)





CA_TRACES{1} = CA_TRACES_RAW;
CA_TRACES{2} = CA_TRACES_FILT;
CA_TRACES{3} = CA_TRACES_BIN;
CA_TRACES{4} = CA_TRACES_ONSET;
CA_TRACES{5} = CA_F;


% Save neural data
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{18} = CA_TRACES;'])





fprintf('Finished calculate neural activity (daily cells).\n');

%% Calculate neural activity (global cells)



% Initialize variables
eval(['CA_TRACES = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{18};'])

CA_TRACES_RAW = CA_TRACES{1};
CA_TRACES_FILT = CA_TRACES{2};
CA_TRACES_BIN = CA_TRACES{3};
CA_TRACES_ONSET = CA_TRACES{4};
CA_F = CA_TRACES{5};
for curr_day = days
    curr_day_ind = find(days==curr_day);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        CA_TRACES_RAW_G{curr_day_ind,curr_ss_ind} = zeros(size(IDX_MAP,1),size(CA_TRACES_RAW{curr_day_ind,curr_ss_ind},2));
        CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind} = zeros(size(IDX_MAP,1),size(CA_TRACES_FILT{curr_day_ind,curr_ss_ind},2));
        CA_TRACES_BIN_G{curr_day_ind,curr_ss_ind} = false(size(IDX_MAP,1),size(CA_TRACES_BIN{curr_day_ind,curr_ss_ind},2));
        CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind} = false(size(IDX_MAP,1),size(CA_TRACES_ONSET{curr_day_ind,curr_ss_ind},2));
        CA_F_G{curr_day_ind,curr_ss_ind} = zeros(size(CA_F{curr_day_ind,curr_ss_ind},1),size(IDX_MAP,1));
    end
end


% Assign activity to commonly identified cells
for k = 1:size(IDX_MAP,1)
    for curr_day = days
        curr_day_ind = find(days==curr_day);
        if IDX_MAP(k,curr_day_ind) ~= 0
            for curr_ss = sessions{curr_day_ind}
                curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
                CA_TRACES_RAW_G{curr_day_ind,curr_ss_ind}(k,:) = CA_TRACES_RAW{curr_day_ind,curr_ss_ind}(IDX_MAP(k,curr_day_ind),:);
                CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind}(k,:) = CA_TRACES_FILT{curr_day_ind,curr_ss_ind}(IDX_MAP(k,curr_day_ind),:);
                CA_TRACES_BIN_G{curr_day_ind,curr_ss_ind}(k,:) = CA_TRACES_BIN{curr_day_ind,curr_ss_ind}(IDX_MAP(k,curr_day_ind),:);
                CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind}(k,:) = CA_TRACES_ONSET{curr_day_ind,curr_ss_ind}(IDX_MAP(k,curr_day_ind),:);
                CA_F_G{curr_day_ind,curr_ss_ind}(:,k) = CA_F{curr_day_ind,curr_ss_ind}(:,IDX_MAP(k,curr_day_ind),:);
            end
        end 
    end   
end






% Calculate cell activity for silent cells
wbar = waitbar(0,'Calculating cell activity for global cells');

for curr_day = days
    curr_day_ind = find(days==curr_day);
    
    
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        
        
        if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            
            % Load frames for each experiment
            load(strcat(BASEDIR,slsh,'Frames',slsh,'Frames_r_D',num2str(cell_type),'M',num2str(mouse_numb),'_D',num2str(curr_day),'S',num2str(curr_ss),'.mat'));
            
            
            
            
            for k = 1:size(IDX_MAP,1)
                if (IDX_MAP(k,curr_day_ind)==0)
      
                    waitbar(((N_ss*(curr_day_ind-1))+curr_ss_ind)/(N_days*N_ss),wbar,'Calculating cell activity for global cells')
                    % One cell at a time
                    CURR_CENT_X = CENT_X_GLOBAL(k);
                    CURR_CENT_Y = CENT_Y_GLOBAL(k);
                    curr_identified_cells = length(CURR_CENT_X); % = 1
                    
                    
                    
                    
                  
                    % Calculate calcium transients
                    [SPIKES FF] = cell_activity3(curr_identified_cells, CURR_CENT_X, CURR_CENT_Y, ROI, CIRCLE, frames, MIN_FRAME{curr_day_ind,curr_ss_ind});
                    
                    SPIKES_FILT = SPIKES';
                    SPIKES_RAW = SPIKES';
                    
                    
                    CURR_TRACES_BIN = false(curr_identified_cells,size(SPIKES_RAW,2));
                    CURR_TRACES_DER_BIN = false(curr_identified_cells,size(SPIKES_RAW,2));
                    CURR_TRACES_ONSET = false(curr_identified_cells,size(SPIKES_RAW,2));
                    CURR_TRACES_DER = zeros(size(SPIKES_RAW));
                    
                    
                    
                    
                    %         % Remove moving average
                    %         for k = 1:size(SPIKES_FILT,1)
                    %             curr_trace = SPIKES_FILT(k,:);
                    %             curr_avg = tsmovavg(curr_trace','s',mov_avg,1);
                    %             curr_avg(isnan(curr_avg)) = curr_avg(mov_avg+1);
                    %             SPIKES_FILT(k,:) = curr_trace - curr_avg';
                    %         end
                    
                    
                    
                    % Remove trace average
                    curr_trace = SPIKES_FILT;
                    %                 % Trace average:
                    %                 curr_avg = mean(curr_trace);
                    % Mode:
                    [tmp1 tmp2] = hist(curr_trace,100);
                    [tmp3 tmp4] = max(tmp1);
                    curr_avg = tmp2(tmp4);
                    % SPIKES_FILT(k,:) = curr_trace;
                    SPIKES_FILT = curr_trace - curr_avg;
                    
                    
                    
                    % Remove trace average
                    curr_trace = SPIKES_RAW;
                    %                 % Trace average:
                    %                 curr_avg = mean(curr_trace);
                    % Mode:
                    [tmp1 tmp2] = hist(curr_trace,100);
                    [tmp3 tmp4] = max(tmp1);
                    curr_avg = tmp2(tmp4);
                    % SPIKES_FILT(k,:) = curr_trace;
                    SPIKES_RAW = curr_trace - curr_avg;
                    
                    
                    % Apply Butterworth filter
                    
                    curr_trace = SPIKES_FILT;
                    SPIKES_FILT = filtfilt(butter_b,butter_a,curr_trace);
                    
                    
                    
                    
                    
                    % Calculate calcium transient statistics
                    [N_spikes, activity, spike_amp, spike_freq, baseline, SPIKES2, SYNC] = neural_stat(SPIKES, bin_size, baseline_thresh, event_thresh);
                    
                    % SPIKES2 = SPIKES;
                    
                    
                    
                    
                    
                    
                    
                    % %  Calculate spike onsets
                    for kk = 1:size(SPIKES_FILT,1)
                        
                        curr_trace = SPIKES_RAW(kk,:);
                        curr_trace(curr_trace<noise_thresh) = 0;
                        SPIKES_FILT2 = SPIKES_FILT;
                        SPIKES_FILT2(SPIKES_FILT2<noise_thresh) = 0;
                        
                        CURR_TRACES_DER(kk,2:end) = diff(SPIKES_FILT2(kk,:));
                        CURR_TRACES_DER(kk,1) = CURR_TRACES_DER(kk,2);
                        
                        % Binary derivative traces
                        n_outliers = find(abs(curr_trace) < 5);
                        curr_baseline = curr_trace(n_outliers);
                        curr_rms = rms(curr_baseline);
                        % curr_thresh = 2.*rms(curr_trace);
                        curr_thresh = bin_thresh*curr_rms; % *** binary threshold
                        CURR_TRACES_BIN(kk,find(curr_trace>curr_thresh)) = true;
                        % SPIKES_FILT(find(curr_trace<curr_thresh)) = 0;
                        
                        % Trace onset
                        curr_rms = rms(CURR_TRACES_DER(kk,n_outliers));
                        curr_thresh = onset_thresh*curr_rms; % *** binary threshold
                        CURR_TRACES_DER_BIN(kk,find(CURR_TRACES_DER(kk,:)>curr_thresh)) = true;
                        % Calculate spike onsets
                        
                        tmp_p = [CURR_TRACES_DER_BIN(kk,2:end) NaN];
                        tmp_p2 = find(CURR_TRACES_DER_BIN(kk,:) < 0.5 & tmp_p > 0.5);
                        % tmp_p2 = tmp_p2+((delay).*ones(length(tmp_p2),1));
                        if ~isempty(tmp_p2)
                            if tmp_p2(1) == 0 % could happen
                                tmp_p2(1) = 1;
                            end
                        end
                        
                        % N_spikes(m) = size(tmp2,1);
                        %if ~isempty(tmp_p2)
                        CURR_TRACES_ONSET(kk,tmp_p2)=true;
                        %end
                    end
                    
                    
                    
            
                    
%                     % Calculate binary data
%                     SPIKES_BIN = false(1,size(SPIKES,1));
%                     curr_trace = SPIKES_FILT;
%                     % curr_thresh = 2.*rms(curr_trace);
%                     % curr_thresh = 5;%
%                     SPIKES_BIN(1,find(curr_trace>curr_thresh)) = true;
%                     % SPIKES_FILT(find(curr_trace<curr_thresh)) = 0;
                    
                    
                    
                    
                    
                    % Binary traces
                    for kk = 1:size(SPIKES_FILT,1)
                        curr_trace = SPIKES_RAW(kk,:);
                        curr_trace(curr_trace<noise_thresh) = 0;
                        n_outliers = find(abs(curr_trace) < 5);
                        curr_baseline = curr_trace(n_outliers);
                        curr_rms = rms(curr_baseline);
                        % curr_thresh = 2.*rms(curr_trace);
                        curr_thresh = bin_thresh*curr_rms; % *** binary threshold
                        CURR_TRACES_BIN(kk,find(curr_trace>curr_thresh)) = true;
                        % SPIKES_FILT(find(curr_trace<curr_thresh)) = 0;
                    end
                    
                    
            
            
            
                    
                    
                    % Filter binary data by minimum duration
                    
                    
                    
                    
%                     % Get onset of spikes
%                     SPIKES_ONSET = false(size(SPIKES_BIN));
%                     
%                     tmp = [SPIKES_BIN(1,2:end) NaN];
%                     tmp2 = find(SPIKES_BIN(1,:) < 0.5 & tmp > 0.5);
%                     N_events = size(tmp2,1);
%                     SPIKES_ONSET(1,tmp2)=1;


                    SPIKES_ONSET = CURR_TRACES_ONSET;
                    SPIKES_BIN = CURR_TRACES_BIN;
                    
                    
                    
                    
                    % Save data
                    
                    CA_TRACES_RAW_G{curr_day_ind,curr_ss_ind}(k,:) = SPIKES_RAW;
                    CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind}(k,:) = SPIKES_FILT;
                    CA_TRACES_BIN_G{curr_day_ind,curr_ss_ind}(k,:) = SPIKES_BIN;
                    CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind}(k,:) = SPIKES_ONSET;
                    CA_F_G{curr_day_ind,curr_ss_ind}(:,k) = FF;
                    clear FF
                    
                    max_value = find(SPIKES_FILT>= max_delta_F);
                    if ~isempty(max_value)
                        
                    end
                end
            end
        end
    end
end

close(wbar)



CA_TRACES_G{1} = CA_TRACES_RAW_G;
CA_TRACES_G{2} = CA_TRACES_FILT_G;
CA_TRACES_G{3} = CA_TRACES_BIN_G;
CA_TRACES_G{4} = CA_TRACES_ONSET_G;
CA_TRACES_G{5} = CA_F_G;


% Save neural data
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{21} = CA_TRACES_G;'])




% %% Normalize traces
% 
% % Calculate total trace length
% tot_length = 0;
% tot_cells = size(CA_TRACES_RAW_G{1,1},1);
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         tot_length = tot_length + size(CA_TRACES_RAW{curr_day_ind,curr_ss_ind},2);
%     end
% end
%  
% % Create TOT_TRACES
% TOT_TRACES = zeros(tot_cells, tot_length);
% curr_ind = 1;
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         waitbar(((N_ss*(curr_day_ind-1))+curr_ss_ind)/(N_days*N_ss),wbar,'Normalizing calcium traces')
%         if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})  
%             TOT_TRACES(:,curr_ind:curr_ind + size(CA_TRACES_RAW{curr_day_ind,curr_ss_ind},2)-1) = CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind};
%             curr_ind = curr_ind + size(CA_TRACES_RAW{curr_day_ind,curr_ss_ind},2);
%         end
%     end
% end
% 
% ZERO_PERC = zeros(tot_cells,1);
% ONE_PERC = zeros(tot_cells,1);
% 
% % Calculate normalization parameters for each trace     
% wbar = waitbar(0,'Calculating normalization parameters');
% for k = 1:tot_cells
%     waitbar(k/tot_cells,wbar,'Calculating normalization parameters')
%     curr_trace = TOT_TRACES(k,:); 
%     [tmp1 tmp2] = hist(curr_trace,100);
%     
%     tmp3 = cumsum(tmp1);
%     tmp3 = tmp3/max(tmp3);
%     tmp4 = find(tmp3<= 0.05);
%     ind1 = tmp4(end);
%     tmp5 = find(tmp3<=0.8);
%     ind2 = tmp5(end);
%     ZERO_PERC(k) = tmp2(ind1);
%     ONE_PERC(k) = tmp2(ind2);
%     
%     TOT_TRACES(k,:) = (curr_trace-ZERO_PERC(k))/ONE_PERC(k);
%     
%     % [tmp3 tmp4] = max(tmp1);
%     % curr_avg = tmp2(tmp4);
% end
% close(wbar);
% 
% % Normalize traces
% curr_ind = 1;
% wbar = waitbar(0,'Normalizing calcium traces');
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         waitbar(((N_ss*(curr_day_ind-1))+curr_ss_ind)/(N_days*N_ss),wbar,'Normalizing calcium traces')   
%         if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})  
%             CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind} = TOT_TRACES(:,curr_ind:curr_ind + size(CA_TRACES_RAW{curr_day_ind,curr_ss_ind},2)-1);
%             curr_ind = curr_ind + size(CA_TRACES_RAW{curr_day_ind,curr_ss_ind},2);
%         end
%     end
% end
% close(wbar);
% 
% CA_TRACES_G{2} = CA_TRACES_FILT_G;
% 
% % Save neural data
% eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{21} = CA_TRACES_G;'])
% 
% clear TOT_TRACES;

fprintf('Finished calculate neural activity (global cells).\n');







%% Calculate calcium trace statistics

decay_time = 50;       % Frames after peak time to use for calculation of decay time
gof_thresh = 10;       % goodness of fit threshold, exclude multiple spikes and noisy signals
peak_thresh = 3;        % Minimum trace value to be detected after peak_time frames from onset to keep each time onset
peak_time = 50;         % Number of frames after spike onset to detect spike peak and calculate avg spike peak
use_filt_onset = 1;     % Use filtered onset (remove peaks below peak_thresh), or find all maxima above peak_thresh
min_spikes = 5;         % Minimum number of total gof spikes to calculate amplitude and decay time

% decay_time = 50;       % Frames after peak time to use for calculation of decay time
% gof_thresh = 50;       % goodness of fit threshold, exclude multiple spikes and noisy signals
% peak_thresh = 5;        % Minimum trace value to be detected after peak_time frames from onset to keep each time onset
% peak_time = 50;         % Number of frames after spike onset to detect spike peak and calculate avg spike peak
% use_filt_onset = 1;     % Use filtered onset (remove peaks below peak_thresh), or find all maxima above peak_thresh
% min_spikes = 5;         % Minimum number of total gof spikes to calculate amplitude and decay time




% Load neural data

eval(['CA_TRACES_G = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{21};'])


CA_TRACES_RAW_G = CA_TRACES_G{1};
CA_TRACES_FILT_G = CA_TRACES_G{2};
CA_TRACES_BIN_G = CA_TRACES_G{3};
CA_TRACES_ONSET_G = CA_TRACES_G{4};
CA_F_G = CA_TRACES_G{5};






format long

tot_cells = length(CENT_X_GLOBAL);

excluded_spikes = nan(1,3); % 1 = no peak detected above peak_thresh; 2 = gof threshold not met; 3 = tot_spikes





wbar = waitbar(0,strcat('Processing...'));


tot_cell_progress = 0;


ex_spikes = 0; % Total spikes excluded for gof
tot_spikes = 0; % Total number of spikes
tot_spikes2 = 0; % Total spikes filtered spikes
noisy_cells = [];

tot_spikes_per_cell = zeros(1,tot_cells);


col_ind = 0;

AVG_FREQ = nan(tot_cells,1);
AVG_AMPL = nan(tot_cells,1);
AVG_ACT = nan(tot_cells,1);
AVG_DECAY = nan(tot_cells,1);








% TOTAL DATA

% Create total time vector
tot_time = [];
for curr_day = days
    curr_day_ind = find(days==curr_day);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
        if ~isempty(curr_ss_ind)
            if ~isempty(tot_time)
                tot_time = [tot_time; tot_time(end)+T{curr_day_ind, curr_ss_ind}];
            else
                tot_time = T{curr_day_ind, curr_ss_ind};
            end
        else
            fprintf(strcat('Warning: missing D',num2str(curr_cell_type),'M',num2str(curr_mouse_numb),' day ',num2str(curr_day),' session ',num2str(curr_ss),'\n'));
        end
    end
end

% Merge traces into one
CURR_TRACES_TOT = zeros(tot_cells, length(tot_time));
CURR_TRACES_RAW_TOT = zeros(tot_cells, length(tot_time));
CURR_TRACES_ONSET_TOT = zeros(tot_cells, length(tot_time));
% CURR_TRACES_BIN_TOT = zeros(tot_cells, length(tot_time));
% SPEED_TOT = zeros(1, length(tot_time));
curr_time_ind = 1;
for curr_day = days
    curr_day_ind = find(days==curr_day);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
        if ~isempty(curr_ss_ind)
            CURR_TRACES_ONSET_TOT(:,curr_time_ind:curr_time_ind + length(T{curr_day_ind,curr_ss_ind})-1) = CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind};
            % CURR_TRACES_BIN_TOT(:,curr_time_ind:curr_time_ind + length(T{curr_day_ind,curr_ss_ind})-1) = CA_TRACES_BIN_G{curr_day_ind,curr_ss_ind};
            CURR_TRACES_TOT(:,curr_time_ind:curr_time_ind + length(T{curr_day_ind,curr_ss_ind})-1) = CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind};
            CURR_TRACES_RAW_TOT(:,curr_time_ind:curr_time_ind + length(T{curr_day_ind,curr_ss_ind})-1) = CA_TRACES_RAW_G{curr_day_ind,curr_ss_ind};
            % SPEED_TOT(curr_time_ind:curr_time_ind + length(T{curr_day_ind,curr_ss_ind})-1) = 0.0947 + ((0.6271-0.0947)/300).*T{curr_day_ind,curr_ss_ind};
            curr_time_ind = curr_time_ind + length(T{curr_day_ind,curr_ss_ind});
        end
    end
end




T_s = (tot_time(end)-tot_time(1))/length(tot_time);


% Remove onset with peakes below peak_thresh
for curr_cell = 1:tot_cells
    curr_onset0 = CURR_TRACES_ONSET_TOT(curr_cell,1:end-peak_time);
    tot_spikes = tot_spikes + sum(CURR_TRACES_ONSET_TOT(curr_cell,1:end-peak_time));
    for curr_onset = find(curr_onset0)
        CURR_TRACES_ONSET_TOT(curr_cell,curr_onset) = max(CURR_TRACES_TOT(curr_cell,curr_onset:curr_onset+peak_time)) >= peak_thresh;
    end
    tot_spikes2 = tot_spikes2 + sum(CURR_TRACES_ONSET_TOT(curr_cell,1:end-peak_time));
end


if ~use_filt_onset
    tot_spikes = 0;
    tot_spikes2 = 0;
end


% Calculate avg amplitude/decay
AVG_PEAK = nan(tot_cells,1);
AVG_DEC = nan(tot_cells,1);
AVG_PEAK_STD = nan(tot_cells,1);
AVG_DEC_STD = nan(tot_cells,1);
for curr_cell = 1:tot_cells
    waitbar(tot_cell_progress/tot_cells,wbar,strcat('Processing...'));
    tot_cell_progress = tot_cell_progress + 1;
    fprintf('Tot progress cell %d/%d: ',tot_cell_progress,tot_cells);
    
    
    if use_filt_onset
        % Use filtered trace onset vector
        curr_onset0 = CURR_TRACES_ONSET_TOT(curr_cell,1:end-peak_time);
        
    else
        % Use local maxima
        curr_trace = CURR_TRACES_TOT(curr_cell,1:end-peak_time);  
        % figure; hold on; plot(CURR_TRACES_RAW_TOT(curr_cell,1:end-peak_time)); 
        curr_trace(curr_trace < peak_thresh) = 0;
        [pks,locs] = findpeaks(curr_trace);
        curr_onset0 = false(size(curr_trace));
        curr_onset0(locs) = true;
        tot_spikes = tot_spikes + sum(curr_onset0);
        tot_spikes2 = tot_spikes2 + sum(curr_onset0);
        
        % plot(locs,(max(curr_trace)+1).*ones(1,length(locs)),'v'); % ylim([-3 10]);
    end
    
  
    tot_spikes_per_cell(curr_cell) = sum(curr_onset0);
    
    
    
    avg_peak = zeros(1,length(find(curr_onset0)));
    decays = nan(1, length(find(curr_onset0)));
    
    
    % % %figure; hold on; plot(CURR_TRACES_RAW_TOT(curr_cell,:)); stem(curr_onset0.*max(CURR_TRACES_RAW_TOT(curr_cell,:)));
    
    
    
    % figure; hold on;
    
    curr_ind = 0;
    curr_ind1 = 0;
    curr_onset1 = find(curr_onset0);
    
    for curr_onset = curr_onset1
        
        curr_ind = curr_ind + 1;
        
        [curr_max, curr_max_ind] = max(CURR_TRACES_RAW_TOT(curr_cell,curr_onset:curr_onset+peak_time));
        peak_t = curr_max_ind + curr_onset - 1;
        if (peak_t > 10) && (peak_t < size(CURR_TRACES_RAW_TOT,2)-decay_time)
            curr_time_v = tot_time(peak_t:peak_t + decay_time)-tot_time(peak_t);
            
            warning('off','all');
            tmp = CURR_TRACES_RAW_TOT(curr_cell,peak_t:peak_t+decay_time)';
            tmp(isnan(tmp)) = 0;
            [curr_fit, gof] = fit(curr_time_v, tmp,'exp1');
            warning('on','all');
            
      
            
            if (gof.sse <= gof_thresh)  && (-1/curr_fit.b)>0 && (-1/curr_fit.b)<60
                % plot(curr_time_v,CURR_TRACES_TOT(curr_cell,peak_t:peak_t+decay_time),'color',[0,0,0,0.2],'lineWidth',0.5);
                % plot(curr_time_v,curr_fit.a*exp(curr_fit.b.*curr_time_v),'color',[0,1,0,0.2],'lineWidth',0.5);
                decays(curr_ind) = -1/curr_fit.b;
                curr_ind1 = curr_ind1 + 1;
            else
                % plot(curr_time_v,CURR_TRACES_TOT(curr_cell,peak_t:peak_t+decay_time),'color',[0,0,0,0.2],'lineWidth',0.5);
                % plot(curr_time_v,curr_fit.a*exp(curr_fit.b.*curr_time_v),'color',[1,0,0,0.2],'lineWidth',0.5);
                ex_spikes = ex_spikes + 1;
            end
        end
        avg_peak(curr_ind) = curr_max;
    end
    
    avg_decay= nan;
    avg_peak0 = nan;
    if (curr_ind1 >= min_spikes) && (curr_ind1/sum(curr_onset0) >= 0.2)
        avg_decay = nanmean(decays);
        avg_peak0 = nanmean(avg_peak);
        AVG_PEAK(curr_cell) = avg_peak0;
        AVG_DEC(curr_cell) = avg_decay;
        AVG_PEAK_STD(curr_cell) = nanstd(avg_peak);
        AVG_DEC_STD(curr_cell) = nanstd(decays);
    end
    
    
    if  curr_ind1 >= min_spikes
        fprintf('%d filtered onsets (%3.2f%% gof excluded), ampl %3.2f (%3.2fSD), decay %3.2fs (%3.2fs SD)',sum(curr_onset0),100-100*curr_ind1/length(find(curr_onset0)),avg_peak0,nanstd(avg_peak),avg_decay,nanstd(decays));
    else
        if ~use_filt_onset
            tot_spikes2 = tot_spikes2-sum(curr_onset0);
        end
        fprintf('%d filtered onsets, %d passed gof (not enough for statistics)',sum(curr_onset0), curr_ind1);
    end
    
    if (sum(curr_onset0) >= 10) && (curr_ind1/sum(curr_onset0) < 0.2)
       noisy_cells = [noisy_cells; curr_cell]; 
       fprintf(', gof < 20%%, potential noisy cell');
    end
    fprintf('\n');
end


% Get spike statistics
AVG_FREQ = 60.*nansum(CURR_TRACES_ONSET_TOT,2)/tot_time(end);
AVG_AMPL = AVG_PEAK;
AVG_DECAY = AVG_DEC;

AVG_AMPL_STD = AVG_PEAK_STD;
AVG_DECAY_STD = AVG_DEC_STD;

CURR_TRACES_TOT2 = CURR_TRACES_TOT;
CURR_TRACES_TOT2(CURR_TRACES_TOT2<2) = 0;
AVG_ACT = nanmean(CURR_TRACES_TOT2,2);
AVG_ACT_STD = nanstd(CURR_TRACES_TOT2,[],2);


excluded_spikes(1) = tot_spikes;
excluded_spikes(2) = tot_spikes-tot_spikes2;
excluded_spikes(3) = ex_spikes;
fprintf('D%d tot spikes: %d (excluded for no peak %3.2f%%, excluded for goodness of fit %3.2f%%)\n',cell_type,excluded_spikes(1),100.*excluded_spikes(2)/excluded_spikes(1),100.*excluded_spikes(3)/excluded_spikes(1));



close(wbar);


% Save data

NEURAL_STAT = cell(5,1);
NEURAL_STAT{1} = [decay_time, gof_thresh, excluded_spikes(1), excluded_spikes(2), excluded_spikes(3)];
NEURAL_STAT{2} = AVG_FREQ;
NEURAL_STAT{3} = AVG_AMPL;
NEURAL_STAT{4} = AVG_ACT;
NEURAL_STAT{5} = AVG_DECAY;
NEURAL_STAT{6} = noisy_cells;
NEURAL_STAT{7} = tot_spikes_per_cell;

% Save neural data
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{22} = NEURAL_STAT;'])



%% Plot neural stat results 

selected_cells = [1:tot_cells];
% selected_cells = setdiff([1:tot_cells],noisy_cells');

dec_max = 10;

amp_max = 30;

aa = find(isnan(AVG_AMPL));


figure;

subplot(1,3,1); hold on; histogram(AVG_DECAY(selected_cells),[0:dec_max/100:dec_max],'normalization','probability');
title(strcat('Decay'));
xlabel('Decay time [s]');
ylabel('Number of cells');
ylim([0 .2]);


subplot(1,3,2); hold on; histogram(AVG_AMPL(selected_cells),[0:amp_max/100:amp_max],'normalization','probability');
title(strcat('Amplitude'));
xlabel('Amplitude [\DeltaF/F]');
ylabel('Number of cells');
ylim([0 .2]);


subplot(1,3,3); hold on; histogram(AVG_FREQ(selected_cells),[0:0.1:4.5],'normalization','probability');
title(strcat('Frequency'));
xlabel('Frequency [transients/minute]');
ylabel('Number of cells');
ylim([0 .2]);



figure; hold on;

plot3(AVG_AMPL(selected_cells),AVG_DECAY(selected_cells),AVG_FREQ(selected_cells),'.');
xlim([0 amp_max]);
ylim([0 dec_max]);
zlim([0 1.5]);
xlabel('Amplitude \DeltaF/F');
ylabel('Decay time [s]');
zlabel('Frequency [minute_{-1}]');







%% Normalize neural data based on average peak amplitude

% Load neural data

eval(['CA_TRACES_G = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{21};'])

CA_TRACES_RAW_G = CA_TRACES_G{1};
CA_TRACES_FILT_G = CA_TRACES_G{2};
CA_TRACES_BIN_G = CA_TRACES_G{3};
CA_TRACES_ONSET_G = CA_TRACES_G{4};
CA_F_G = CA_TRACES_G{5};



% Day/session, session data

for curr_day = days
    curr_day_ind = find(days==curr_day);
    
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
        
        % Replace NaNs with average
        AVG_AMPL2 = AVG_AMPL;
        AVG_AMPL2(find(isnan(AVG_AMPL2))) = nanmean(AVG_AMPL2);
        
        
        
        
        CURR_TRACES = CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind};
        CURR_TRACES = CURR_TRACES./repmat(AVG_AMPL2,1,size(CURR_TRACES,2));
        
        
        % CURR_TRACES2 = CA_TRACES_RAW_G{curr_day_ind,curr_ss_ind};
        % CURR_TRACES2 = CURR_TRACES2./repmat(AVG_AMPL2(tot_cell_ind:tot_cell_ind+tot_cells-1),1,size(CURR_TRACES2,2));
        
        
        
        % Save session data
        CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind} = CURR_TRACES;
        % CA_TRACES_RAW_G{curr_day_ind,curr_ss_ind} = CURR_TRACES2;
        
    end
end





% Replace data
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{21}{2} = CA_TRACES_FILT_G;'])







   








%% Calculate average fluorescence field
  


AVG_F = cell(N_days,N_ss);

frame_circ_radius = round(width/2)-2;    
FRAME_MASK = zeros(height,width);
for m = 1:2*(frame_circ_radius)+1
    for n = 1:2*(frame_circ_radius)+1
         if ((m-frame_circ_radius-1)^2 + (n-frame_circ_radius-1)^2 <= (frame_circ_radius+1)^2)
            FRAME_MASK(m+round((height/2)-frame_circ_radius),n+round((width/2)-frame_circ_radius)) = 1;
         end
    end
end
figure; imagesc(FRAME_MASK); colormap(gray); axis equal;

FRAME_MASK(FRAME_MASK==0) = nan;



wbar = waitbar(0,'Calculating avg fluorescence for each session');
for curr_day = days
    curr_day_ind = find(days==curr_day);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        waitbar(((N_ss*(curr_day_ind-1))+curr_ss_ind)/(N_days*N_ss),wbar,'Calculating avg fluorescence for each session')
        
        if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            
            
            
            % Load frames for each experiment
            load(strcat(BASEDIR,slsh,'Frames',slsh,'Frames_r_D',num2str(cell_type),'M',num2str(mouse_numb),'_D',num2str(curr_day),'S',num2str(curr_ss),'.mat'));
            
            frames = double(frames);
            for k = 1:size(frames,3)
                frames(:,:,k) = frames(:,:,k) .*(FRAME_MASK);
            end
            
            
            curr_trace = squeeze(nanmean(nanmean(frames,1),2));
            
            % Save data
            
            AVG_F{curr_day_ind,curr_ss_ind} = curr_trace;
            
            clear FF
        end
    end
end
close(wbar)


CA_TRACES_G{6} = AVG_F;


% Save neural data
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{21}{6} = CA_TRACES_G{6};'])



fprintf('Finished calculate average fluoresence field.\n');














%% Plot calcium traces

normalization = 0; % Normalize traces (for visualization only)


 curr_fig_cells = 0;

cell_step = 50;
curr_fig_cells = curr_fig_cells + 1;



% selected_days = days;
% selected_sessions = sessions;

% selected_days = 3;
% selected_sessions = sessions;

selected_days = 6;
selected_sessions = cell(N_days,1);
% selected_sessions{find(days==selected_days)} = sessions{selected_days}%(1:5);
% selected_days = 1;
% selected_sessions = cell(N_days,1);
selected_sessions{find(days==selected_days)} = [1:3];


plot_binary = 0;

tmp_step = 1.5;

tot_cells = length(CENT_X_GLOBAL);



% selected_cells = [13 14 15 16 17 18 19 23 27 28 30 35 36 39 44 46 62 66 83 85 109 127 134 138 154];
% selected_cells = (1:tot_cells);
% selected_cells = [rows' cols'];

selected_cells = [cell_step*(curr_fig_cells-1)+1:min(tot_cells, cell_step*curr_fig_cells)];


% sessions{1} = (1:10);
% sessions{1} = (1:11);



% % Rasterplot
% cmin = 0;
% cmax = 5;
% 
% figure; imagesc(0.1.*[1:size(CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind},2)],[1:size(CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind},1)],CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind}); colormap(gray); caxis([cmin cmax])
% xlabel('Time [s]');
% ylabel('Cell #');
% cbar = colorbar('EastOutside');
% 
% figure; imagesc(0.1.*[1:size(CA_TRACES_BIN_G{curr_day_ind,curr_ss_ind},2)],[1:size(CA_TRACES_BIN_G{curr_day_ind,curr_ss_ind},1)],CA_TRACES_BIN_G{curr_day_ind,curr_ss_ind}); colormap(gray); caxis([0 1.5])
% xlabel('Time [s]');
% ylabel('Cell #');




% Create total time vector
tot_time = [];
for curr_day = selected_days
    curr_day_ind = find(days==curr_day);
    for curr_ss = selected_sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
        if ~isempty(tot_time)
            tot_time = [tot_time; tot_time(end)+T{curr_day_ind, curr_ss_ind}];
        else
            tot_time = T{curr_day_ind, curr_ss_ind};
        end
    end
end



% Plot continuous traces
figure; hold on;
curr_time_ind = 1;
for curr_day = selected_days
    curr_day_ind = find(days==curr_day);
    for curr_ss = selected_sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
        
        % % Choose which trace to plot
      
        curr_trace_filt = CA_TRACES_FILT_G{curr_day_ind, curr_ss_ind};
        
        curr_trace_filt(curr_trace_filt<0) = 0;
        
        curr_trace_raw = CA_TRACES_RAW_G{curr_day_ind, curr_ss_ind};
        % curr_trace_filt(curr_trace<3)=0;
        % curr_trace_raw(curr_trace<3)=0;
        
        
        tmp = length(selected_cells)*tmp_step+tmp_step;
        for kk = 1:length(selected_cells)
            color_ind = circshift([1:10],[10, -kk+1]);
            if normalization == 1
                
                % curr_plot = plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp + curr_trace_raw(selected_cells(kk),:)./max(curr_trace_raw(selected_cells(kk),:)), 'color', COLORS(color_ind(1),:));
                plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp + curr_trace_filt(selected_cells(kk),:)./curr_trace_filt(selected_cells(kk),:), 'color', curr_plot.Color);
                
                
                % curr_plot2 = plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp.*ones(size(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1))),'--k');
            
            curr_onset = nan(1,size(CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind},2));
            curr_onset(find(CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind}(selected_cells(kk),:))) = 1;
            plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp + curr_onset,'v','MarkerSize',5,'markerfacecolor',curr_plot.Color,'markeredgecolor',curr_plot.Color);
            
            
            plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp.*ones(size(CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind}(selected_cells(kk),:))),'--k');
            plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),(tmp+5).*ones(size(CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind}(selected_cells(kk),:))),':k');
    
    
            % text(-length(tot_time)*0.001,tmp+0.5,num2str(selected_cells(kk)),'color',curr_plot.Color);
            text(-length(tot_time)*0.001,tmp+0.5,num2str(kk),'color',curr_plot.Color);
            
            tmp = tmp - tmp_step;
            
            else
                
                
                % curr_col = COLORS(color_ind(1),:);
                
                
                
                % Red/black color for noisy/low activity cells
                
                curr_col = [0 0 1];
                if any(selected_cells(kk)==noisy_cells)
                    curr_col = [1 0 0];
                end
              
                if any(selected_cells(kk)==find(tot_spikes_per_cell<8))
                    curr_col = [0 0 0];
                end
                
                
                
                % curr_plot = plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp + curr_trace_raw(selected_cells(kk),:), 'color', COLORS(color_ind(1),:));
                curr_plot = plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp + curr_trace_filt(selected_cells(kk),:), 'color', curr_col);
                
                
                % curr_plot2 = plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp.*ones(size(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1))),'--k');
                curr_onset = nan(1,size(CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind},2));
                curr_onset(find(CA_TRACES_ONSET_G{curr_day_ind,curr_ss_ind}(selected_cells(kk),:))) = tmp_step-0.2.*tmp_step;
                plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp + curr_onset,'v','MarkerSize',5,'markerfacecolor',curr_plot.Color,'markeredgecolor',curr_plot.Color);
                
                
                % plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp.*ones(size(CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind}(selected_cells(kk),:))),'--k');
                % plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),(tmp+5).*ones(size(CA_TRACES_FILT_G{curr_day_ind,curr_ss_ind}(selected_cells(kk),:))),':k');
                
                
                text(-length(tot_time)*0.001,tmp+tmp_step/6,num2str(selected_cells(kk)),'color',curr_plot.Color);
                
                tmp = tmp - tmp_step;
                
            end
            
            
        end
        
        
        
        curr_time_ind = curr_time_ind + length(T{curr_day_ind, curr_ss_ind});
    end
end 

xlim([-round(0.02*tot_time(end)) tot_time(end)]);
ylim([-tmp_step length(selected_cells)*tmp_step+8*tmp_step]);
xlabel('T [s]');


% Plot day/session separation line
curr_time_ind = 1;
for curr_day = selected_days
    curr_day_ind = find(days==curr_day);
    for curr_ss = selected_sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
        if curr_ss_ind == 1
            line([tot_time(curr_time_ind), tot_time(curr_time_ind)],get(gca,'YLim'),'Color','k','LineStyle','--');
        else
            line([tot_time(curr_time_ind), tot_time(curr_time_ind)],get(gca,'YLim'),'Color','k','LineStyle',':');
        end
        curr_time_ind = curr_time_ind + length(T{curr_day_ind, curr_ss_ind});
    end
end
        



if plot_binary
    % Plot binary traces
    bin_fig = figure;
    curr_time_ind = 1;
    for curr_day = selected_days
        curr_day_ind = find(days==curr_day);
        for curr_ss = selected_sessions{curr_day_ind}
            curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
            
            % Choose which trace to plot
            curr_trace = CA_TRACES_BIN_G{curr_day_ind, curr_ss_ind};
            
            POS_EDGES = false(length(selected_cells),size(curr_trace,2));
            NEG_EDGES = false(length(selected_cells),size(curr_trace,2));
            for m = 1:length(selected_cells)
                tmp = [curr_trace(selected_cells(m),2:end) NaN];
                
                tmp2 = find(curr_trace(selected_cells(m),:) < 0.5 & tmp > 0.5);
                tmp3 = find(curr_trace(selected_cells(m),:) > 0.5 & tmp < 0.5);
                POS_EDGES(m,tmp2)=true;
                NEG_EDGES(m,tmp3)=true;
                if curr_trace(selected_cells(m),1) == 1
                    POS_EDGES(m,1) = 1;
                end
            end
            
            
            
            
            tmp = length(selected_cells)+1;
            for kk = 1:length(selected_cells)
                color_ind = circshift((1:10),[10, -kk+1]);
                
                pe = find(POS_EDGES(kk,:));
                ne = find(NEG_EDGES(kk,:));
                
                
                if length(pe)~=length(ne)
                    
                    ne = [ne size(curr_trace,2)];
                    % ne = [ne size(curr_trace,2)];
                end
                
                for k = 1:length(pe)
                    curr_pol_x = [tot_time(curr_time_ind+pe(k)), tot_time(curr_time_ind+ne(k)-1), tot_time(curr_time_ind+ne(k)-1), tot_time(curr_time_ind+pe(k))];
                    curr_pol_y = [tmp, tmp, tmp+.8, tmp+.8];
                    curr_patch = patch(curr_pol_x,curr_pol_y,COLORS(1,:),'LineStyle','-','EdgeColor',COLORS(1,:));
                    % curr_line = line([tot_time(curr_time_ind+pe(k)), tot_time(curr_time_ind+ne(k)-1)],tmp.*[1 1],'LineWidth',20,'color','k');
                    
                end
                
                % curr_plot2 = plot(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1),tmp.*ones(size(tot_time(curr_time_ind:curr_time_ind + length(T{curr_day_ind, curr_ss_ind})-1))),'--k');
                
                
                text(-round(tot_time(end)*0.08),tmp+0.5,num2str(selected_cells(kk)),'color','k');
                
                tmp = tmp - 1;
                
                
            end
            
            curr_time_ind = curr_time_ind + length(T{curr_day_ind, curr_ss_ind});
        end
    end
    
    xlim([-round(0.02*tot_time(end)) tot_time(end)]);
    ylim([-.5 length(selected_cells)+.5]);
    xlabel('T [s]');
    
    % Plot day/session separation line
    curr_time_ind = 1;
    for curr_day = selected_days
        curr_day_ind = find(days==curr_day);
        for curr_ss = selected_sessions{curr_day_ind}
            curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
            if curr_ss_ind == 1
                line([tot_time(curr_time_ind), tot_time(curr_time_ind)],get(gca,'YLim'),'Color','k','LineStyle','--');
            else
                line([tot_time(curr_time_ind), tot_time(curr_time_ind)],get(gca,'YLim'),'Color','k','LineStyle',':');
            end
            curr_time_ind = curr_time_ind + length(T{curr_day_ind, curr_ss_ind});
        end
    end
    
end
















%% Clear unused variables
all_var = who;
for k = 1:length(all_var)
    if isempty(strfind('all_var slsh cell_type mouse_numb eartag BASEDIR',all_var{k})) && ~strcmp('MOUSE_',all_var{k}(1:min(end,6)))
        eval(['clear ',all_var{k}]);
    end
end

clear all_var k 







%% Save full workspace
fprintf('Saving workspace to disk...\n');
eval(['save(''D:' slsh,'Giovanni',slsh,'data',slsh,'D',num2str(cell_type),'M',num2str(mouse_numb),'_rfull.mat'',''-v7.3'');'])
fprintf('Done saving process_register workspace to disk.\n');














