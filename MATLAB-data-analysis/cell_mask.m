function [CENT_X_COMMON, CENT_Y_COMMON, IDX_MAP] = cell_mask(dist_threshold, MOUSE)
%CELL_MASK Generate single cell mask from different days. Activity for
%cells which are not active for some days needs to be recalculated
%
%   giobarbera@neuralmapper.com


[FRAMES, N_days, N_ss, T, width, height, days, sessions, cell_mask_size, cell_size, min_cell_distance, rms_fact, area_min, T_cell, N_det, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, X_SHIFT, Y_SHIFT, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CENT_X_DAY, CENT_Y_DAY, CA_TRACES_RAW, CA_TRACES_FILT, CA_TRACES_BIN, CA_TRACES_ONSET clust_idx clust_idx_s clust_cell_s  CLUST_CONN_FULL clust_conn_thresh clust_iterations clust_bin_size clust_curr_days clust_curr_sessions clust_metric clust_tot_kmeans_iter clust_max_iter global_clusters] = load_mouse_full(MOUSE);
colors;




MAP = [];


CURR_X = CENT_X_DAY;
CURR_Y = CENT_Y_DAY;



for curr_day1 = days
    curr_day_ind1 = find(days==curr_day1);
   
    curr_x1 = CURR_X{curr_day_ind1};
    curr_y1 = CURR_Y{curr_day_ind1};
    
    matches = zeros(1,length(curr_x1));

    
    for curr_day2 = days
        curr_day_ind2 = find(days==curr_day2);
        if (curr_day2 ~= curr_day1)
            curr_x2 = CURR_X{curr_day_ind2};
            curr_y2 = CURR_Y{curr_day_ind2};
            
            for m = 1:length(curr_x1)
                min_dist = 10000;
                closest_idx = [];
                for n = 1:length(curr_x2)
                    curr_dist = norm([curr_x1(m); curr_y1(m)]-[curr_x2(n); curr_y2(n)],2);
                    if curr_dist <= dist_threshold && curr_dist < min_dist
                        closest_idx = n;
                        min_dist = curr_dist;
                    end
                end
                
                
                
                % Remove candidate cells that are closer to another cell of the
                % reference day
                if ~isempty(closest_idx)
                    min_dist = norm([curr_x1(m); curr_y1(m)]-[curr_x2(closest_idx); curr_y2(closest_idx)],2);
                    for mm = 1:length(curr_x1)
                        curr_dist = norm([curr_x1(mm); curr_y1(mm)]-[curr_x2(closest_idx); curr_y2(closest_idx)],2);
                        if curr_dist < min_dist
                            closest_idx = [];
                            break;
                        end
                    end
                end
                
                
                
                
                if ~isempty(closest_idx)
                    matches(m) = matches(m)+1;
                    
                    
                    tmp2 = zeros(1,N_days);
                    tmp2(curr_day_ind1) = m; 
                    tmp2(curr_day_ind2) = closest_idx;
                    MAP = [MAP; tmp2];
                end
                   
                
                
                % CURR_X{curr_day_ind2}(closest_idx) = [];
                % CURR_Y{curr_day_ind2}(closest_idx) = [];
     
            end
            
            
        
            
        end
    end
    
    isolated = find(matches==0);
    for nn = 1:length(isolated)
        tmp2 = zeros(1,N_days);
        tmp2(curr_day_ind1) = isolated(nn); 
        MAP = [MAP; tmp2];
    end
    
end





% Find common cells indices
IDX_MAP = zeros(1,N_days);
for k = 1:size(MAP,1)
    tmp = find(MAP(k,:));
    
    

    tmp1 = find(IDX_MAP(:,tmp(1))==MAP(k,tmp(1)));
    if length(tmp)>1 % Non isolated cell
        tmp2 = find(IDX_MAP(:,tmp(2))==MAP(k,tmp(2)));
    else 
        tmp2 = [];
    end
    
    if ~isempty(tmp1)
        IDX_MAP(tmp1,tmp(1)) = MAP(k,tmp(1));
        IDX_MAP(tmp1,tmp(2)) = MAP(k,tmp(2));
    elseif ~isempty(tmp2)
        IDX_MAP(tmp1,tmp(1)) = MAP(k,tmp(1));
        IDX_MAP(tmp1,tmp(2)) = MAP(k,tmp(2));
    else
        IDX_MAP = [IDX_MAP; MAP(k,:)];
    end
        
end
IDX_MAP(1,:) = [];



% Find cell centroid for paired maps common to all days (COMMON) and
% partial days (PARTIAL)

CENT_X_COMMON = zeros(size(IDX_MAP,1),1);
CENT_Y_COMMON = zeros(size(IDX_MAP,1),1);

IDX2 = IDX_MAP;
IDX2(IDX2==0) = nan;



for k = 1:size(IDX_MAP,1)
    tmpx = 0;
    tmpy = 0;
    cnt = 0;
    for curr_day = days
        curr_day_ind = find(days == curr_day);
        if IDX_MAP(k,curr_day_ind)~=0
            tmpx = tmpx + CENT_X_DAY{curr_day_ind}(IDX_MAP(k,curr_day_ind));
            tmpy = tmpy + CENT_Y_DAY{curr_day_ind}(IDX_MAP(k,curr_day_ind));
            cnt = cnt + 1;
        end
    end
    CENT_X_COMMON(k) = tmpx/cnt;
    CENT_Y_COMMON(k) = tmpy/cnt;
end
            
            
        



% figure; plot(CENT_X_DAY{1},CENT_Y_DAY{1},'or'); hold on; plot(CENT_X_DAY{2},CENT_Y_DAY{2},'og'); plot(CENT_X_DAY{3},CENT_Y_DAY{3},'ob'); plot(CENT_X_DAY{4},CENT_Y_DAY{4},'oc'); plot(CENT_X_DAY{5},CENT_Y_DAY{5},'om');
% 
% curr_cell = 1;
% for d = 1:5
%     if(IDX_MAP(curr_cell,d) ~= 0)
%         plot(CENT_X_DAY{d}(IDX_MAP(curr_cell,d)),CENT_Y_DAY{d}(IDX_MAP(curr_cell,d)),'xk');
%     end
% end
% 
% plot(CENT_X_COMMON(curr_cell),CENT_Y_COMMON(curr_cell),'dr');








end