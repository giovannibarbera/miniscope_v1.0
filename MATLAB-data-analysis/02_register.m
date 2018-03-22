% close all; 
% clear; 
 clc

%% Load single mouse data

BASEDIR             = '..\data'; 

cell_type = 1;
mouse_numb = 3;

%workspace_name = 'D736_042716';  % ***
% workspace_name = 'D736';

eval(['MOUSE = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL;'])
[FRAMES, N_days, N_ss, T, width, height, days, sessions, cell_mask_size, cell_size, min_cell_distance, rms_fact, area_min, T_cell, N_det, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, X_SHIFT, Y_SHIFT, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CENT_X_DAY, CENT_Y_DAY, CENT_X_GLOBAL, CENT_Y_GLOBAL, IDX_MAP, CA_TRACES_RAW, CA_TRACES_FILT, CA_TRACES_BIN, CA_TRACES_ONSET, CA_TRACES_RAW_G, CA_TRACES_FILT_G, CA_TRACES_BIN_G, CA_TRACES_ONSET_G, clust_idx, clust_idx_s, sort_idx,  CLUST_CONN_FULL, clust_conn_thresh, clust_iterations, clust_bin_size, clust_curr_days, clust_curr_sessions, clust_metric, clust_tot_kmeans_iter, clust_max_iter, global_clusters, eartag, ref_day, ref_angle] = load_mouse_full(MOUSE);

colors;
fprintf('Loaded! Ready to go.\n');





%% Manually split overlapping cells for each sub-session

% % 500um GRIN lens
% min_cell_size = 4;          
% large_cell_thresh = 40;     


% 1000um GRIN lens
min_cell_size = 3;          % Typically 3
large_cell_thresh = 40;     % Typically 16


sub_size = 30; %cell masks are created for all cells active in 30 frames

CENT_X_DAY = cell(N_days,1);
CENT_Y_DAY = cell(N_days,1);

CENT_XY_DAY = cell(4,1); %  CENT_X_DAY, CENT_Y_DAY, CENT_X_DAY_R, CENT_Y_DAY_R



% figure; hold on;
for curr_day = days
    curr_rgb = zeros(height, width, 3);
    curr_map = zeros(height, width);
    curr_day_ind = find(days==curr_day);
    wbar1 = waitbar(0,'Unpacking binary data');
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        
        curr_cells_dyn = logical(zeros(height,width,size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)));
  
        waitbar(curr_ss_ind/N_ss,wbar1);
        
        for jj = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
            
            tmp = bwunpack(CELLS_DYN{curr_day_ind,curr_ss_ind}(:,:,jj));
            curr_cells_dyn(:,:,jj) = tmp(1:height,1:width);
        end
        
        
        if ~isempty(curr_cells_dyn)
            for kk = 1:floor(size(curr_cells_dyn,3)/sub_size)
                curr_col = [1 1 1];%COLORS(curr_ss_ind,:);
                % curr_rgb(:,:,1) = curr_rgb(:,:,1)+1*curr_col(1).*CELLS_FILT{curr_day_ind,curr_ss_ind};
                % curr_rgb(:,:,2) = curr_rgb(:,:,2)+1*curr_col(2).*CELLS_FILT{curr_day_ind,curr_ss_ind};
                % curr_rgb(:,:,3) = curr_rgb(:,:,3)+1*curr_col(3).*CELLS_FILT{curr_day_ind,curr_ss_ind};     
                curr_map = curr_map+(1/(length(sessions{curr_day_ind})*ceil(size(curr_cells_dyn,3)/sub_size))).*double(any(curr_cells_dyn(:,:,(kk-1)*sub_size+1:min(end,kk*sub_size)),3));
            end
        end
    end
    close(wbar1);
    
    curr_rgb(curr_rgb > 1) = 1;
    curr_map(curr_map > 1) = 1;
    
    
    curr_map0 = curr_map;
    % curr_rgb(curr_rgb<1)=0;
    curr_map0(curr_map0>0)=1;
    % curr_rgb(curr_rgb>=0.5)=1; curr_rgb(curr_rgb<0.5)=0;
    % subplot(1, N_days, curr_day_ind);
    
    
    % Remove cells with smaller than min_cell_size (typically 3)
    curr_map2 = bwareaopen(curr_map0,min_cell_size);
    
 
    curr_map = curr_map.*double(curr_map2); 

    
    % Identify all cells
    CELL_MASK = [0 1 1 0;
             1 1 1 1;
             1 1 1 1;
             0 1 1 0];  
    area_min = 3;
    min_cell_distance = 5;
     
    [curr_identified_cells, CURR_CELL_AREAS, CURR_CENT_X, CURR_CENT_Y, CURR_CELL_REGIONS, CURR_CELLS_FILT] = find_cells(logical(curr_map), CELL_MASK, area_min, min_cell_distance);
    CURR_CELL_AREAS = struct2array(CURR_CELL_AREAS);
  
    
    % Augmented x and y location including overlapping cells
    NEW_CENT_X = CURR_CENT_X;
    NEW_CENT_Y = CURR_CENT_Y;
    
    
    % Identify large cell areas which could include more than 1 cell
    large_cells = find(CURR_CELL_AREAS>large_cell_thresh);
    
    
    
    
    % Apply watershed to large cell areas to identify overlapping cells
    
    if length(large_cells) > 0
        fprintf('Found %d large cells on day %d.\n',length(large_cells),curr_day);
        for k = 1:length(large_cells)
            curr_region = CURR_CELL_REGIONS{large_cells(k),1};
            curr_area = curr_map(min(curr_region(:,2)):max(curr_region(:,2)),min(curr_region(:,1)):max(curr_region(:,1)));

            % Find identified centroids within the large cell region
            idx_x = find( CURR_CENT_X>min(curr_region(:,1)) & CURR_CENT_X<max(curr_region(:,1)) & CURR_CENT_Y>min(curr_region(:,2)) & CURR_CENT_Y<max(curr_region(:,2)));
            idx_y = find( CURR_CENT_X>min(curr_region(:,1)) & CURR_CENT_X<max(curr_region(:,1)) & CURR_CENT_Y>min(curr_region(:,2)) & CURR_CENT_Y<max(curr_region(:,2)));
            curr_x = CURR_CENT_X(idx_x) - min(curr_region(:,1))+1;
            curr_y = CURR_CENT_Y(idx_y) - min(curr_region(:,2))+1;
            
            winsize = [ 90 90 30*(max(curr_region(:,1))-min(curr_region(:,1))) 30*(max(curr_region(:,2))-min(curr_region(:,2)))];
            tmp1 = figure('Position',winsize); 
            % imagesc(curr_area); caxis([0 length(sessions{curr_day_ind})]);  hold on;
            imagesc(curr_area); 
            caxis([0 0.3.*max(max(curr_map))]); 
            % colormap(gray);
            hold on;
            old_cells = plot(curr_x,curr_y, '.r','MarkerSize',50);
            
            title(strcat(num2str(k),'/',num2str(length(large_cells))));
            
            
            
            
            
            
            while(1)
                [x_pos,y_pos,button] = ginput(1);
                
                if isempty(x_pos)
                    break;
                end
              

                switch button
                    case 1 % Add new cells
                        new_x = x_pos + min(curr_region(:,1))-1;
                        new_y = y_pos + min(curr_region(:,2))-1;
                        NEW_CENT_X = [NEW_CENT_X new_x'];
                        NEW_CENT_Y = [NEW_CENT_Y new_y'];
                        plot(x_pos,y_pos,'.m','markersize',70);
                    case 3 % Remove cells
                        del_x = x_pos + min(curr_region(:,1))-1;
                        del_y = y_pos + min(curr_region(:,2))-1;
                        plot(x_pos,y_pos,'.k','markersize',70);
                        r_idx = find(CURR_CENT_X>del_x-0.5 & CURR_CENT_X<del_x+0.5 & CURR_CENT_Y >del_y-0.5 & CURR_CENT_Y <del_y+0.5);
                        NEW_CENT_X(r_idx) = nan;
                        NEW_CENT_Y(r_idx) = nan;
                        
                        % old_cells.Marker = 'none';
                end
            end
            

    


            
            
            
            
            
%             old_cells.Marker = 'none';
%             refresh(tmp1);
%             plot(new_x, new_y, '.m','MarkerSize',50);

         
            
            
            
            
            
            
            
            % Old code (no refresh):
            
%             [x_pos,y_pos,button] = ginput;
%             
%             new_x = x_pos(find(button==1)) + min(curr_region(:,1))-1;
%             new_y = y_pos(find(button==1)) + min(curr_region(:,2))-1;
%             del_x = x_pos(find(button==3)) + min(curr_region(:,1))-1;
%             del_y = y_pos(find(button==3)) + min(curr_region(:,2))-1;
%             
% %             old_cells.Marker = 'none';
% %             refresh(tmp1);
% %             plot(new_x, new_y, '.m','MarkerSize',50);
% 
%             % Remove closest cell to the selected points (+- 0.5 pixel)
%             for m = 1:length(del_x)
%                 r_idx = find(CURR_CENT_X>del_x(m)-0.5 & CURR_CENT_X<del_x(m)+0.5 & CURR_CENT_Y >del_y(m)-0.5 & CURR_CENT_Y <del_y(m)+0.5);
%                 NEW_CENT_X(r_idx) = nan;
%                 NEW_CENT_Y(r_idx) = nan;
%             end   
%             % Add new cells
%             NEW_CENT_X = [NEW_CENT_X new_x'];
%             NEW_CENT_Y = [NEW_CENT_Y new_y'];
            
            
            
            
            
            close(tmp1);
          
        end
  
        NEW_CENT_X(isnan(NEW_CENT_X))=[];
        NEW_CENT_Y(isnan(NEW_CENT_Y))=[];

    
    end
    
    
    
    CENT_X_DAY{curr_day_ind} = NEW_CENT_X;
    CENT_Y_DAY{curr_day_ind} = NEW_CENT_Y;
   

%     figure; hold on;
%     imagesc(curr_map); colormap(gray); xlim([0, width]); ylim([0, height]);
%     plot(CURR_CENT_X,CURR_CENT_Y,'or')
%     set(gca, 'YDir', 'reverse');
%     axis equal
%     xlim([0, MOUSE{1}{5}(1)]);
%     ylim([0, MOUSE{1}{5}(1)]);
%     
    figure; hold on;
    imagesc(curr_map); colormap(gray); xlim([0, width]); ylim([0, height]);
    plot(NEW_CENT_X,NEW_CENT_Y,'og')
    % plot(CURR_CENT_X,CURR_CENT_Y,'or')
    set(gca, 'YDir', 'reverse');
    axis equal
    xlim([0, MOUSE{1}{5}(1)]);
    ylim([0, MOUSE{1}{5}(1)]);
    caxis([0 0.005]);
%     
    figure; hold on;
    imagesc(curr_map); colormap(gray); xlim([0, width]); ylim([0, height]);
    set(gca, 'YDir', 'reverse');
    axis equal
    xlim([0, MOUSE{1}{5}(1)]);
    ylim([0, MOUSE{1}{5}(1)]);
    caxis([0 0.005]);
    
%      figure; hold on;
%     imagesc(curr_map0); colormap(gray); xlim([0, width]); ylim([0, height]);
%     set(gca, 'YDir', 'reverse');
%     axis equal
%     xlim([0, MOUSE{1}{5}(1)]);
%     ylim([0, MOUSE{1}{5}(1)]);
    
end

% Save data
CENT_XY_DAY{1} = CENT_X_DAY;
CENT_XY_DAY{2} = CENT_Y_DAY;
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{17} = CENT_XY_DAY;'])















fprintf('You finished splitting cells.\n');






%% Merge overlapping cells



for curr_day = days
    curr_day_ind = find(days==curr_day);
    CURR_CENT_X = CENT_X_DAY{curr_day_ind};
    CURR_CENT_Y = CENT_Y_DAY{curr_day_ind};
    
    
    curr_rgb = zeros(height, width, 3);
    curr_map = zeros(height, width);
    for curr_ss = sessions{curr_day_ind}  
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            curr_col = [1 1 1];%COLORS(curr_ss_ind,:);
            % curr_rgb(:,:,1) = curr_rgb(:,:,1)+1*curr_col(1).*CELLS_FILT{curr_day_ind,curr_ss};
            % curr_rgb(:,:,2) = curr_rgb(:,:,2)+1*curr_col(2).*CELLS_FILT{curr_day_ind,curr_ss};
            % curr_rgb(:,:,3) = curr_rgb(:,:,3)+1*curr_col(3).*CELLS_FILT{curr_day_ind,curr_ss};     
            curr_map = curr_map+(1/N_ss).*CELLS{curr_day_ind,curr_ss_ind};
        end
    end
    curr_rgb(curr_rgb > 1) = 1;
    curr_map(curr_map > 1) = 1;
    
    
    curr_map0 = curr_map;
    % curr_rgb(curr_rgb<1)=0;
    curr_map0(curr_map0>0)=1;
    % curr_rgb(curr_rgb>=0.5)=1; curr_rgb(curr_rgb<0.5)=0;
    % subplot(1, N_days, curr_day_ind);
    
    
    % % Remove cells with smaller than min_cell_size (typically 3)
    %  curr_map2 = bwareaopen(curr_map0,min_cell_size);
    % curr_map = curr_map.*double(curr_map2); 
    
    
    curr_map = ones(size(curr_map))-curr_map;


    
    
    NEW_CENT_X = CURR_CENT_X;
    NEW_CENT_Y = CURR_CENT_Y;
    

    
    
   
    winsize = [ 0 0 1500 1500];
    figure('Position',winsize);        
    hold on;
    imagesc(curr_map); colormap(gray); xlim([0, width]); ylim([0, height]);
    plot(NEW_CENT_X,NEW_CENT_Y,'or')
    % plot(CURR_CENT_X,CURR_CENT_Y,'or')
    set(gca, 'YDir', 'reverse');
    axis equal
    xlim([0, MOUSE{1}{5}(1)]);
    ylim([0, MOUSE{1}{5}(1)]);
            
    
    
    
    
    
    for k = 1:length(NEW_CENT_X)-1
        for kk = k + 1 : length(NEW_CENT_X)
            curr_dist = norm([NEW_CENT_X(k); NEW_CENT_Y(k)] - [NEW_CENT_X(kk); NEW_CENT_Y(kk)],2);
            if (curr_dist <= 3)
                plot(NEW_CENT_X(k),NEW_CENT_Y(k),'ok')
                plot(NEW_CENT_X(kk),NEW_CENT_Y(kk),'ok')
                NEW_CENT_X(k) = (NEW_CENT_X(k)+NEW_CENT_X(kk))/2;
                NEW_CENT_Y(k) = (NEW_CENT_Y(k)+NEW_CENT_Y(kk))/2;
                NEW_CENT_X(kk) = nan;
                NEW_CENT_Y(kk) = nan;
                plot(NEW_CENT_X(k),NEW_CENT_Y(k),'og')
            end
        end
    end
    

    


   
    NEW_CENT_X(isnan(NEW_CENT_X))=[];
    NEW_CENT_Y(isnan(NEW_CENT_Y))=[];

        
  
    
    
    
    
    
    

    
    
    
    CENT_X_DAY{curr_day_ind} = NEW_CENT_X;
    CENT_Y_DAY{curr_day_ind} = NEW_CENT_Y;
   

end

% Save data
CENT_XY_DAY{1} = CENT_X_DAY;
CENT_XY_DAY{2} = CENT_Y_DAY;
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{17} = CENT_XY_DAY;'])






fprintf('You finished merging cells.\n');










%%  Calculate optimal image shift/rotation

ref_day_ind = 1; % % Should always start from 1
ref_day = days(ref_day_ind); 

contrast_ref = 15;
contr = 3;

% Set initial values
transparency_red = 1;
transparency_green = 1; 
overlay_ref = .5; 
overlay_bkg = .5;

ref_angle = 90;    % Absolute angle [deg] for reference day
                    % This is the angle between the head line of the mouse
                    % (head to tip) and the side of the image sensor board
                    % (last pixel row to first pixel row). It is positive
                    % clockwise.
rot_matrix_ref = [cos(ref_angle*pi/180) -sin(ref_angle*pi/180); sin(ref_angle*pi/180) cos(ref_angle*pi/180)];                    
                    


XSHIFT = cell(N_days, N_ss);
YSHIFT = cell(N_days, N_ss);
ROT = cell(N_days, N_ss);

% Set first day's parameters to 0
for curr_ss = sessions{1}
    curr_ss_ind = find(sessions{1}==curr_ss);
    if ~isempty(CELLS_DYN{1,curr_ss_ind})
        XSHIFT{1,curr_ss_ind} = 0;
        YSHIFT{1,curr_ss_ind} = 0;
        ROT{1,curr_ss_ind} = 0;
    end
end


CELLS_REF2 = zeros(height,width);
CELLS_REF = zeros(height,width);

for curr_ss = sessions{ref_day_ind}
    curr_ss_ind = find(sessions{ref_day_ind}==curr_ss);
    if ~isempty(CELLS_DYN{ref_day_ind,curr_ss_ind})
        CELLS_REF = CELLS_REF+contrast_ref*(1/N_ss).*CELLS{ref_day_ind,curr_ss_ind};
    end
end
CELLS_REF(CELLS_REF>1)=1;

% REF_MAP = zeros(height,width);
% for curr_ss = 1%sessions{ref_day_ind}
%     curr_ss_ind = find(sessions{ref_day_ind}==curr_ss);
%     if ~isempty(CELLS_DYN{ref_day_ind,curr_ss_ind})
%         REF_MAP = REF_MAP + double(MIN_FRAME{ref_day_ind,curr_ss_ind});
%     end
% end
% REF_MAP = 0.3.*REF_MAP./max(max(REF_MAP));


X_REF = CENT_XY_DAY{1}{ref_day_ind};
Y_REF = CENT_XY_DAY{2}{ref_day_ind};

% curr_ref_day = ref_day;
for curr_day = days(2:end)
    curr_day_ind = find(days==curr_day);
    
    
    % curr_ref_day = days(curr_day_ind-1);
    % ref_day_ind = find(days==curr_ref_day);


    CELLS_MAP = zeros(height,width);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            CELLS_MAP = CELLS_MAP+contr*(1/length(sessions{curr_day_ind})).*CELLS{curr_day_ind,curr_ss_ind};
        end
    end
    CELLS_MAP(CELLS_MAP>1)=1;
    
    
%     CURR_MAP = zeros(height,width);
%     for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
%             CURR_MAP = CURR_MAP + double(MIN_FRAME{curr_day_ind,curr_ss_ind});
%         end
%     end
%     CURR_MAP = 0.3.*CURR_MAP./max(max(CURR_MAP));
 
    
    % Use full day
    
    [xshift,yshift,rot, transparency_red, transparency_green, overlay_ref, overlay_bkg] = rotationGUI(CELLS_MAP, CELLS_REF, CENT_XY_DAY{1}{curr_day_ind}, CENT_XY_DAY{2}{curr_day_ind}, X_REF, Y_REF, MIN_FRAME{curr_day_ind,1}, MIN_FRAME{ref_day_ind,1},curr_day_ind,1, transparency_red, transparency_green, overlay_ref, overlay_bkg); 
    % [xshift,yshift,rot] = rotationGUI(CURR_MAP, REF_MAP, CENT_XY_DAY{1}{curr_day_ind}, CENT_XY_DAY{2}{curr_day_ind}, X_REF, Y_REF, MIN_FRAME{curr_day_ind,1}, MIN_FRAME{ref_day_ind,1},curr_day_ind,1, transparency); 
    
    % rot = rot + ref_angle;
    
  
    
    
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        % if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            % CELLS_REF = (CELLS{ref_day,1} | CELLS{ref_day,2} | CELLS{ref_day,3});
            % % Use each experiment
            % [xshift,yshift,rot] = rotationGUI(CELLS{curr_day_ind,curr_ss_ind}, CELLS_REF, MIN_FRAME{curr_day_ind,curr_ss_ind}, MIN_FRAME{curr_day_ind,curr_ss_ind},curr_day_ind,curr_ss_ind); 
            
            XSHIFT{curr_day_ind,curr_ss_ind} = xshift;
            YSHIFT{curr_day_ind,curr_ss_ind} = yshift;
            ROT{curr_day_ind,curr_ss_ind} = rot;
            
           
        % end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    X_REF0 = CENT_XY_DAY{1}{curr_day_ind};
    Y_REF0 = CENT_XY_DAY{2}{curr_day_ind};
    
    tot_cells = length(X_REF0);
    rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
    
    for k = 1:tot_cells
        coord_x = X_REF0(k);
        coord_y = Y_REF0(k);
        
        
        % Image at the bottom of the GRIN lens is flipped l/r and u/d by
        % the GRIN lens, then l/r and u/d by the miniscope, then l/r only
        % by the image sensor, then l/r again by the NeuView viewer (but
        % not in the saved tiff).
        
        
        % Flip LR (image is flipped on the back of the image sensor, first pixel is on the right)
        % The image on the NeuView viewer is flipped lr (so it looks
        % straight at the bottom of the GRIN lens), but the saved image needs to be flipped lr.
        % When plotting images in Matlab do not reverse y direction (first pixel should be
        % on top)
        
      % coord_x = -coord_x + width +1; 
        
        % Rotate
        [coord] = [width/2; height/2] + rot_matrix*([coord_x; coord_y]-[width/2; height/2]);
        coord_x = coord(1)-xshift;
        coord_y = coord(2)-yshift;
        

        % Save
        X_REF0(k) = coord_x;
        Y_REF0(k) = coord_y;

    end
    
    X_REF = [X_REF, X_REF0];
    Y_REF = [Y_REF, Y_REF0];
    
    
    
    
    
    
    
    
    % xshift = -xshift;
    % yshift = -yshift;
    % rot = -rot;
    
    CELLS_REF2 = CELLS_MAP;
    
    % Rotate clockwise (because of the l/r flip)
    if(rot ~= 0)
        CELLS_REF2 = imrotate(CELLS_REF2,rot,'bilinear','crop');
    end
    
    % Shift x
    if (xshift > 0)
        CELLS_REF2(:,1:end-xshift) = CELLS_REF2(:,xshift+1:end);
    elseif (xshift < 0)
        CELLS_REF2(:,-xshift+1:end) = CELLS_REF2(:,1:end+xshift);
    end
    
    % Shift y
    if (yshift > 0)
        CELLS_REF2(1:end-yshift,:) = CELLS_REF2(yshift+1:end,:);
    elseif (yshift < 0)
        CELLS_REF2(-yshift+1:end,:) = CELLS_REF2(1:end+yshift,:);
    end
    CELLS_REF = CELLS_REF + CELLS_REF2;
    CELLS_REF(CELLS_REF>1)=1;
    
    
    % xshift = -xshift;
    % yshift = -yshift;
    % rot = -rot; 
    
    
            
    
end








IM_REG = cell(4,1);
IM_REG{1} = XSHIFT;
IM_REG{2} = YSHIFT;
IM_REG{3} = ROT;
IM_REG{4} = [ref_day, ref_angle];
        
        

% Save data

eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{16} = IM_REG;'])


for jj = 1:length(days)
    fprintf('Absolute rotation day %d: %d\n',days(jj),ROT{jj,1}+ref_angle)
end




fprintf('Calculate optimal image shift/rotation is Done.\n');

% tmp = cell2mat(ROT);
% tmp(:,1)




















%% Create registered cell mask

CELLS_MAP = zeros(height,width);
for curr_day = days
    curr_day_ind = find(days==curr_day);

    rot = -ROT{curr_day_ind,1};
    xshift = -XSHIFT{curr_day_ind,1};
    yshift = YSHIFT{curr_day_ind,1};
    
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            CURR_CELLS = CELLS{curr_day_ind,curr_ss_ind};
            
            CURR_CELLS = flip(CURR_CELLS,2);
             % Rotate clockwise (because of the l/r flip)
            if(rot ~= 0)     
                CURR_CELLS = imrotate(CURR_CELLS,rot,'bilinear','crop');
            end
              % Shift x
            if (xshift > 0)
                CURR_CELLS(:,1:end-xshift) = CURR_CELLS(:,xshift+1:end);
            elseif (xshift < 0)
                CURR_CELLS(:,-xshift+1:end) = CURR_CELLS(:,1:end+xshift);
            end
            % Shift y
            if (yshift > 0)
                CURR_CELLS(1:end-yshift,:) = CURR_CELLS(yshift+1:end,:);
            elseif (yshift < 0)      
                CURR_CELLS(-yshift+1:end,:) = CURR_CELLS(1:end+yshift,:);
            end
            
            % Rotate by reference angle (rotate clockwise)
            if(ref_angle ~= 0)     
                CURR_CELLS = imrotate(CURR_CELLS,-ref_angle,'bilinear','crop');
            end
            
            
            CELLS_MAP = CELLS_MAP+0.2*(1/length(days))*(1/N_ss).*CURR_CELLS;
        end
    end
    
end

CELLS_MAP = CELLS_MAP./max(max(CELLS_MAP));

figure; imagesc(1-CELLS_MAP); colormap(gray); caxis([0.88 1]);
axis equal;







%%  Register different days (rotate cell mask and cell position)


area_min = 8;
area_max = 80;
merge_dist = 4;
merge_dist2 = 4;
manual_select = 1;



% % Rotate and shift cell position for each day

CENT_X_DAY_R = cell(N_days,1);
CENT_Y_DAY_R = cell(N_days,1);
for curr_day = days
    curr_day_ind = find(days==curr_day);
    
    
    % Rotate clockwise (because of the l/r flip)
    % Use negative rotation and xshift as the image will be flipped lr
    rot = -ROT{curr_day_ind,1};
    xshift = -XSHIFT{curr_day_ind,1};
    yshift = YSHIFT{curr_day_ind,1};
    
    rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
    

    tot_cells = length(CENT_XY_DAY{1}{curr_day_ind});
    for k = 1:tot_cells
        coord_x = CENT_XY_DAY{1}{curr_day_ind}(k);
        coord_y = CENT_XY_DAY{2}{curr_day_ind}(k);
        
        
        % Image at the bottom of the GRIN lens is flipped l/r and u/d by
        % the GRIN lens, then l/r and u/d by the miniscope, then l/r only
        % by the image sensor, then l/r again by the NeuView viewer (but
        % not in the saved tiff).
        
        
        % Flip LR (image is flipped on the back of the image sensor, first pixel is on the right)
        % The image on the NeuView viewer is flipped lr (so it looks
        % straight at the bottom of the GRIN lens), but the saved image needs to be flipped lr.
        % When plotting images in Matlab do not reverse y direction (first pixel should be
        % on top)
   
        coord_x = -coord_x + width +1; 

        % Rotate
        [coord] = [width/2; height/2] + rot_matrix*([coord_x; coord_y]-[width/2; height/2]);
        coord_x = coord(1)-xshift;
        coord_y = coord(2)-yshift;
        [coord] = [width/2; height/2] + rot_matrix_ref*([coord_x; coord_y]-[width/2; height/2]);
        coord_x = coord(1);
        coord_y = coord(2);
     
        % Save
        CENT_X_DAY_R{curr_day_ind}(k) = coord_x;
        CENT_Y_DAY_R{curr_day_ind}(k) = coord_y;
      
    end
end    




% % Rotate and shift cell position for each experiment

CENT_X_R = cell(N_days,N_ss);
CENT_Y_R = cell(N_days,N_ss);
for curr_day = days
    curr_day_ind = find(days==curr_day);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
    
        % Rotate clockwise (because of the l/r flip)
        % Use negative rotation and xshift as the image will be flipped lr
        rot = -ROT{curr_day_ind,curr_ss_ind};
        xshift = -XSHIFT{curr_day_ind,curr_ss_ind};
        yshift = YSHIFT{curr_day_ind,curr_ss_ind};

        rot_matrix = [cos(-rot*pi/180) -sin(-rot*pi/180); sin(-rot*pi/180) cos(-rot*pi/180)];
        tot_cells = length(CENT_X{curr_day_ind,curr_ss_ind});
        for k = 1:tot_cells
            coord_x = CENT_X{curr_day_ind,curr_ss_ind}(k);
            coord_y = CENT_Y{curr_day_ind,curr_ss_ind}(k);


            % Image at the bottom of the GRIN lens is flipped l/r and u/d by
            % the GRIN lens, then l/r and u/d by the miniscope, then l/r only
            % by the image sensor, then l/r again by the NeuView viewer (but
            % not in the saved tiff).


            % Flip LR (image is flipped on the back of the image sensor, first pixel is on the right)
            % The image on the NeuView viewer is flipped lr (so it looks
            % straight at the bottom of the GRIN lens), but the saved image needs to be flipped lr.
            % When plotting images in Matlab do not reverse y direction (first pixel should be
            % on top)

           coord_x = -coord_x + width +1; 

            % Rotate
            [coord] = [width/2; height/2] + rot_matrix*([coord_x; coord_y]-[width/2; height/2]);
            coord_x = coord(1)-xshift;
            coord_y = coord(2)-yshift;
            [coord] = [width/2; height/2] + rot_matrix_ref*([coord_x; coord_y]-[width/2; height/2]);
            coord_x = coord(1);
            coord_y = coord(2);

            % Save
            CENT_X_R{curr_day_ind,curr_ss_ind}(k) = coord_x;
            CENT_Y_R{curr_day_ind,curr_ss_ind}(k) = coord_y;
        end
    end
end    
% figure; hold on;
% % plot(CENT_X_DAY_R{1},CENT_Y_DAY_R{1},'r.','MarkerSize',10);
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind} == curr_ss);
%         plot(CENT_X_R{curr_day_ind,curr_ss_ind},CENT_Y_R{curr_day_ind,curr_ss_ind},'o','color',COLORS(curr_day_ind,:));
%     end
% end
% axis ij
% axis equal
% xlim([0, width]);
% ylim([0, height]);






% figure; hold on;
% imagesc(1-CELLS_MAP); colormap(gray); caxis([0.75 1]);
% % plot(CENT_X_DAY_R{1},CENT_Y_DAY_R{1},'r.','MarkerSize',10);
% for k = days
%     curr_day_ind = find(days==k);
% 	plot(CENT_X_DAY_R{curr_day_ind},CENT_Y_DAY_R{curr_day_ind},'o','color',COLORS(curr_day_ind,:),'markersize',4);
% end
% axis ij
% axis equal
% xlim([0, width]);
% ylim([0, height]);
% plot(global_cells(:,1), global_cells(:,2),'om','markersize',20);





figure; imagesc(1-CELLS_MAP); colormap(gray); caxis([0.9 1]);
hold on;
% plot(global_cells(:,1), global_cells(:,2),'or','markersize',8,'linewidth',2);
for k = days
    curr_day_ind = find(days==k);
    curr_cells_ind = 1:length(CENT_X_DAY_R{curr_day_ind});
    curr_cells_ind(curr_cells_ind == 0) = [];
	plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'.','color',COLORS(curr_day_ind,:),'markersize',15);

end

axis equal




%% Register different days: identify global cells

   



% Find threshold for global cells
thresh = findThresh(CELLS_MAP, area_min, area_max);
% thresh = 0.01;


CELLS_MAP2 = CELLS_MAP;
CELLS_MAP2(CELLS_MAP2 <= thresh) = 0;
CELLS_MAP2 = logical(CELLS_MAP2);
CELLS_MAP2 = double(CELLS_MAP2);

figure; hold on;
imagesc(1-CELLS_MAP2); colormap(gray); caxis([-3 1]);
axis ij
axis equal
xlim([0, width]);
ylim([0, height]);

% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
[labeledImage, numberOfBlobs] = bwlabel(CELLS_MAP2);
blobMeasurements = regionprops(labeledImage, 'area','Centroid','PixelList');
% Get all the areas
allAreas = [blobMeasurements.Area]; % No semicolon so it will print to the command window.


[B,L] = bwboundaries(CELLS_MAP2,'noholes');

hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)
end



% Create global cell vector

tmp = find(allAreas>=area_max);
fprintf('Found %d large cells.\n',length(tmp));

global_cells = [];
for k = 1 : numberOfBlobs           % Loop through all blobs.
    thisCentroid = [blobMeasurements(k).Centroid(1), blobMeasurements(k).Centroid(2)];
    message = sprintf('%d', allAreas(k));
    if allAreas(k) >= area_min
        text(thisCentroid(1)+2, thisCentroid(2), message, 'Color', 'r');
        if allAreas(k) <= area_max
            
            global_cells = [global_cells; [thisCentroid(1), thisCentroid(2)]];
        else
            curr_region = blobMeasurements(k).PixelList;
            
            curr_cells_x = [];
            curr_cells_y = [];
             
            curr_area = CELLS_MAP(min(curr_region(:,2)):max(curr_region(:,2)),min(curr_region(:,1)):max(curr_region(:,1)));
            
            % curr_x = CURR_CENT_X(idx_x) - min(curr_region(1,:))+1;
            % curr_y = CURR_CENT_Y(idx_y) - min(curr_region(2,:))+1;
            
            winsize = [ 40 40 30*(max(curr_region(:,1))-min(curr_region(:,1))) 30*(max(curr_region(:,2))-min(curr_region(:,2)))];
            
            tmp1 = figure;%figure('Position',winsize);
            
            % imagesc(curr_area); caxis([0 length(sessions{curr_day_ind})]);  hold on;
            imagesc(curr_area);
            caxis([0 0.6.*max(max(CELLS_MAP))]);
            % colormap(gray);
            hold on;
            % old_cells = plot(curr_x,curr_y, '.r','MarkerSize',50);
            
            [x_pos,y_pos,button] = ginput;
            
            new_x = x_pos(find(button==1)) + min(curr_region(:,1))-1;
            new_y = y_pos(find(button==1)) + min(curr_region(:,2))-1;
            del_x = x_pos(find(button==3)) + min(curr_region(:,1))-1;
            del_y = y_pos(find(button==3)) + min(curr_region(:,2))-1;
            
            %             old_cells.Marker = 'none';
            %             refresh(tmp1);
            %             plot(new_x, new_y, '.m','MarkerSize',50);

            %                     % Remove closest cell to the selected points (+- 0.5 pixel)
            %                     for m = 1:length(del_x)
            %                         r_idx = find(CURR_CENT_X>del_x(m)-0.5 & CURR_CENT_X<del_x(m)+0.5 & CURR_CENT_Y >del_y(m)-0.5 & CURR_CENT_Y <del_y(m)+0.5);
            %                         NEW_CENT_X(r_idx) = nan;
            %                         NEW_CENT_Y(r_idx) = nan;
            %                     end 
            
            % Add new cells

            global_cells = [global_cells; [new_x, new_y]];
            close(tmp1);
        end
    end
end
plot(global_cells(:,1), global_cells(:,2),'og');








% % Match cells

fprintf('Found %d global cells.\n',size(global_cells,1));

CELLS_IDX = zeros(size(global_cells,1),length(days));
remove_cells = cell(length(days),1);
for curr_day = days
    curr_day_ind = find(days==curr_day);
    curr_cells = [CENT_X_DAY_R{curr_day_ind}', CENT_Y_DAY_R{curr_day_ind}'];
    curr_remove_cells = [];
    curr_ind_cum = [];
    for k = 1:size(global_cells,1)
        curr_match_ind = [];
        curr_match_norm = [];
        for kk = 1:size(curr_cells,1)
            curr_norm = norm(curr_cells(kk,:)-global_cells(k,:),2);
            if  (curr_norm <= merge_dist) && (~any(kk==curr_ind_cum))
                curr_match_ind = [curr_match_ind; kk];
                curr_match_norm = [curr_match_norm; curr_norm];
            end
        end
        
        if length(curr_match_ind) > 1
            [aa bb] = sort(curr_match_norm,'ascend');
            curr_remove_cells = [curr_remove_cells; curr_match_ind(bb(2:end))]; 
        else
            bb = 1;
        end
        
        if length(curr_match_ind) >= 1
            curr_ind_cum = [curr_ind_cum; curr_match_ind(bb(1))];
            CELLS_IDX(k,curr_day_ind) = curr_match_ind(bb(1)); 
        end
        
    end
    remove_cells{curr_day_ind} = curr_remove_cells;
end

% % Remove additional cells close to global cells
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     CENT_X_DAY_R{curr_day_ind}(remove_cells{curr_day_ind}) = [];
%     CENT_Y_DAY_R{curr_day_ind}(remove_cells{curr_day_ind}) = [];
% end
% 
% % Identify cells in a merge_dist radius from global cells as the same
% CELLS_IDX = zeros(size(global_cells,1),length(days));
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     curr_cells = [CENT_X_DAY_R{curr_day_ind}', CENT_Y_DAY_R{curr_day_ind}'];
%     curr_remove_cells = [];
%     curr_ind_cum = [];
%     for k = 1:size(global_cells,1)
%         % curr_match = [];
%         curr_match_ind = [];
%         curr_match_norm = [];
%         for kk = 1:size(curr_cells,1)
%             curr_norm = norm(curr_cells(kk,:)-global_cells(k,:),2);
%             if  (curr_norm <= merge_dist) && (~any(kk==curr_ind_cum))
%                 % curr_match = [curr_match; curr_cells(kk,:)];
%                 curr_match_ind = [curr_match_ind; kk];
%                 curr_match_norm = [curr_match_norm; curr_norm];
%             end
%         end
%         if length(curr_match_ind) == 1
%             CELLS_IDX(k,curr_day_ind) = curr_match_ind; 
%             curr_ind_cum = [curr_ind_cum; curr_match_ind];
%         end
%     end
% end










figure; imagesc(1-CELLS_MAP); colormap(gray); caxis([0.7 1]);
hold on;
plot(global_cells(:,1), global_cells(:,2),'or','markersize',8,'linewidth',2);
for k = days
    curr_day_ind = find(days==k);
    curr_cells_ind = CELLS_IDX(1:size(global_cells,1),curr_day_ind);
    curr_cells_ind(curr_cells_ind == 0) = [];
	plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'.','color',COLORS(curr_day_ind,:),'markersize',15);
   
    % curr_cells_ind = remove_cells{curr_day_ind};
	% plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'x','color',COLORS(curr_day_ind,:),'markersize',4);
    
    curr_cells_ind = CELLS_IDX(1:size(global_cells,1),curr_day_ind);
    curr_cells_ind(curr_cells_ind == 0) = [];
    curr_cells_ind = setdiff([1:length(CENT_X_DAY_R{curr_day_ind})],curr_cells_ind);
	plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'x','color',COLORS(curr_day_ind,:),'markersize',5);
    
end






% Remove global cells with no match
rm_idx = [];
for k = 1:size(CELLS_IDX,1)
    if sum(CELLS_IDX(k,:)) == 0
        rm_idx = [rm_idx; k];
    end
end
global_cells(rm_idx,:) = [];
CELLS_IDX(rm_idx,:) = [];

plot(global_cells(:,1), global_cells(:,2),'ob','markersize',8,'linewidth',2);
axis ij;
axis equal;


%% Manual selection


if manual_select   
    NEW_CENT_X = global_cells(:,1);
    NEW_CENT_Y = global_cells(:,2);

    winsize = [ 0 0 1500 1500];
    curr_fig = figure('Position',winsize);        
    hold on;
    imagesc(1-CELLS_MAP); colormap(gray); caxis([0.86   1]);  
    xlim([0, width]); ylim([0, height]);
    plot(NEW_CENT_X,NEW_CENT_Y,'or')
    for k = days
        
        curr_day_ind = find(days==k);
        
        
        curr_cells_ind = CELLS_IDX(1:size(global_cells,1),curr_day_ind);
        curr_cells_ind(curr_cells_ind == 0) = [];
        % plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'.','color',0.5.*COLORS(curr_day_ind,:),'markersize',15);
        
        
        
        curr_cells_ind2 = setdiff([1:length(CENT_X_DAY_R{curr_day_ind})]', curr_cells_ind);
        
        % plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'.','color',COLORS(curr_day_ind,:),'markersize',15);
        
  
        plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind2),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind2),'.','color',COLORS(curr_day_ind,:),'markersize',15);
        
    end
    axis ij
    axis equal

    % plot(CURR_CENT_X,CURR_CENT_Y,'or')
   
            
    [x_pos,y_pos,button] = ginput;

    

    new_x = x_pos;
    new_y = y_pos;




    % Add new cells

    NEW_CENT_X = [NEW_CENT_X; new_x];
    NEW_CENT_Y = [NEW_CENT_Y; new_y];

  
    
    new_global_x = NEW_CENT_X;
    new_global_y = NEW_CENT_Y;
    
    
    fprintf('Added %d cells.\n', length(new_x));
    
    
    global_cells = zeros(length(new_global_x),2);
    global_cells(:,1) = new_global_x;
    global_cells(:,2) = new_global_y;
    
    
 end






close(curr_fig);


CELLS_IDX = zeros(size(global_cells,1),length(days));
remove_cells = cell(length(days),1);
for curr_day = days
    curr_day_ind = find(days==curr_day);
    curr_cells = [CENT_X_DAY_R{curr_day_ind}', CENT_Y_DAY_R{curr_day_ind}'];
    curr_remove_cells = [];
    curr_ind_cum = [];
    for k = 1:size(global_cells,1)
        curr_match_ind = [];
        curr_match_norm = [];
        for kk = 1:size(curr_cells,1)
            curr_norm = norm(curr_cells(kk,:)-global_cells(k,:),2);
            if  (curr_norm <= merge_dist) && (~any(kk==curr_ind_cum))
                curr_match_ind = [curr_match_ind; kk];
                curr_match_norm = [curr_match_norm; curr_norm];
            end
        end
        
        if length(curr_match_ind) > 1
            [aa bb] = sort(curr_match_norm,'ascend');
            curr_remove_cells = [curr_remove_cells; curr_match_ind(bb(2:end))]; 
        else
            bb = 1;
        end
        
        if length(curr_match_ind) >= 1
            curr_ind_cum = [curr_ind_cum; curr_match_ind(bb(1))];
            CELLS_IDX(k,curr_day_ind) = curr_match_ind(bb(1)); 
        end
        
    end
    remove_cells{curr_day_ind} = curr_remove_cells;
end

% % Remove additional cells close to global cells
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     CENT_X_DAY_R{curr_day_ind}(remove_cells{curr_day_ind}) = [];
%     CENT_Y_DAY_R{curr_day_ind}(remove_cells{curr_day_ind}) = [];
% end
% 
% % Identify cells in a merge_dist radius from global cells as the same
% CELLS_IDX = zeros(size(global_cells,1),length(days));
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     curr_cells = [CENT_X_DAY_R{curr_day_ind}', CENT_Y_DAY_R{curr_day_ind}'];
%     curr_remove_cells = [];
%     curr_ind_cum = [];
%     for k = 1:size(global_cells,1)
%         % curr_match = [];
%         curr_match_ind = [];
%         curr_match_norm = [];
%         for kk = 1:size(curr_cells,1)
%             curr_norm = norm(curr_cells(kk,:)-global_cells(k,:),2);
%             if  (curr_norm <= merge_dist) && (~any(kk==curr_ind_cum))
%                 % curr_match = [curr_match; curr_cells(kk,:)];
%                 curr_match_ind = [curr_match_ind; kk];
%                 curr_match_norm = [curr_match_norm; curr_norm];
%             end
%         end
%         if length(curr_match_ind) == 1
%             CELLS_IDX(k,curr_day_ind) = curr_match_ind; 
%             curr_ind_cum = [curr_ind_cum; curr_match_ind];
%         end
%     end
% end




fprintf('Total %d global cells.\n', size(CELLS_IDX,1));

%% Plot

min_thresh = 5; % Highlight cells closer than min_thresh pixels

figure; hold on;
imagesc(1-CELLS_MAP); colormap(gray); caxis([0.75 1]);
% plot(CENT_X_DAY_R{1},CENT_Y_DAY_R{1},'r.','MarkerSize',10);
for k = days
    
    curr_day_ind = find(days==k);
    
    
    curr_cells_ind = CELLS_IDX(1:size(global_cells,1),curr_day_ind);
    curr_cells_ind(curr_cells_ind == 0) = [];
    plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'.','color',COLORS(curr_day_ind,:),'markersize',15);
    
    % curr_cells_ind = remove_cells{curr_day_ind};
    % plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'x','color',COLORS(curr_day_ind,:),'markersize',4);
    
    curr_cells_ind = CELLS_IDX(size(global_cells,1):end,curr_day_ind);
    curr_cells_ind(curr_cells_ind == 0) = [];
    curr_cells_ind = setdiff([1:length(CENT_X_DAY_R{curr_day_ind})],curr_cells_ind);
    plot(CENT_X_DAY_R{curr_day_ind}(curr_cells_ind),CENT_Y_DAY_R{curr_day_ind}(curr_cells_ind),'x','color',COLORS(curr_day_ind,:),'markersize',5);
    
end

plot(global_cells(:,1), global_cells(:,2),'om','markersize',18,'linewidth',1);


for k = 1:size(global_cells,1)
    text(global_cells(k,1)+2, global_cells(k,2)-1,num2str(k),'color','g');
end
% plot(global_cellsx(size(global_cells,1):end), global_cellsy(size(global_cells,1):end),'ok','markersize',17,'linewidth',2);


% Highlight neurons that are too close
dist = nan(size(global_cells,1),size(global_cells,1));
for m = 1:size(global_cells,1)-1
    for n = m+1:size(global_cells,1)
        dist(n,m) = norm([global_cells(m,1); global_cells(m,2)]-[global_cells(n,1); global_cells(n,2)],2);
    end
end

dist(dist>min_thresh) = nan;
[cells1,cells2] = find(~isnan(dist));

for k = 1:length(cells1)
    
    plot(global_cells(cells1(k),1),global_cells(cells1(k),2),'o','color',COLORS(2,:),'markersize',18,'linewidth',3);
    plot(global_cells(cells2(k),1),global_cells(cells2(k),2),'o','color',COLORS(2,:),'markersize',18,'linewidth',3);
    
end




axis ij
axis equal
xlim([0, width]);
ylim([0, height]);


%% Identify same cells across different days and build global cell map

GLOBAL_CELLS = cell(3,1);
CENT_X_GLOBAL = [];
CENT_Y_GLOBAL = [];
IDX_MAP = [];


if (~manual_select)
    
    
    % Add global cells where no daily cells were identified
    tmp = sum(CELLS_IDX,2);
    tmp2 = find(tmp == 0);
    for k = 1:length(tmp2)
        for curr_day = days
            curr_day_ind = find(days==curr_day);
            CENT_X_DAY_R{curr_day_ind} = [CENT_X_DAY_R{curr_day_ind}, global_cells(tmp2(k),1)];
            CENT_Y_DAY_R{curr_day_ind} = [CENT_Y_DAY_R{curr_day_ind}, global_cells(tmp2(k),2)];
            CELLS_IDX(tmp2(k),curr_day_ind) = length(CENT_X_DAY_R{curr_day_ind});
        end
    end
    
    
    % Save temp
    eval(['mouse_tmp = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL;'])
    
    other_cent_x = CENT_X_DAY_R;
    other_cent_y = CENT_Y_DAY_R;
    other_cent_x0 = CENT_X_DAY;
    other_cent_y0 = CENT_Y_DAY;
    global_indices = cell(length(days),1);
    for curr_day = days
        curr_day_ind = find(days==curr_day);
        global_indices{curr_day_ind} = CELLS_IDX(:,curr_day_ind);
        tmp = global_indices{curr_day_ind};
        tmp(tmp==0) = [];
        other_cent_x0{curr_day_ind}(tmp) = nan;
        other_cent_y0{curr_day_ind}(tmp) = nan;
        other_cent_x{curr_day_ind}(tmp) = nan;
        other_cent_y{curr_day_ind}(tmp) = nan;
    end
    
    CENT_XY_DAY2 = cell(4,1);
    CENT_XY_DAY2{1} = other_cent_x0;
    CENT_XY_DAY2{2} = other_cent_y0;
    CENT_XY_DAY2{3} = other_cent_x;
    CENT_XY_DAY2{4} = other_cent_y;
    mouse_tmp{17} = CENT_XY_DAY2;
    
    
    
    [CENT_X_GLOBAL, CENT_Y_GLOBAL, IDX_MAP] = cell_mask(merge_dist2, mouse_tmp);
    clear mouse_tmp
    % remove nans
    tmp_ind = IDX_MAP>0;
    tmp_ind = sum(tmp_ind,2);
    remove_ind = [];
    for k = 1:length(tmp_ind)
        if tmp_ind(k)==1
            if isnan(CENT_X_GLOBAL(k))
                remove_ind = [remove_ind; k];
            end
        end
    end
    CENT_X_GLOBAL(remove_ind) = [];
    CENT_Y_GLOBAL(remove_ind) = [];
    IDX_MAP(remove_ind,:) = [];
    
    
   
    
end
    








% merge all
global_cellsx = [global_cells(:,1); CENT_X_GLOBAL];
global_cellsy = [global_cells(:,2); CENT_Y_GLOBAL];
CELLS_IDX = [CELLS_IDX; IDX_MAP];




GLOBAL_CELLS{1} = global_cellsx;
GLOBAL_CELLS{2} = global_cellsy;
GLOBAL_CELLS{3} = CELLS_IDX;


fprintf('Found total of %d global cells across %d days.\n',size(CELLS_IDX,1),length(days));




%% Save registered data

eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{12} = CENT_X_R;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{13} = CENT_Y_R;'])

CENT_XY_DAY{3} = CENT_X_DAY_R;
CENT_XY_DAY{4} = CENT_Y_DAY_R;
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{17} = CENT_XY_DAY;'])

% Save global cell position
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{20} = GLOBAL_CELLS;'])     

fprintf('Rotate and shift each cell position for each day is Done.\n');

         








 
%% Rotate and shift frames 


wbar = waitbar(0,'Performing image registration');
for curr_day = days
    curr_day_ind = find(days==curr_day);
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        waitbar(((N_ss*(curr_day_ind-1))+curr_ss_ind)/(N_days*N_ss),wbar,'Performing image registration')
         
        
        if ~isempty(CELLS_DYN{curr_day_ind,curr_ss_ind})
            
            % Load frames for each experiment
            load(strcat(BASEDIR,slsh,'Frames',slsh,'Frames_D',num2str(cell_type),'M',num2str(mouse_numb),'_D',num2str(curr_day),'S',num2str(curr_ss),'.mat'));

            curr_cells_dyn = logical(zeros(height,width,size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)));
            for jj = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
                tmp = bwunpack(CELLS_DYN{curr_day_ind,curr_ss_ind}(:,:,jj));
                curr_cells_dyn(:,:,jj) = tmp(1:height,1:width);
            end
            
            curr_cells_dyn_sum = logical(zeros(height,width,size(CELLS_DYN_SUM{curr_day_ind,curr_ss_ind},3)));
            for jj = 1:size(CELLS_DYN_SUM{curr_day_ind,curr_ss_ind},3)
                tmp = bwunpack(CELLS_DYN_SUM{curr_day_ind,curr_ss_ind}(:,:,jj));
                curr_cells_dyn_sum(:,:,jj) = tmp(1:height,1:width);
            end

            
            rot = -ROT{curr_day_ind,curr_ss_ind};
            xshift = -XSHIFT{curr_day_ind,curr_ss_ind};
            yshift = YSHIFT{curr_day_ind,curr_ss_ind};




            % Flip LR (image is flipped on the back of the image sensor, first pixel is on the right)
            % When plotting images in Matlab reverse y direction (first pixel should be
            % on top)

            frames = flip(frames,2);

            curr_cells_dyn = flip(curr_cells_dyn,2);
            curr_cells_dyn_sum = flip(curr_cells_dyn,2);
            % CELLS_DYN{curr_day_ind,curr_ss_ind} = flip(CELLS_DYN{curr_day_ind,curr_ss_ind},2);

            CELLS{curr_day_ind,curr_ss_ind} = flip(CELLS{curr_day_ind,curr_ss_ind},2);
            CELLS_FILT{curr_day_ind,curr_ss_ind} = flip(CELLS_FILT{curr_day_ind,curr_ss_ind},2);
            MIN_FRAME{curr_day_ind,curr_ss_ind} = flip(MIN_FRAME{curr_day_ind,curr_ss_ind},2);
            MIN_FRAME_FILT{curr_day_ind,curr_ss_ind} = flip(MIN_FRAME_FILT{curr_day_ind,curr_ss_ind},2);
            AVG_FRAME{curr_day_ind,curr_ss_ind} = flip(AVG_FRAME{curr_day_ind,curr_ss_ind},2);





            % Rotate clockwise (because of the l/r flip)
            if(rot ~= 0)     
                for k = 1:size(frames,3)
                    frames(:,:,k) = imrotate(frames(:,:,k),rot,'bilinear','crop');

                end 
                for k = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
                    curr_cells_dyn(:,:,k) = imrotate(curr_cells_dyn(:,:,k),rot,'bilinear','crop');
                    curr_cells_dyn_sum(:,:,k) = imrotate(curr_cells_dyn_sum(:,:,k),rot,'bilinear','crop');
                    % CELLS_DYN{curr_day_ind,curr_ss_ind}(:,:,k) = imrotate(CELLS_DYN{curr_day_ind,curr_ss_ind}(:,:,k),-rot,'bilinear','crop');
                end

                CELLS{curr_day_ind,curr_ss_ind} = imrotate(CELLS{curr_day_ind,curr_ss_ind},rot,'bilinear','crop');
                CELLS_FILT{curr_day_ind,curr_ss_ind} = imrotate(CELLS_FILT{curr_day_ind,curr_ss_ind},rot,'bilinear','crop');

                MIN_FRAME{curr_day_ind,curr_ss_ind} = imrotate(MIN_FRAME{curr_day_ind,curr_ss_ind},rot,'bilinear','crop');
                MIN_FRAME_FILT{curr_day_ind,curr_ss_ind} = imrotate(MIN_FRAME_FILT{curr_day_ind,curr_ss_ind},rot,'bilinear','crop');
                AVG_FRAME{curr_day_ind,curr_ss_ind} = imrotate(AVG_FRAME{curr_day_ind,curr_ss_ind},rot,'bilinear','crop');
            end



            % Shift x
            if (xshift > 0)
                for k = 1:size(frames,3)
                    frames(:,1:end-xshift,k) = frames(:,xshift+1:end,k);    
                end 
                for k = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
                    curr_cells_dyn(:,1:end-xshift,k) = curr_cells_dyn(:,xshift+1:end,k);
                    curr_cells_dyn_sum(:,1:end-xshift,k) = curr_cells_dyn_sum(:,xshift+1:end,k);
                    % CELLS_DYN{curr_day_ind,curr_ss_ind}(:,1:end-xshift,k) = CELLS_DYN{curr_day_ind,curr_ss_ind}(:,xshift+1:end,k);
                end
                CELLS{curr_day_ind,curr_ss_ind}(:,1:end-xshift) = CELLS{curr_day_ind,curr_ss_ind}(:,xshift+1:end);
                CELLS_FILT{curr_day_ind,curr_ss_ind}(:,1:end-xshift) = CELLS_FILT{curr_day_ind,curr_ss_ind}(:,xshift+1:end);
                MIN_FRAME{curr_day_ind,curr_ss_ind}(:,1:end-xshift) = MIN_FRAME{curr_day_ind,curr_ss_ind}(:,xshift+1:end);
                MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(:,1:end-xshift) = MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(:,xshift+1:end);
                AVG_FRAME{curr_day_ind,curr_ss_ind}(:,1:end-xshift) = AVG_FRAME{curr_day_ind,curr_ss_ind}(:,xshift+1:end);
            elseif (xshift < 0)
                for k = 1:size(frames,3)
                    frames(:,-xshift+1:end,k) = frames(:,1:end+xshift,k);   
                end 
                for k = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
                    curr_cells_dyn(:,-xshift+1:end,k) = curr_cells_dyn(:,1:end+xshift,k);
                    curr_cells_dyn_sum(:,-xshift+1:end,k) = curr_cells_dyn_sum(:,1:end+xshift,k);
                    % CELLS_DYN{curr_day_ind,curr_ss_ind}(:,-xshift+1:end,k) = CELLS_DYN{curr_day_ind,curr_ss_ind}(:,1:end+xshift,k);
                end
                CELLS{curr_day_ind,curr_ss_ind}(:,-xshift+1:end) = CELLS{curr_day_ind,curr_ss_ind}(:,1:end+xshift);
                CELLS_FILT{curr_day_ind,curr_ss_ind}(:,-xshift+1:end) = CELLS_FILT{curr_day_ind,curr_ss_ind}(:,1:end+xshift);
                MIN_FRAME{curr_day_ind,curr_ss_ind}(:,-xshift+1:end) = MIN_FRAME{curr_day_ind,curr_ss_ind}(:,1:end+xshift);
                MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(:,-xshift+1:end) = MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(:,1:end+xshift);
                AVG_FRAME{curr_day_ind,curr_ss_ind}(:,-xshift+1:end) = AVG_FRAME{curr_day_ind,curr_ss_ind}(:,1:end+xshift);
            end


            % Shift y
            if (yshift > 0)
                for k = 1:size(frames,3)
                    frames(1:end-yshift,:,k) = frames(yshift+1:end,:,k);
                end 
                for k = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
                    curr_cells_dyn(1:end-yshift,:,k) = curr_cells_dyn(yshift+1:end,:,k);
                    curr_cells_dyn_sum(1:end-yshift,:,k) = curr_cells_dyn_sum(yshift+1:end,:,k);
                    % CELLS_DYN{curr_day_ind,curr_ss_ind}(1:end-yshift,:,k) = CELLS_DYN{curr_day_ind,curr_ss_ind}(yshift+1:end,:,k);
                end
                CELLS{curr_day_ind,curr_ss_ind}(1:end-yshift,:) = CELLS{curr_day_ind,curr_ss_ind}(yshift+1:end,:);
                CELLS_FILT{curr_day_ind,curr_ss_ind}(1:end-yshift,:) = CELLS_FILT{curr_day_ind,curr_ss_ind}(yshift+1:end,:);
                MIN_FRAME{curr_day_ind,curr_ss_ind}(1:end-yshift,:) = MIN_FRAME{curr_day_ind,curr_ss_ind}(yshift+1:end,:);
                MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(1:end-yshift,:) = MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(yshift+1:end,:);
                AVG_FRAME{curr_day_ind,curr_ss_ind}(1:end-yshift,:) = AVG_FRAME{curr_day_ind,curr_ss_ind}(yshift+1:end,:);
            elseif (yshift < 0)      
                for k = 1:size(frames,3)
                    frames(-yshift+1:end,:,k) = frames(1:end+yshift,:,k);
                end 
                for k = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
                    curr_cells_dyn(-yshift+1:end,:,k) = curr_cells_dyn(1:end+yshift,:,k);
                    curr_cells_dyn_sum(-yshift+1:end,:,k) = curr_cells_dyn_sum(1:end+yshift,:,k);
                    % CELLS_DYN{curr_day_ind,curr_ss_ind}(-yshift+1:end,:,k) = CELLS_DYN{curr_day_ind,curr_ss_ind}(1:end+yshift,:,k);
                end
                CELLS{curr_day_ind,curr_ss_ind}(-yshift+1:end,:) = CELLS{curr_day_ind,curr_ss_ind}(1:end+yshift,:);
                CELLS_FILT{curr_day_ind,curr_ss_ind}(-yshift+1:end,:) = CELLS_FILT{curr_day_ind,curr_ss_ind}(1:end+yshift,:);
                MIN_FRAME{curr_day_ind,curr_ss_ind}(-yshift+1:end,:) = MIN_FRAME{curr_day_ind,curr_ss_ind}(1:end+yshift,:);
                MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(-yshift+1:end,:) = MIN_FRAME_FILT{curr_day_ind,curr_ss_ind}(1:end+yshift,:);
                AVG_FRAME{curr_day_ind,curr_ss_ind}(-yshift+1:end,:) = AVG_FRAME{curr_day_ind,curr_ss_ind}(1:end+yshift,:);
            end
            
            
            % Rotate by reference angle (rotate clockwise)
            if(ref_angle ~= 0)     
                for k = 1:size(frames,3)
                    frames(:,:,k) = imrotate(frames(:,:,k),-ref_angle,'bilinear','crop');

                end 
                for k = 1:size(CELLS_DYN{curr_day_ind,curr_ss_ind},3)
                    curr_cells_dyn(:,:,k) = imrotate(curr_cells_dyn(:,:,k),-ref_angle,'bilinear','crop');
                    curr_cells_dyn_sum(:,:,k) = imrotate(curr_cells_dyn_sum(:,:,k),-ref_angle,'bilinear','crop');
                    % CELLS_DYN{curr_day_ind,curr_ss_ind}(:,:,k) = imrotate(CELLS_DYN{curr_day_ind,curr_ss_ind}(:,:,k),-ref_angle,'bilinear','crop');
                end

                CELLS{curr_day_ind,curr_ss_ind} = imrotate(CELLS{curr_day_ind,curr_ss_ind},-ref_angle,'bilinear','crop');
                CELLS_FILT{curr_day_ind,curr_ss_ind} = imrotate(CELLS_FILT{curr_day_ind,curr_ss_ind},-ref_angle,'bilinear','crop');

                MIN_FRAME{curr_day_ind,curr_ss_ind} = imrotate(MIN_FRAME{curr_day_ind,curr_ss_ind},-ref_angle,'bilinear','crop');
                MIN_FRAME_FILT{curr_day_ind,curr_ss_ind} = imrotate(MIN_FRAME_FILT{curr_day_ind,curr_ss_ind},-ref_angle,'bilinear','crop');
                AVG_FRAME{curr_day_ind,curr_ss_ind} = imrotate(AVG_FRAME{curr_day_ind,curr_ss_ind},-ref_angle,'bilinear','crop');
            end
            
           
            % Compress binary data
            CELLS_DYN_PACKED = uint32(zeros(ceil(height/32),width,size(curr_cells_dyn,3)));
            for jj = 1:size(curr_cells_dyn,3)
                CELLS_DYN_PACKED(:,:,jj) = bwpack(curr_cells_dyn(:,:,jj));
            end
            CELLS_DYN{curr_day_ind,curr_ss_ind} = CELLS_DYN_PACKED;
            
            CELLS_DYN_SUM_PACKED = uint32(zeros(ceil(height/32),width,size(curr_cells_dyn,3)));
            for jj = 1:size(curr_cells_dyn,3)
                CELLS_DYN_SUM_PACKED(:,:,jj) = bwpack(curr_cells_dyn(:,:,jj));
            end
            CELLS_DYN_SUM{curr_day_ind,curr_ss_ind} = CELLS_DYN_SUM_PACKED;
      
            
            
            
            
            
            % Save frames to a separate file, one per experiment
            save(strcat(BASEDIR,slsh,'Frames',slsh,'Frames_r_D',num2str(cell_type),'M',num2str(mouse_numb),'_D',num2str(curr_day),'S',num2str(curr_ss),'.mat'),'frames','-v7.3');
            
        end
    end
end
close(wbar)



% Save registered data (overwrite)
% eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{2} = FRAMES;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{3} = MIN_FRAME;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{4} = MIN_FRAME_FILT;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{5} = AVG_FRAME;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{8} = CELLS_DYN;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{9} = CELLS_DYN_SUM;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{10} = CELLS;'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{11} = CELLS_FILT;'])



fprintf('Rotate and shift frames is Done.\n');








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
eval(['save(''..' slsh,'data',slsh,'D',num2str(cell_type),'M',num2str(mouse_numb),'_rfull.mat'',''-v7.3'');'])
fprintf('Done saving process_register workspace to disk.\n');

