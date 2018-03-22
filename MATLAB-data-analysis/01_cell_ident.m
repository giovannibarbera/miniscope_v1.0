close all;
clc
clear


%% Set parameters

    
cell_type           = 1;    
mouse_numb          = 3;

no_cell_id          = 0;    % When true, do not perform cell identification (use for CNMFE)
save_mat_frames     = 0;    % Save frames in Matlab format (one file per session)
save_tiff           = 0;    % Save sessions in tiff format too (stacked, one file per session)


BASEDIR             = '..\data';


DARK_FRAMES_FOLDER = 'b';


N_det               = 99999;    % N_frames 3000 is 5 minutes; % Number of frames used for cell identification
MacOS               = ~ispc;
use_wbar            = true;
remove_bad_frames   = true;
dark_frames         = true;    
register_image_flag = true;
load_timestamps     = true;
show_movie          = false;
height              = 400;  
width               = 400;  
col_depth           = 1024;

cell_size           = 11;   % 9 size of the buffer between positive and negative peak in the gradient
cell_mask_size      = 2;    % 2 x and y cell mask size
T_cell              = 3;    % 5 Cell persistence parameter        
area_min            = 1;    % 1 cell area filtering (may not reflect actual cell area)
min_cell_distance   = 2;    % 2 minimum euclidean distance of centroids to be recognized as cells (1 pixel = 2.75um)
rms_fact            = 4;    % 5 rms factor threshold for gradient peaks identification 


if MacOS
    slsh = '/';
else
    slsh = '\';
end



eartag =  '';





% Days

dir_content = dir(MOUSE_PATH);
dir_count = 0;
if length(dir_content)<=2
    error('ERROR: no data found in specified path.');
else
    for k = 1:length(dir_content)
        if(dir_content(k).isdir) && ~strcmp(dir_content(k).name,'.') && ~strcmp(dir_content(k).name,'..') && ~strcmp(dir_content(k).name,DARK_FRAMES_FOLDER)
            dir_count = dir_count + 1;
        end
    end 
end

N_days = dir_count;


DAY_FOLDERS         = cell(N_days,1);
dir_content = dir(MOUSE_PATH);
dir_count = 0;




% % Sequential days
% days = [1:N_days];



% Or day number reflects folder name

for kk = 1:length(dir_content)
    if (dir_content(kk).isdir) && ~strcmp(dir_content(kk).name,'.') && ~strcmp(dir_content(kk).name,'..') && ~strcmp(dir_content(kk).name,DARK_FRAMES_FOLDER)
        dir_count = dir_count + 1;
        DAY_FOLDERS{dir_count,1} = strcat(dir_content(kk).name,slsh);
        
        % Use folder name to order sessions
        curr_str = dir_content(kk).name;
        numb_ind = regexp(curr_str,'\d');
        days(dir_count) = str2num(curr_str(numb_ind));
    end
end

% Exclude empty folders
rm_days = [];
for k = 1:N_days
    if length(dir(strcat(MOUSE_PATH,DAY_FOLDERS{k}))) == 2
        rm_days = [rm_days; k];
    end
end

if ~isempty(rm_days)
    DAY_FOLDERS2 = cell(N_days-length(rm_days),1);
    tot_ind = 0;
    for k = 1:N_days
        if sum(days(k) == rm_days) == 0
            tot_ind = tot_ind + 1;
            DAY_FOLDERS2{tot_ind} = DAY_FOLDERS{k};
        end
    end
    DAY_FOLDERS = DAY_FOLDERS2;
    days(rm_days) = [];
    N_days = length(days);
end









% Sessions

sessions = cell(N_days,1);
sessions_i = cell(N_days,1);
for k = 1:length(days)
    ss_count = 0;
    dir_content = dir(strcat(MOUSE_PATH, DAY_FOLDERS{k}));
    for kk = 1:length(dir_content)
        if (dir_content(kk).isdir) && ~strcmp(dir_content(kk).name,'.') && ~strcmp(dir_content(kk).name,'..') && ~strcmp(dir_content(kk).name,DARK_FRAMES_FOLDER)
            ss_count = ss_count + 1;
            
            % Use folder name to order sessions
            curr_str = dir_content(kk).name;
            numb_ind = regexp(curr_str,'\d');
            sessions{k} = [sessions{k}, str2num(curr_str(numb_ind))];
            
%             % Or order sessions sequentially
%             sessions{k} = [sessions{k}, str2double(curr_str(numb_ind))];
        end
    end 
    % sessions{k} = [1:ss_count];
    [ss_s, ss_i] = sort(sessions{k},'ascend');
    sessions{k} = ss_s;
    sessions_i{k} = ss_i;
end



N_ss = 0;           % Maximum number of session in a day
max_day = 0;            % Day with maximum number of sessions
for k = 1:length(days)
   N_ss = max(N_ss,length(sessions{k}));
   if length(sessions{k})==N_ss
       max_day = k;
   end
end

SESSION_FOLDERS     = cell(N_days,N_ss);
for k = 1:length(days)
    ss_count = 0;
    dir_content = dir(strcat(MOUSE_PATH, DAY_FOLDERS{k}));
    for kk = 1:length(dir_content)
        if(dir_content(kk).isdir) & ~strcmp(dir_content(kk).name,'.') & ~strcmp(dir_content(kk).name,'..') & ~strcmp(dir_content(kk).name,DARK_FRAMES_FOLDER)
            ss_count = ss_count + 1;
            ss_count_s = find(sessions_i{k}==ss_count);
            SESSION_FOLDERS{k,ss_count_s} = strcat(dir_content(kk).name,slsh);
        end
    end        
end


fprintf('Current mouse data structure:\n');

for k = 1:N_days
    
    fprintf('\t Day %d (%s):\n',days(k),DAY_FOLDERS{k}(1:end-1));
    for kk = 1:length(sessions{k})
        fprintf('\t \t Session %d: %s\n',kk,SESSION_FOLDERS{k,kk}(1:end-1));
    end
end





clc





% Calculate total number of sessions
tot_sessions = 0;
for curr_day = days
   curr_day_ind = find(days==curr_day);
   for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        tot_sessions = tot_sessions + 1;
   end
end




for curr_day = days
   fprintf('Day %02d:',curr_day);
   curr_day_ind = find(days==curr_day);
   for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        fprintf(' %02d -',curr_ss);
   end
   fprintf('\n');
end




%% Initialize workspace variables

% Save metadata
METADATA = cell(10);
METADATA{1} = days;
METADATA{2} = sessions;
METADATA{3} = [cell_mask_size; cell_size; min_cell_distance; rms_fact; area_min; T_cell];
METADATA{4} = [N_det; register_image_flag; remove_bad_frames];
METADATA{5} = [width; height; col_depth];
METADATA{6} = zeros(N_days,N_ss);   % cnt bad frames
METADATA{7} = zeros(N_days,N_ss);   % Identified cells
METADATA{8} = cell(N_days,N_ss);    % Time stamps
METADATA{9} = cell(N_days,N_ss);    % Time strings
METADATA{10} = eartag;              % Mouse eartag number

eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL = cell(30,1);'])
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{1} = METADATA;'])                 % Metadata
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{2} = cell(N_days,N_ss);'])        % Background subtracted registered frames
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{3} = cell(N_days,N_ss);'])        % Minimum frame
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{4} = cell(N_days,N_ss);'])        % Gaussian filtered minimum frame
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{5} = cell(N_days,N_ss);'])        % Average frame
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{6} = cell(N_days,N_ss);'])        % x shift [pixels]
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{7} = cell(N_days,N_ss);'])        % y shift [pixels]
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{8} = cell(N_days,N_ss);'])        % Dynamic cell map
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{9} = cell(N_days,N_ss);'])        % Dynamic cell map (cumulative)
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{10} = cell(N_days,N_ss);'])       % Static cell map
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{11} = cell(N_days,N_ss);'])       % Filtered cell map
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{12} = cell(N_days,N_ss);'])       % Identified cells x locations [pixels]
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{13} = cell(N_days,N_ss);'])       % Identified cells y locations [pixels]
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{14} = cell(N_days,N_ss);'])       % Identified cell areas [pixels]
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{15} = cell(N_days,N_ss);'])       % Identified cell regions after segmentation [pixel]
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{16} = cell(4,1);'])               % FULL{16} IM_REG: XSHIFT, YSHIFT, ROT, [ref_day, ref_angle]; %  Interday image registration
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{17} = cell(4,1);'])               % FULL{17} CENT_XY_DAY: CENT_X_DAY, CENT_Y_DAY, CENT_X_DAY_R, CENT_Y_DAY_R
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{18} = cell(5,1);'])               % FULL{18} CA_TRACES: CA_TRACES_RAW, CA_TRACES_FILT, CA_TRACES_BIN, CA_TRACES_OFFSET, CA_F
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{19} = cell(13,1);'])              % FULL{19} CLUSTERS: clust_idx clust_idx_s sort_idx CONN_FULL conn_thresh iterations bin_size curr_days curr_sessions metric tot_kmeans_iter max_iter global_clusters
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{20} = cell(3,1);'])               % FULL{20} GLOBAL_CELLS: CENT_X_GLOBAL, CENT_Y_GLOBAL, IDX_MAP
eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{21} = cell(5,1);'])               % FULL{21} CA_TRACES_G: CA_TRACES_RAW_G, CA_TRACES_FILT_G, CA_TRACES_BIN_G, CA_TRACES_ONSET_G, CA_F_G

   


% For rotarod, skip first and last frames
FIRST = cell(8,1);
LAST = cell(8,1);








% % Compress binary data (for old data format)
% 
% eval(['a = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{8};'])
% eval(['b = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{9};'])
% 
% eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{8} = cell(N_days,N_ss);'])
% eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{9} = cell(N_days,N_ss);'])
% 
% for curr_day = days
%    curr_day_ind = find(days==curr_day);
%    for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         
%         fprintf('Day %d Session %d\n',curr_day,curr_ss);
%         CELLS_DYN = a{curr_day_ind,curr_ss_ind};
%         CELLS_DYN_SUM = b{curr_day_ind,curr_ss_ind};
%         
%         
%         CELLS_DYN_PACKED = uint32(zeros(ceil(height/32),width,size(CELLS_DYN,3)));
%         for jj = 1:size(CELLS_DYN,3)
%             CELLS_DYN_PACKED(:,:,jj) = bwpack(CELLS_DYN(:,:,jj));
%         end
%         
%         CELLS_DYN_SUM_PACKED = uint32(zeros(ceil(height/32),width,size(CELLS_DYN_SUM,3)));
%         for jj = 1:size(CELLS_DYN_SUM,3)
%             CELLS_DYN_SUM_PACKED(:,:,jj) = bwpack(CELLS_DYN_SUM(:,:,jj));
%         end
%         
%         eval(['MOUSE_',num2str(cell_type),'_',num2str(mouse_numb),'_FULL{8}{curr_day_ind,curr_ss_ind}= CELLS_DYN_PACKED;'])
%         eval(['MOUSE_',num2str(cell_type),'_',num2str(mouse_numb),'_FULL{9}{curr_day_ind,curr_ss_ind}= CELLS_DYN_SUM_PACKED;'])
%    end
% end

     


%%  Begin batch processing

curr_ss_tot = 0;

for curr_day = days
   curr_day_ind = find(days==curr_day);
   large_filename = strcat(BASEDIR,slsh,'Frames',slsh,'img_day_',num2str(curr_day,'%02d'),'.tif');
   for curr_ss = sessions{curr_day_ind}
        tic;
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        
        curr_ss_tot = curr_ss_tot + 1;
        
        fprintf('Mouse %d D%d, D%dS%d (%d/%d). ',mouse_numb, cell_type, curr_day, curr_ss, curr_ss_tot, tot_sessions);
        
        % Set working directory
        
        DAY_FOLDER = DAY_FOLDERS{curr_day_ind};
        SESSION_FOLDER = SESSION_FOLDERS{curr_day_ind,curr_ss_ind};
        filePath = strcat(MOUSE_PATH, DAY_FOLDER, SESSION_FOLDER);  
        
        % % Dark frames defined for all sessions
        % DARK_FRAMES_DIR = strcat(MOUSE_PATH, DAY_FOLDER, SESSION_FOLDER,DARK_FRAMES_FOLDER, slsh);
        
        % % One dark frame calculated for each day
        DARK_FRAMES_DIR = strcat(MOUSE_PATH, DAY_FOLDER, DARK_FRAMES_FOLDER, slsh);
        
        
        % DARK_FRAMES_DIR = strcat(BASEDIR,slsh,DAY_FOLDER(1:end-1),'_background\');
        % DARK_FRAMES_DIR = strcat(MOUSE_PATH,DAY_FOLDER,SESSION_FOLDER,DARK_FRAMES_FOLDER,slsh);
        

        
    
        
        % Load dark frames and calculate average
        
        if dark_frames == true
            
            dir_content = dir(DARK_FRAMES_DIR);
            fileName = '';
            str_ind = 0;
            frame_ind = [];
            
            % Find the first .tif, remove the frame number to find file name
            % Frame do not need be sequential nor continuous
            for k = 1:length(dir_content)
                fileName = dir_content(k).name;
                if strcmp(fileName(max(1,length(fileName)-3):length(fileName)),'.tif') == true
                    str_ind = str_ind + 1;
                    frame_ind = [frame_ind k];
                end
            end
            
            frame_numbers = zeros(length(frame_ind),1);
            frames = cell(length(frame_ind),1);
            str_ind = 1;
            for k = frame_ind
                fileName = dir_content(k).name;
                frames{str_ind} = fileName;
                fileName = fileName(1:end-4);
                fileName2 = fileName;
                numb_l = 0;
                while strcmp(fileName2(end), '_') == 0
                    fileName2 = fileName2(1:end-1);
                    numb_l = numb_l + 1;
                end
                fileName = fileName(end-numb_l+1:end);
                frame_numbers(str_ind) = str2double(fileName);
                str_ind = str_ind + 1;
            end
            
            % Sort frames
            [sort_numb, sort_idx ] = sort(frame_numbers,'ascend');
            frames = frames(sort_idx);
            N_frames = length(frames);
            
            
            
            
            DARK_FRAMES = uint16(zeros(height, width, N_frames));
            
            if N_frames == 0
                fprintf('WARNING: no dark frames found in the specified folder. Using black mask.');
                DARK_FRAME = zeros(height,width);
            else
                
                for k = 1:N_frames
                    curr_str = strcat(DARK_FRAMES_DIR,frames{k});
                    DARK_FRAMES(:,:,k) = uint16(imread(curr_str)/64);
                end
                DARK_FRAME = uint16(mean(DARK_FRAMES,3));
                
                % figure; imagesc(DARK_FRAME); colormap(gray); caxis([0 100]);
            end
        else
            DARK_FRAME = zeros(height,width);
        end
        
        
        
        
        
        
        
        
        
        
%         first_frame = FIRST{curr_day_ind}(curr_ss_ind); 
%         last_frame = LAST{curr_day_ind}(curr_ss_ind);
%         if first_frame == 0
%             first_frame = [];
%         end
%         if last_frame == 0
%             last_frame = [];
%         end
        
        first_frame = [];
        last_frame = [];
       
        
        
        
         % Perform calculations
        
        if no_cell_id % Do not perform cell identification
            [FRAMES, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, CNT_BAD_FRAMES, X_SHIFT, Y_SHIFT] = cell_ident_min(MacOS, use_wbar, remove_bad_frames, register_image_flag, show_movie, filePath, DARK_FRAME, N_det, width, height, first_frame, last_frame);
        
        else
            
           
            [FRAMES, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, CNT_BAD_FRAMES, X_SHIFT, Y_SHIFT, IDENTIFIED_CELLS, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CELL_AREAS, CELL_REGIONS] = cell_ident(MacOS, use_wbar, remove_bad_frames, register_image_flag, show_movie, filePath, DARK_FRAME, N_det, width, height, cell_size, cell_mask_size, T_cell, area_min, min_cell_distance, rms_fact, first_frame, last_frame);
            
            % Compress binary data
            CELLS_DYN_PACKED = uint32(zeros(ceil(height/32),width,size(CELLS_DYN,3)));
            for jj = 1:size(CELLS_DYN,3)
                CELLS_DYN_PACKED(:,:,jj) = bwpack(CELLS_DYN(:,:,jj));
            end
            
            CELLS_DYN_SUM_PACKED = uint32(zeros(ceil(height/32),width,size(CELLS_DYN_SUM,3)));
            for jj = 1:size(CELLS_DYN_SUM,3)
                CELLS_DYN_SUM_PACKED(:,:,jj) = bwpack(CELLS_DYN_SUM(:,:,jj));
            end
        end
        
        
        % Load time stamps
        if load_timestamps
            [t, time_str] = get_timestamps(filePath);
            if ~isempty(first_frame) && ~isempty(last_frame)
                t = t(first_frame:last_frame);
            end
            t = t-t(1); % Start from 0
        else
            t = 0.1.*[0:size(FRAMES,3)-1]';
            time_str = '';
        end
        METADATA{8}{curr_day_ind,curr_ss_ind} = t; % Time stamps
        METADATA{9}{curr_day_ind,curr_ss_ind} = time_str;
        eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{1}{8}{curr_day_ind,curr_ss_ind} = t;'])                 % Metadata
        eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{1}{9}{curr_day_ind,curr_ss_ind} = time_str;'])

        
        if min(N_det, length(t)) ~=  size(FRAMES,3)
            error('ERROR: number of frames does not match timestamps!');
        end
        
        
        
%          % Save frames only (old code)
%         FRAMES = get_frames(MacOS, use_wbar, remove_bad_frames, register_image_flag, filePath, DARK_FRAME, N_det, width, height);
        
        
        
        if size(FRAMES,1)>0
            % Save data for current mouse
            eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{1}{6}(curr_day_ind,curr_ss_ind) = CNT_BAD_FRAMES;'])
            
            % eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{2}{curr_day_ind,curr_ss_ind} = FRAMES;'])
            eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{3}{curr_day_ind,curr_ss_ind} = MIN_FRAME;'])
            eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{4}{curr_day_ind,curr_ss_ind} = MIN_FRAME_FILT;'])
            eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{5}{curr_day_ind,curr_ss_ind} = AVG_FRAME;'])
            eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{6}{curr_day_ind,curr_ss_ind} = X_SHIFT;'])
            eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{7}{curr_day_ind,curr_ss_ind} = Y_SHIFT;'])
            
            
            if ~no_cell_id
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{1}{7}(curr_day_ind,curr_ss_ind) = IDENTIFIED_CELLS;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{8}{curr_day_ind,curr_ss_ind} = CELLS_DYN_PACKED;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{9}{curr_day_ind,curr_ss_ind} = CELLS_DYN_SUM_PACKED;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{10}{curr_day_ind,curr_ss_ind} = CELLS;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{11}{curr_day_ind,curr_ss_ind} = CELLS_FILT;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{12}{curr_day_ind,curr_ss_ind} = CENT_X;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{13}{curr_day_ind,curr_ss_ind} = CENT_Y;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{14}{curr_day_ind,curr_ss_ind} = CELL_AREAS;'])
                eval(['MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL{15}{curr_day_ind,curr_ss_ind} = CELL_REGIONS;'])
            end
            
            
            % Save frames to a separate file, one per experiment
            frames = FRAMES;
            if save_mat_frames
                save(strcat(BASEDIR,slsh,'Frames',slsh,'Frames_D',num2str(cell_type),'M',num2str(mouse_numb),'_D',num2str(curr_day),'S',num2str(curr_ss),'.mat'),'frames','-v7.3');
            end
        
            if save_tiff
               fprintf('Saving tiff...');
 
               skip_frames = 3; % skip first frame(s) of each session and replace with first frame
    
               %     curr_folder = strcat(tmpdir,'day_',num2str(curr_day,'%02d'));
               %     curr_cmd = strcat('mkdir',{' '},curr_folder);
               %     curr_cmd = curr_cmd{1};
               %     dos(curr_cmd);

               % Load frames/time vector for each experiment
      
               N_frames = size(frames,3);
               
               for k = 1:skip_frames
                   frames(:,:,k) = frames(:,:,skip_frames+1);
               end
               
               
               % fprintf('Day %d/%d, session %d/%d ...\n',curr_day,length(days),curr_ss_ind,length(sessions{curr_day_ind}));
               
               
               % % Uncomment this line for separate files per session
               % large_filename = strcat(tmpdir,'img_day_',num2str(curr_day,'%02d'),'_ss',num2str(curr_ss),'.tif');
               
               
               
               % img = (double(frames(:,:,skip_frames+ii)))/1024;
               % imwrite(img,large_filename);
               
               
               
               % curr_filename = strcat(BASEDIR,slsh,'Frames',slsh,'D',num2str(cell_type),'M',num2str(mouse_numb),'_day_',num2str(curr_day,'%02d'),'_ss_',num2str(curr_ss,'%02d'),'.tif');
               curr_filename = strcat('D:',slsh,'Giovanni',slsh,'data',slsh,'Yan',slsh,'16_Nac-1.2018',slsh,'Frames',slsh,'D',num2str(cell_type),'M',num2str(mouse_numb),'_day_',num2str(curr_day,'%02d'),'_ss_',num2str(curr_ss,'%02d'),'.tif');
               
               
               for ii = 1:N_frames
                   % curr_filename = strcat(curr_folder,'\','img_day_',num2str(curr_day,'%02d'),'_ss_',num2str(curr_ss,'%02d'),'_',num2str(ii,'%04d'),'.tif');
                   
                   % img = (double(frames(:,:,skip_frames+ii)))./1024;
                   img = (double(frames(:,:,ii)))/1024;
                   
                   %     img2 = zeros(size(frames,1),size(frames,2),3);
                   %       img2 = zeros(448,448,3);
                   %     img2(1:size(frames,1),1:size(frames,2),1) = img;
                   %     img2(1:size(frames,1),1:size(frames,2),1) = img;
                   %     img2(1:size(frames,1),1:size(frames,2),1) = img;
                   %     img2(1:size(frames,1),1:size(frames,2),1) = img;
                   %     img2(1:size(frames,1),1:size(frames,2),1) = img;
                   %     img2(1:size(frames,1),1:size(frames,2),1) = img;
                   % curr_filename = strcat(tmpdir,'img_',num2str(ii,'%04d'),'.tif');
                   % imwrite(img, curr_filename);
                   
                   
                   
                   
                   % imwrite(img,large_filename,'compression','none','WriteMode','append');
                   imwrite(img,curr_filename,'compression','none','WriteMode','append');
                   
                   
                   
                   % imwrite(img,curr_filename,'compression','none');
                   
               end
               
               % nams = cell(N_frames,1);
               % for ii = 1:N_frames
               %     curr_filename = strcat(tmpdir,'img_',num2str(ii,'%04d'),'.tif');
               %     nams{ii} = curr_filename;
               % end
               % fprintf('Conversion to .tiff done.\n');
               %
               


               
               fprintf('Done\n');
            end
        
        end  
        
        
        clear cell_ident
        all_var = who;
        for k = 1:length(all_var)
            if isempty(strfind('all_var slsh cell_type mouse_numb eartag days sessions curr_day curr_ss curr_day_ind curr_ss_ind curr_ss_tot tot_sessions filePath FIRST LAST METADATA BASEDIR DARK_FRAMES_FOLDER DAY_FOLDERS SESSION_FOLDERS N_det MacOS use_wbar remove_bad_frames dark_frames register_image_flag load_timestamps show_movie height width  col_depth   cell_size cell_mask_size T_cell area_min min_cell_distance rms_fact no_cell_id save_mat_frames save_tiff',all_var{k})) && ~strcmp('MOUSE_',all_var{k}(1:min(length(all_var{k}),6)))
                eval(['clear ',all_var{k}]);
            end
        end
        
        clear all_var

        % clear DARK_FRAMES FRAMES frames MIN_FRAME MIN_FRAME_FILT AVG_FRAME  CNT_BAD_FRAMES  X_SHIFT  Y_SHIFT  IDENTIFIED_CELLS  CELLS_DYN  CELLS_DYN_SUM  CELLS  CELLS_FILT  CENT_X  CENT_Y  CELL_AREAS  CELL_REGIONS
        
        % clear -except MOUSE_* slsh cell_type mouse_numb eartag BASEDIR DARK_FRAMES_FOLDER DAY_FOLDERS SESSION_FOLDERS N_det MacOS use_wbar remove_bad_frames dark_frames register_image_flag load_timestamps show_movie height width  col_depth   cell_size cell_mask_size T_cell area_min min_cell_distance rms_fact
        
        curr_time = round(toc);
        curr_min = floor(curr_time/60);
        curr_s = rem(curr_time,60);
        
        fprintf('Elapsed time is %d min and %d s.\n', curr_min, curr_s)
   end
end



%% Clear unused variables
all_var = who;
for k = 1:length(all_var)
    if isempty(strfind('all_var slsh cell_type mouse_numb eartag BASEDIR',all_var{k})) && ~strcmp('MOUSE_',all_var{k}(1:min(length(all_var{k}),6)))
        eval(['clear ',all_var{k}]);
    end
end

clear all_var k 


%% Save full workspace
fprintf('Saving workspace to disk...\n');
eval(['save(''D:' slsh,'Giovanni',slsh,'data',slsh,'D',num2str(cell_type),'M',num2str(mouse_numb),'.mat'',''-v7.3'');'])
fprintf('Done.\n');



%% Show min frames


eval(['MOUSE = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL;'])
[FRAMES, N_days, N_ss, T, width, height, days, sessions, cell_mask_size, cell_size, min_cell_distance, rms_fact, area_min, T_cell, N_det, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, X_SHIFT, Y_SHIFT, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CENT_X_DAY, CENT_Y_DAY, CA_TRACES_RAW, CA_TRACES_FILT, CA_TRACES_BIN, CA_TRACES_ONSET clust_idx clust_idx_s clust_cell_s  CLUST_CONN_FULL clust_conn_thresh clust_iterations clust_bin_size clust_curr_days clust_curr_sessions clust_metric clust_tot_kmeans_iter clust_max_iter global_clusters] = load_mouse_full(MOUSE);


% % One figure for all
% 
% N_col = 16;
% N_rows = ceil((N_days*N_ss)/N_col);
% figure; hold on;
% 
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         subplot(N_rows, N_col, N_ss*(curr_day_ind-1)+curr_ss_ind);
%         
%         imagesc(MIN_FRAME{curr_day_ind, curr_ss_ind});
%         colormap(gray); 
%         set(gca,'Ydir','reverse');
%         axis equal
%         xlim([0, MOUSE{1}{5}(1)]);
%         ylim([0, MOUSE{1}{5}(2)]);
%         set(gca,'xtick',[]);
%         set(gca,'ytick',[]);
%         title(strcat('D',num2str(curr_day),'S',num2str(curr_ss)));
%        
%     end
% end

% Or one figure per day

N_col = 5;
N_rows = ceil((N_ss)/N_col);


for curr_day0 = 1:length(days)
    curr_day = days(length(days)-curr_day0+1);
    curr_day_ind = find(days==curr_day);
    curr_fig = figure('Position',[30,30,1600,1200]); hold on;
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        subplot(N_rows, N_col, curr_ss);
        
        imagesc(MIN_FRAME{curr_day_ind, curr_ss_ind});
        colormap(gray); 
        set(gca,'Ydir','reverse');
        axis equal
        xlim([0, MOUSE{1}{5}(1)]);
        ylim([0, MOUSE{1}{5}(2)]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        title(strcat('D',num2str(curr_day),'S',num2str(curr_ss)));
       
    end
end



clear MOUSE;














%% Plot x and y shift


eval(['MOUSE = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL;'])
[FRAMES, N_days, N_ss, T, width, height, days, sessions, cell_mask_size, cell_size, min_cell_distance, rms_fact, ~, T_cell, N_det, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, X_SHIFT, Y_SHIFT, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CENT_X_DAY, CENT_Y_DAY, CA_TRACES_RAW, CA_TRACES_FILT, CA_TRACES_BIN, CA_TRACES_ONSET clust_idx clust_idx_s clust_cell_s  CLUST_CONN_FULL clust_conn_thresh clust_iterations clust_bin_size clust_curr_days clust_curr_sessions clust_metric clust_tot_kmeans_iter clust_max_iter global_clusters] = load_mouse_full(MOUSE);





% figure; hold on;
% 
% N_col = 8;
% N_rows = ceil((N_days*N_ss)/N_col);
% 
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     for curr_ss = sessions{curr_day_ind}
%         
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         subplot(N_rows, N_col, N_ss*(curr_day_ind-1)+curr_ss_ind);
%         plot(X_SHIFT{curr_day_ind,curr_ss_ind}); hold on;
%         plot(Y_SHIFT{curr_day_ind,curr_ss_ind}); hold on;
%         ylim([-6 6]);
%         % title('X and Y sahift [pixels]');
%         title(strcat('D',num2str(curr_day),'S',num2str(curr_ss)));
%         set(gca,'xtick',[]);
%         % set(gca,'ytick',[]);
%     end
% end



% Or one figure per day

N_col = 5;
N_rows = ceil((N_ss)/N_col);

for curr_day0 = 1:length(days)
    curr_day = days(length(days)-curr_day0+1);
    curr_day_ind = find(days==curr_day);
    curr_fig = figure('Position',[30,30,1600,1200]); hold on;
    for curr_ss = sessions{curr_day_ind}
        
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        subplot(N_rows, N_col, curr_ss);
        plot(X_SHIFT{curr_day_ind,curr_ss_ind}); hold on;
        plot(Y_SHIFT{curr_day_ind,curr_ss_ind}); hold on;
        ylim([-6 6]);
        % title('X and Y sahift [pixels]');
        title(strcat('D',num2str(curr_day),'S',num2str(curr_ss)));
        set(gca,'xtick',[]);
        % set(gca,'ytick',[]);
    end
end









clear MOUSE;





%% Show all cells masks

eval(['MOUSE = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL;'])
[FRAMES, N_days, N_ss, T, width, height, days, sessions, cell_mask_size, cell_size, min_cell_distance, rms_fact, area_min, T_cell, N_det, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, X_SHIFT, Y_SHIFT, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CENT_X_DAY, CENT_Y_DAY, CA_TRACES_RAW, CA_TRACES_FILT, CA_TRACES_BIN, CA_TRACES_ONSET clust_idx clust_idx_s clust_cell_s  CLUST_CONN_FULL clust_conn_thresh clust_iterations clust_bin_size clust_curr_days clust_curr_sessions clust_metric clust_tot_kmeans_iter clust_max_iter global_clusters] = load_mouse_full(MOUSE);




% figure; hold on;
% 
% N_col = 8;
% N_rows = ceil((N_days*N_ss)/N_col);
% 
% for curr_day = days
%     curr_day_ind = find(days==curr_day);
%     for curr_ss = sessions{curr_day_ind}
%         curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
%         subplot(N_rows, N_col, N_ss*(curr_day_ind-1)+curr_ss_ind);
%         imagesc(MOUSE{10}{curr_day_ind, curr_ss_ind});
%         colormap(gray); 
%         set(gca,'Ydir','reverse');
%         axis equal
%         xlim([0, MOUSE{1}{5}(1)]);
%         ylim([0, MOUSE{1}{5}(1)]);
%         set(gca,'xtick',[]);
%         set(gca,'ytick',[]);
%         title(strcat('D',num2str(curr_day),'S',num2str(curr_ss)));
%     end
% end





% Or one figure per day

N_col = 5;
N_rows = ceil((N_ss)/N_col);

for curr_day0 = 1:length(days)
    curr_day = days(length(days)-curr_day0+1);
    curr_day_ind = find(days==curr_day);
    curr_fig = figure('Position',[30,30,1600,1200]); hold on;
    for curr_ss = sessions{curr_day_ind}
        curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
        subplot(N_rows, N_col, curr_ss);
        imagesc(MOUSE{10}{curr_day_ind, curr_ss_ind});
        colormap(gray); 
        set(gca,'Ydir','reverse');
        axis equal
        xlim([0, MOUSE{1}{5}(1)]);
        ylim([0, MOUSE{1}{5}(1)]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        title(strcat('D',num2str(curr_day),'S',num2str(curr_ss)));
    end
end

        
        
        
clear MOUSE;



%% Show cells mask for specific day/session...D736 days = [1,2,22,23,29,41]

eval(['MOUSE = MOUSE_', num2str(cell_type),'_',num2str(mouse_numb),'_FULL;'])
curr_day = 1;
curr_ss = 1;
curr_day_ind = find(days==curr_day);
curr_ss_ind = find(sessions{curr_day_ind}==curr_ss);
figure; hold on;
imagesc(MOUSE{10}{curr_day_ind, curr_ss_ind});
colormap(gray); 
set(gca,'Ydir','reverse');
axis equal
xlim([0, MOUSE{1}{5}(1)]);
xlim([0, MOUSE{1}{5}(1)]);


% figure; imagesc(MOUSE{9}{1,1} + MOUSE{9}{1,2} + MOUSE{9}{1,3} ); colormap(gray);


clear MOUSE;






