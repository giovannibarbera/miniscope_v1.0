function [FRAMES, MIN_FRAME, MIN_FRAME_FILT,  AVG_FRAME, CNT_BAD_FRAMES, X_SHIFT, Y_SHIFT, IDENTIFIED_CELLS, CELLS_DYN, CELLS_DYN_SUM, CELLS, CELLS_FILT, CENT_X, CENT_Y, CELL_AREAS, CELL_REGIONS] = cell_ident(MacOS, use_wbar, remove_bad_frames, register_image_flag, show_movie, filePath, DARK_FRAME, N_det, width, height, cell_size, cell_mask_size, T_cell, area_min, min_cell_distance, rms_fact, first, last)
%CELL_IDENT Load raw images and apply image registration and automatic cell
%detection algorithm
%
% MacOS:                    (bool) 1 for OSX, 0 for Windows        
% use_wbar                  (bool) use status bar
% remove_bad_frames         (bool) remove and mark bad frames
% register_image_flag       (bool) perform image registration
% show_movie                (bool) display movie at the end of processing
% filePath                  (string array) full path
% N_det                     (scalar) numer of frames used for cell detection
% width                     (scalar) image width
% height                    (scalar) image height
% cell_size                 (scalar) size of the buffer between positive and negative peak in the gradient
% cell_mask_size            (scalar) x and y cell mask size
% T_cell                    (scalar) Cell persistence parameter 
% area_min                  (scalar) cell area filtering (may not reflect actual cell area)
% min_cell_distance         (scalar) minimum euclidean distance of centroids to be recognized as cells   
% rms_fact                  (scalar) rms factor threshold for gradient peaks identification 
% FRAMES                    (height,width,N_frames) frame data, registered and background subtracted
% MIN_FRAME                 (height,width) minimum pixel value
% MIN_FRAME_FILT            (height,width) minimum pixel value (used for background subtraction)
% AVG_FRAME                 (height, width) average pixel value
% CNT_BAD_FRAMES            (scalar) number of bad frames detected and replaced
% X_SHIFT                   (1,N_frames) x shift of frames before registration
% Y_SHIFT                   (1,N_frames) y shift of frames before registration
% IDENTIFIED_CELLS          (scalar) number of identified cells 
% CELLS_DYN                 (height,width,N_frames) dynamic incremental cell mask
% CELLS                     (height,width) cell mask for all considered frames N_det
% CELLS_FILT                (height,widht) filtered cell mask
% CENT_X                    (1,IDENTIFIED_CELLS) x coordinates [pixels] of detected cells
% CENT_Y                    (1,IDENTIFIED_CELLS) y coordinates [pixels] of detected cells
%
%   giobarbera@neuralmapper.com


% To batch rename files in terminal: for f in *.tif; do echo mv "$f" "${f/_4_/_6_}"; done


FRAMES = []; 
MIN_FRAME = []; 
MIN_FRAME_FILT= [];  
AVG_FRAME = [];  
CNT_BAD_FRAMES = []; 
X_SHIFT = []; 
Y_SHIFT = []; 
IDENTIFIED_CELLS = []; 
CELLS_DYN = []; 
CELLS_DYN_SUM = []; 
CELLS = []; 
CELLS_FILT = []; 
CENT_X = []; 
CENT_Y = []; 
CELL_AREAS = []; 
CELL_REGIONS = [];


% (m, m) binary mask to associate with each cell centroid

CELL_MASK = [0 1 1 0;
             1 1 1 1;
             1 1 1 1;
             0 1 1 0];  

CELL_MASK_SMALL = [1 1;
    1 1];


% discrete gaussian fiter (FRAMES)
GAUSS = (gauss_filt(7,1));

% continuos gaussiang filter (MIN FRAME)
GAUSS_C = (gauss_filt(5,1));





tic







%% Load frames from current folder

dir_content = dir(filePath);
fileName = '';
str_ind = 0;
frame_ind = [];



% Find the first .tif, remove the frame number to find file name
% Frames do not need be sequential nor continuous
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
[sort_numb sort_idx ] = sort(frame_numbers,'ascend');

frames = frames(sort_idx);

N_frames = length(frames);


if N_frames == 0
    error('ERROR: the specified directory does not contain any tiff image.');
    % fprintf('No frames were found in %s, skipping.\n', filePath);
    return
end






if N_frames <= N_det
    % fprintf('WARNING: number of frames is less than N_det. N_det set to %d.\n',min(N_det, N_frames));
    N_det = min(N_det, N_frames);
   
else
    FRAMES(:,:,N_det:end) = [];
    N_frames = N_det;
end








 
if use_wbar
    wbar = waitbar(0,'Loading frames');
end

fprintf('Loading %d frames from %s. ',N_frames, filePath);


FRAMES = uint16(zeros(height, width, N_frames));


for k = 1:N_frames
    waitbar(k/N_frames,wbar,'Loading frames')
    curr_str = strcat(filePath,frames{k});

  
    
    % % For NeuralMapper data
    FRAMES(:,:,k) = uint16(imread(curr_str)/64);
    
% %     For simulated data
%     FRAMES(:,:,k) = uint16(imread(curr_str));
end



% [FRAMES, N_frames] = load_data2(filePath, frames, height, width);
% [FRAMES, N_frames] = load_data('../data/D2/Mouse6/Image5/Ijco3/', 'img_data_12_', height, width);
% [FRAMES, N_frames] = load_data('../data/D2/Mouse5/Image5/Ijco3/', 'img_data_6_', height, width);

% [FRAMES, N_frames] = load_data('/Users/gio/Desktop/3/', 'img_data_6_', height, width);


% Number of frames used for cell detection
% avoid errors if less than 3000 frames (still process what's left)




% FRAMES = double(FRAMES);









% Use first to last frames only
if ~isempty(first) && ~isempty(last)
    FRAMES = FRAMES(:,:,first:last);
    N_frames = last-first+1;
    N_det = N_frames;
end
    
    













FILT_IMG = zeros(height, width, N_det,'uint16');
GRAD_X = zeros(height, width, N_det,'int16');
GRAD_Y = zeros(height, width, N_det,'int16');

% GRAD_XT = zeros(height, width, N_det,'int8');
% GRAD_YT = zeros(height, width, N_det,'int8');
CONV_X = false(height, width, N_det);
CONV_Y = false(height, width, N_det);
CELLS = false(height, width, N_det);
CELLS_FILT = false(height, width);







%% Find bad frames and replace with previous one

cnt_bad_frames = 0;
bad_frame = 0;
   
if remove_bad_frames
    tmp_pix = 0;
    for k = 1:N_frames
        waitbar(k/N_frames,wbar,'Removing bad frames')
        bad_frame = 0;
        for m = 1:height
           for n = 1:40:width-10
               if FRAMES(m,n,k) == tmp_pix
                  indx = n;
                  while ((indx <= (n + 10)) && (FRAMES(m,indx,k) == tmp_pix ))
                      indx = indx+1;
                  end
                  if indx-n>10
                      bad_frame = 1;
                  end
               end
               tmp_pix = FRAMES(m,n,k);
           end
        end

        if bad_frame == 1
            cnt_bad_frames = cnt_bad_frames + 1;
            FRAMES(:,:,k) = FRAMES(:,:,k-1);
        end
    end
end






%% Remove dark frame


if sum(sum(DARK_FRAME)) > 0 
    for k = 1:size(FRAMES,3)
    
        if use_wbar == 1
           waitbar(k/size(FRAMES,3),wbar,'Removing dark frame')
        end

        FRAMES(:,:,k) = FRAMES(:,:,k)-DARK_FRAME;
    
    end
    
else
    % Apply median filter (for registration only) to remove salt and pepper noise when no dark frame is available
    FRAMES_FILT = FRAMES;

    for k = 1:size(FRAMES,3)
        if use_wbar == 1
           waitbar(k/size(FRAMES,3),wbar,'Applying median filter')
        end
        FRAMES_FILT(:,:,k) = medfilt2(FRAMES(:,:,k));
    end

end











%% Calculate min frame

 MIN_FRAME = uint16(zeros(height,width));
 
% % Use most frequent value for very large datasets (mostly dark, N_frames>3000)
% 
% for m = 1:width
%     if use_wbar == 1
%         waitbar(m/width,wbar,'Calculating background for subtraction')
%     end
%     for n = 1:height
%         [tmp1, tmp2] = hist(squeeze(double(FRAMES(n,m,:))),100);
%         [~, tmp4] = max(tmp1);
%         MIN_FRAME(n,m) = tmp2(tmp4);
%     end
% end









% Use min pix value (if there are no bad frames)

skip_frames = min(30,size(FRAMES,3));
% skip_frames = min(170,size(FRAMES,3));
% skip_frames = min(1,size(FRAMES,3));

% Exclude frames with exstreme intensity changes (dark LED)
FRAMES_2 = double(FRAMES);
frame_idx = zeros(1,size(FRAMES,3));
mn = mean(mean(mean(FRAMES,1),2),3);
for k = 1:size(FRAMES,3)
    frame_idx(k) = mean(mean(FRAMES(:,:,k),1),2);
end
frame_mean = frame_idx;
% Remove dark frames
frame_idx1 = find(frame_idx < mn-25);
frame_idx = [diff(frame_idx) max(frame_idx)];
% Remove frames with large intensity changes
frame_idx2 = find(abs(frame_idx)>20);
frame_idx = union(frame_idx1, frame_idx2);
if ~isempty(frame_idx([frame_idx<size(FRAMES,3)] & [frame_idx>min(skip_frames,size(FRAMES,3))]))
    fprintf('WARNING: large intensity changes detected, possible LED issue.');
    figure; plot(frame_mean); ylabel('Average pixel intensity'); xlabel('Frame number');
    hold on;
    frame_err = nan(1,size(FRAMES,3));
    frame_err(frame_idx([frame_idx<size(FRAMES,3)] & [frame_idx>min(skip_frames,size(FRAMES,3))])) = frame_mean(frame_idx([frame_idx<size(FRAMES,3)] & [frame_idx>min(skip_frames,size(FRAMES,3))]));
    plot(frame_err,'r','LineWidth',2);
    for k = frame_idx
        FRAMES_2(:,:,k) = nan(size(FRAMES,1),size(FRAMES,2));
    end
end
for m = 1:width
    if use_wbar == 1
        waitbar(m/width,wbar,'Calculating background for subtraction')
    end
    for n = 1:height
        %         [tmp1, tmp2] = hist(squeeze(double(FRAMES(n,m,:))),128);
        %         [~, tmp4] = max(tmp1);
        
        % MIN_FRAME(n,m) = nanmin(FRAMES_2(n,m,min(skip_frames,size(FRAMES,3)-1):max(min(skip_frames,size(FRAMES,3)-1),size(FRAMES,3)-5)));
        MIN_FRAME(n,m) = nanmin(FRAMES(n,m,min(skip_frames,size(FRAMES,3)-1):max(min(skip_frames,size(FRAMES,3)-1),size(FRAMES,3)-5)));
        
        
    end
end
MIN_FRAME = uint16(MIN_FRAME);



if size(FRAMES,3)<20
    MIN_FRAME = uint16(zeros(height,width));
end

    
    


%% Perform image registration

% Make a copy just for debugging

if sum(sum(DARK_FRAME)) > 0
    FRAMES2 = FRAMES;
else
    FRAMES2 = FRAMES_FILT;
end




x_shift = zeros(1, N_frames);
y_shift = zeros(1, N_frames); 
x_shift_1 = zeros(1, N_frames);
y_shift_1 = zeros(1, N_frames); 

if register_image_flag

    b_width = round(1.6*width);
    b_height = round(1.6*height);
    blackman_filt_w = blackman(b_width,'symmetric');
    blackman_filt_h = blackman(b_height,'symmetric');
    ff=(blackman_filt_w*blackman_filt_h'); 

    ff = ff( floor((b_height-height)/2):floor((b_height-height)/2) + height-1, floor((b_width-width)/2):floor((b_width-width)/2) + width-1);
    


    [filt_b,filt_a] = butter(16,0.2,'low');


    NFFTX = width;%1024; %2^nextpow2(width); % Next power of 2 from width
    freqx = 1/2*linspace(0,1,NFFTX/2+1);
    NFFTY = height; %1024; %2^nextpow2(width); % Next power of 2 from width
    freqy = 1/2*linspace(0,1,NFFTY/2+1);

    corr = zeros(NFFTY, NFFTX, N_frames);

    int_min = 199; 
    int_max = 203; 
    int_step = 0.03;

    shift = 0;

    shiftx = zeros(1,N_frames);
    shifty = zeros(1,N_frames);

    for k = 1:N_frames
        shiftx(k) = 10.*(sin(2*pi*k*0.001)+sin(2*pi*k*0.0033)+sin(2*pi*k*0.001)+sin(2*pi*k*0.01)+sin(2*pi*k*0.0001));
        shifty(k) = 10.*(sin(2*pi*k*0.003)+sin(2*pi*k*0.002)+sin(2*pi*k*0.0066)+sin(2*pi*k*0.04)+sin(2*pi*k*0.0002));
    end
    shiftx = round(abs(shiftx));
    shifty = round(abs(shifty)); %abs need to work on negative shift code


    x_shift_f = zeros(1, N_frames);
    y_shift_f = zeros(1, N_frames);

    x_shift_1 = zeros(1, N_frames);
    y_shift_1 = zeros(1, N_frames);

  

     low_cut = 3;
     high_cut = 50;


    FFTF = ones(height, width);
    for m = 1:height
        for n = 1:width
            if (m <= low_cut || m > width-low_cut || n <= low_cut || n > height-low_cut) 
                FFTF(m,n) = 0;
            end

            if ( m > high_cut && m < width - high_cut)
                FFTF(m,n) = 0;
            end

            if ( n > high_cut && n < height - high_cut)
                FFTF(m,n) = 0;
            end

        end
    end


    img1 = double(MIN_FRAME);
    % img1 = double(AVG_FRAME);


    img1 = img1.*ff;
    F1 = fft2(img1,NFFTY,NFFTX);

    
    F1 = F1 .* FFTF;
    
    
    R = zeros(NFFTX, NFFTY);


    % %

    % FRAMES_SHIFT = zeros(height,width,N_frames);
    % FRAMES_SHIFT2 = zeros(height,width,N_frames);

    for k = 2:N_frames
        if use_wbar == 1
            waitbar(k/N_frames,wbar,'Performing image registration')
        end

        % Define current pair of images to register
        % img1 = double(FRAMES2(:,:,k-1));
        % img1 = double(FRAMES2(:,:,1875));

%         img1 = double(FRAMES2(:,:,k-1)+MIN_FRAME);
%         img1 = img1.*ff;
%         F1 = fft2(double(img1),NFFTX,NFFTY);
%         F1 = F1 .* FFTF;


     % img2 = zeros(height,width);  img2(1+shift:end,1+shift:end) = FRAMES2(1:end-shift,1:end-shift,k)+MIN_FRAME(1:end-shift,1:end-shift);
     img2 = zeros(height,width);  img2(1+shift:end,1+shift:end) = double(FRAMES2(1:end-shift,1:end-shift,k));


        % img2 = double(FRAMES2(:,:,k)+MIN_FRAME);
        % Test on artificially shifted image


%         % Artificially shift image
%         img2 = zeros(height,width);  img2(1+shiftx(k):end,1+shifty(k):end) = FRAMES2(1:end-shiftx(k),1:end-shifty(k),k)+MIN_FRAME(1:end-shiftx(k),1:end-shifty(k));
%         img2b = img2;
%         img2b(1+shiftx(k):end,1+shifty(k):end) = FRAMES2(1:end-shiftx(k),1:end-shifty(k),k);
%         FRAMES_SHIFT(:,:,k) = img2;
%         FRAMES_SHIFT2(:,:,k) =img2b;



        

        % Normalize images
%         img1 = img1 - min(min(img1));%-50;
%         img1(img1<0) = 0;
%         img1 = img1/max(max(img1));
% 
%         img2 = img2 - min(min(img2));%-50;
%         img2(img2<0) = 0;
%         img2 = img2/max(max(img2));



        %     % whiten images
        %     fudgefactor = 1;
        %     
        %     X = bsxfun(@minus, img1, mean(img1));
        %     A = X'*X;
        %     [V,D] = eig(A);
        %     img1 = X*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
        %     
        %     X = bsxfun(@minus, img2, mean(img2));
        %     A = X'*X;
        %     [V,D] = eig(A);
        %     img2 = X*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';


        % img1 = img1.*ff;
        img2 = img2.*ff;

%         H = padarray(2,[2 2]) - fspecial('gaussian' ,[5 5],2); % create unsharp mask
%         img1 = imfilter(img1,H);  % create a sharpened version of the image using that mask
%         img2 = imfilter(img2,H);  % create a sharpened version of the image using that mask


        % Calculate fft
        % F1 = fft2(double(img1),NFFTX,NFFTY);
        F2 = fft2(img2,NFFTY,NFFTX);

        % F1F = hipass_filter(size(img1, 1),size(img1,2)).*abs(F1);  
        % F2F = hipass_filter(size(img2, 1),size(img2,2)).*abs(F2); 

        % [deltaX(k) deltaY(k)] = RegisterFourierMellin(img1, img2);

        % figure; imagesc(abs((ifft2(F1))));

        % Filter fft frequencies
        % F1 = F1 .* FFTF;
        F2 = F2 .* FFTF;

        % Calculate cross correlation
        R = F1.*conj(F2);
        R = R./((abs(R))+0.00000000000001);
        % R = R./((abs(R)));


        
        corr(:,:,k) = abs(fftshift((ifft2(R))));
        r = corr(:,:,k);
        [max_val_i, max_ind_i] = max(r(:));
        [max_x_i, max_y_i] = ind2sub(size(r),max_ind_i);


        % corr(max_y_i,max_x_i,k) = mean([corr(max_y_i-1:max_y_i+1,max_x_i-1,k); corr(max_y_i-1:max_y_i+1,max_x_i+1,k); corr(max_y_i-1,max_x_i,k);corr(max_y_i+1,max_x_i,k)]);
        % corr(max_y_i,max_x_i,k) = max([corr(max_y_i-1:max_y_i+1,max_x_i-1,k); corr(max_y_i-1:max_y_i+1,max_x_i+1,k); corr(max_y_i-1,max_x_i,k);corr(max_y_i+1,max_x_i,k)]);



%         % Interpolate surface
%         [X1,X2] = ndgrid((int_min:int_max));
%         [Xq1,Xq2] = ndgrid((int_min:int_step:int_max));
% 
%         Vq = interpn(X1,X2,corr(int_min:int_max,int_min:int_max,k),Xq1,Xq2,'cubic');
% 
%          %  figure; mesh([int_min:int_step:int_max],[int_min:int_step:int_max],Vq); 
% 
% 
%         [max_val,max_ind] = max(Vq(:));
%         [max_y, max_x] = ind2sub(size(Vq),max_ind);
% 
% 
%         if  k>1
%             x_shift_f(k) = x_shift_f(k-1) + (Xq1(max_x,1)-round(NFFTX/2)-1);
%             y_shift_f(k) = y_shift_f(k-1) + (Xq1(max_y,1)-round(NFFTY/2)-1);
%         end

          x_shift_1(k) = max_x_i-round(width/2) - 1;
          y_shift_1(k) = max_y_i-round(height/2) - 1;


    end


    %  x_shift_1 = filtfilt(filt_b,filt_a,x_shift_f);
    %  y_shift_1 = filtfilt(filt_b,filt_a,y_shift_f);

    x_shift = round(x_shift_1);
    y_shift = round(y_shift_1);
    
    
    
    
    
    
    
    
    % Check maximum shift
    
    if( ((max(abs(x_shift)))>5) || ((max(abs(y_shift)))>5) )
        fprintf('WARNING: shift greater than 5 pixels detected (replaced with value from previous frame).');
        figure; subplot(2,1,1); plot(x_shift); ylabel('X shift [pixel]')
        subplot(2,1,2); plot(y_shift); ylabel('Y shift [pixel]')
    end

    
    
    for k = 2:length(x_shift)
        if abs(x_shift(k)-x_shift(k-1))>=5
            x_shift(k) = x_shift(k-1);
        end
        
        if abs(y_shift(k)-y_shift(k-1))>=5
            y_shift(k) = y_shift(k-1);
        end
    end




    % MIN_FRAME2 = MIN_FRAME;
    % MIN_FRAME_FILT = imfilter(MIN_FRAME2, GAUSS_C, 'same');
   

    % figure; imagesc(MIN_FRAME); axis equal; xlim([1,width]); ylim([1, height]); colormap(gray(256)); caxis([0,col_depth]);

    
    
    % Remove min frame
    
    for k = 1:N_frames
        if use_wbar == 1 
            waitbar(k/N_frames,wbar,'Removing min frame')
        end
        FRAMES(:,:,k) = FRAMES(:,:,k) - MIN_FRAME;
        % FRAMES(:,:,k) = FRAMES(:,:,k) - AVG_FRAME;
        % FRAMES(:,:,k) = FRAMES(:,:,k) - MIN_FRAME_FILT;
    end












    % Shift image

    % Positive x_shift->left, positive y_shift->up




%    for k = 1:N_frames
%        % FRAMES2(:,:,k) = FRAMES(:,:,k) + MIN_FRAME; % Regular frames
%        FRAMES2(:,:,k) = FRAMES(:,:,k); % Regular frames 
%        %  FRAMES2(:,:,k) = FRAMES_SHIFT(:,:,k); % Artifiacilly shifted frames(with background)
%        % FRAMES2(:,:,k) = FRAMES_SHIFT2(:,:,k);  % Artifiacilly shifted frames(no background)
%    end
   
   
   % % Marked frames with lines
   % FRAMES4 = FRAMES;

    for k = 1:N_frames
        
        if use_wbar == 1
            waitbar(k/N_frames,wbar,'Shifting images')
        end

        xshift = -y_shift(k);
        yshift = -x_shift(k);
        
        
        
        curr_frame = FRAMES(:,:,k);


        % % Mark measured displacement with white lines
        % FRAMES4(:,214+xshift,k) = 1023.*ones(width,1);
        % FRAMES4(160+yshift,:,k) = 1023.*ones(1,height);



        if (xshift > 0)
            FRAMES(:,1:end-xshift,k) = curr_frame(:,xshift+1:end);
        elseif (xshift < 0)
            FRAMES(:,-xshift+1:end,k) = curr_frame(:,1:end+xshift);
        end

        curr_frame = FRAMES(:,:,k);

        if (yshift > 0)
            FRAMES(1:end-yshift,:,k) = curr_frame(yshift+1:end,:);
        elseif (yshift < 0)
            FRAMES(-yshift+1:end,:,k) = curr_frame(1:end+yshift,:);
        end

    end

    clear FRAMES4
end







%     for k = 1:N_frames-1
%         if use_wbar == 1
%             waitbar(k/N_det,wbar,'Shifting images')
%         end
%     % for k = 1535:1635
%         x_shift_curr = -x_shift(k);
%         y_shift_curr = -y_shift(k);
% 
%         % Shift x and y:
%         centered_frame = FRAMES(max(y_shift_curr,1):min(y_shift_curr+height,height),max(x_shift_curr,1):min(x_shift_curr+width,width),k); 
%         sz_cent = size(centered_frame);
%         FRAMES2(max(-y_shift_curr,1):sz_cent(1)+max(-y_shift_curr,1)-1,max(-x_shift_curr,1):sz_cent(2)+max(-x_shift_curr,1)-1,k) = centered_frame;
% 
%     %     % Mark measured displacement with black lines
%     %     FRAMES2(:,217+x_shift_curr,k) = zeros(width,1);
%     %     FRAMES2(218+y_shift_curr,:,k) = zeros(1,height);
%     end
%     
%     FRAMES = FRAMES2;



% %% Play video
% 
% save_movie = 1;
% cmin = 0;
% cmax = 400;
% fps_r=15;
% skip=1;
% xlimits = [50,350];
% ylimits = [50,350];
% 
% % xlimits = [1,max_x_shift*2+1];
% % ylimits = [1,max_y_shift*2+1];
% 
% winsize = [ 0 0 2000 1100];
% figure('Position',winsize);
% set(gca,'NextPlot','replaceChildren');
% %set(gcf,'Colormap',greencmap); % caxis([200,500])
% 
% % subplot(1, 2, 1)
% set(gcf, 'color', [0 0 0])
% 
% 
% text_position = [100 100]; % [x y]
% box_color = {'blue'};
% 
% 
% % plot green cell areas
% szfrm = size(FRAMES);
% green = cat(3, zeros(szfrm(1), szfrm(2)), ones(szfrm(1), szfrm(2)), zeros(szfrm(1), szfrm(2)));
% blue = cat(3, zeros(szfrm(1), szfrm(2)), zeros(szfrm(1), szfrm(2)), ones(szfrm(1), szfrm(2)));
% 
% hold on 
% % h = imshow(green);
% % hh = imshow(blue);
% %     set(h, 'AlphaData', 0.2*CELLS_FILT);
% 
% % % plot centroids
% % plot(CENT_X,CENT_Y,'go');
% 
% 
% if save_movie == 1
%     vid = VideoWriter('movie2.avi');
%     vid.Quality = 100;
%     vid.FrameRate = fps_r;
%     open(vid);
% end
% 
% 
% 
% 
% 
% for jj = 500:1000%3000% round(size(FRAMES,3)/skip)-1
% % for jj = 1570:1585
% % for jj = 1:9
%     
% %        subplot(1, 2, 1)
% %     % imagesc(FRAMES_FILT_BIN(:,:,jj*skip));
% %     imagesc(ROW_F(:,:,jj*skip)); 
% %     colormap(gray); axis equal; 
% %     xlim([200-20, 200+20]); ylim([ymin, ymax]);
% %     % xlim(xlimits); ylim(ylimits);
%     
%     subplot(1, 2, 1);
%     
%    
%     
% %      text_str = ['T = ' num2str(jj*skip/10,'%1.1f')];
% %     RGB = insertText(FRAMES4(:,:,jj*skip)+MIN_FRAME, text_position, text_str, 'TextColor' ,'red', 'FontSize', 20, ...
% %         'BoxColor', box_color, 'BoxOpacity', 0);
% % 
% %     imagesc(RGB(:,:,1));  colormap(gray); caxis([cmin, cmax]); hold on;
% %     % caxis([0,col_depth]); hold on;
%     
%    imagesc(FRAMES4(:,:,jj*skip));  hold on;
%    
%    
%     text_str = strcat('T=', num2str(jj*skip/10,'%1.1f'));
%    
%     text(60,67, text_str,'Color' ,'green', 'FontSize', 48);
%     
%     hold off;
%    axis equal; 
%    xlim(xlimits);
%    ylim(ylimits);
%    caxis([cmin cmax]);
%   
%    colormap(gray); 
% 
%     
%     
%     
% %     plot(rx_mean_f(:,jj*skip)) % hold on
% %     xlim([197,207]); ylim([-0.02,0.25]);
% 
% 
% 
%     % ylim([-0.04, 0.04]); % xlim([200-20,200+20]);
%     % plot(rx_mean(:,1),'r') plot(rx_mean(:,11),'g')
% 
% %     subplot(1, 3, 1)
% %     imagesc(ROW_F(:,:,jj*skip)); colormap(gray) % hold on
% %     ylim([180, 250]); xlim([258-20,258+20]);
% %     % plot(rx_mean(:,1),'r') plot(rx_mean(:,11),'g')
% 
% 
% 
% 
%     
%     
%     subplot(1,2,2)
%     % imagesc(FRAMES(:,:,jj*skip));  
%     imagesc(FRAMES3(:,:,jj*skip));  
%     axis equal; 
%     colormap(gray);
%     xlim(xlimits); 
%     ylim(ylimits);
%     caxis([cmin cmax]);
%    
%     if save_movie == 1
%         writeVideo(vid, getframe(gcf));
%     end
%     hold off;
% 
%     pause(1/fps_r);
% end
% 
% 
% if  save_movie ==1 
%     close(vid);
% end
% 
% 
% beep
% 
% 
% 
















%     %%
%     
%     % Shift image
% 
%     % Positive x_shift->left, positive y_shift->up
% 
%    
%    
%     for k = 1:N_frames
%     
%         xshift = -x_shift(k);
%         yshift = -y_shift(k);
%         
%         
%         % Mark measured displacement with black lines
%         FRAMES(:,217-xshift,k) = ones(width,1);
%         FRAMES(218-yshift,:,k) = ones(1,height);
%     
%     
%         if use_wbar == 1
%             waitbar(k/N_det,wbar,'Shifting images')
%         end
% 
% 
%         if (xshift > 0)
%             FRAMES(:,1:end-xshift,k) = FRAMES2(:,xshift+1:end,k);
%         elseif (xshift < 0)
%         	FRAMES(:,-xshift+1:end,k) = FRAMES2(:,1:end+xshift,k);
%         else
%             FRAMES(:,:,k) = FRAMES2(:,:,k);
%         end
% 
% 
%         if (yshift > 0)
%         	FRAMES(1:end-yshift,:,k) = FRAMES(yshift+1:end,:,k);
%         elseif (yshift < 0)
%         	FRAMES(-yshift+1:end,:,k) = FRAMES(1:end+yshift,:,k);
%         else
%             FRAMES(:,:,k) = FRAMES(:,:,k);
%         end
%     
%     end
%     
%     
%     
% %     for k = 1:N_frames-1
% %         if use_wbar == 1
% %             waitbar(k/N_det,wbar,'Shifting images')
% %         end
% %     % for k = 1535:1635
% %         x_shift_curr = -x_shift(k);
% %         y_shift_curr = -y_shift(k);
% % 
% %         % Shift x and y:
% %         centered_frame = FRAMES(max(y_shift_curr,1):min(y_shift_curr+height,height),max(x_shift_curr,1):min(x_shift_curr+width,width),k); 
% %         sz_cent = size(centered_frame);
% %         FRAMES2(max(-y_shift_curr,1):sz_cent(1)+max(-y_shift_curr,1)-1,max(-x_shift_curr,1):sz_cent(2)+max(-x_shift_curr,1)-1,k) = centered_frame;
% % 
% %     %     % Mark measured displacement with black lines
% %     %     FRAMES2(:,217+x_shift_curr,k) = zeros(width,1);
% %     %     FRAMES2(218+y_shift_curr,:,k) = zeros(1,height);
% %     end
% %     
% %     FRAMES = FRAMES2;
%     
%     clear xshift yshift blackman_filt_w blackman_filt_h filt_b filt_a NFFTX NFFTY freqx freqy ROW_F shift corr int_min int_max int_step img1 img2 FFTF F1 F2 R corr r max_val_i max_ind_i max_x_i max_y_i X1 X2 Xq1 Xq2 Vq max_val max_ind max_x max_y x_shift_curr y_shift_curr centered_frame FRAMES2 FRAMES2*
%     
% end
% 
% 
% 
% % figure; plot(x_shift); hold on; plot(y_shift,'r');
% % 
% % figure; subplot(2,2,1); imagesc(img1); colormap(gray); caxis([0 1000]);
% % subplot(2,2,2); imagesc(img2); colormap(gray); caxis([0 1000]);
% % subplot(2,2,3); imagesc(abs(((ifft2(F1))))); colormap(gray); 
% % subplot(2,2,4); imagesc(abs(((ifft2(F2))))); colormap(gray); 
% % 
% % 
% % figure; imagesc(abs(((ifft2(F1))))-abs(((ifft2(F2))))); colormap(gray);
% 
% % play_movie(FRAMES, CELLS, CELLS, CENT_X, CENT_Y, 20, 1, 0, 170,'false');
% 
% % play_movie(FRAMES2, CELLS, CELLS, CENT_X, CENT_Y, 20, 1, 0, 170,'false');
% 
% 
% % play_movie(FRAMES2, zeros(height,width), zeros(height,width), 0, 0, 1, 1, 0, 170,'false');
% 
% 
% % play_movie(FRAMES, zeros(height,width), zeros(height,width), 0, 0, 20, 1, 0, 170,'false');








%% Remove min frame

if ~register_image_flag
   
    % MIN_FRAME2 = MIN_FRAME;
    % MIN_FRAME_FILT = imfilter(MIN_FRAME2, GAUSS_C, 'same');


    % figure; imagesc(MIN_FRAME); axis equal; xlim([1,width]); ylim([1, height]); colormap(gray(256)); caxis([0,col_depth]);

    for k = 1:N_frames
        if use_wbar == 1 
            waitbar(k/N_frames,wbar,'Removing min frame')
        end
        FRAMES(:,:,k) = FRAMES(:,:,k) - MIN_FRAME;
        % FRAMES(:,:,k) = FRAMES(:,:,k) - AVG_FRAME;
        % FRAMES(:,:,k) = FRAMES(:,:,k) - MIN_FRAME_FILT;
    end
end



%% Perform gaussian filtering

for k = 1:N_frames
    if use_wbar
        waitbar(k/N_det,wbar,'Performing gaussian filtering')
    end
    % aaa = 0.1*kern_conv(GAUSS, FRAMES(:,:,k)); % slow
    FILT_IMG(:,:,k) = imfilter(FRAMES(:,:,k), GAUSS, 'same');
end



%% Calculate image gradient
for k = 1:N_det
    if use_wbar
        waitbar(k/N_det,wbar,'Calculating image gradient')
    end

    % Overloaded method ('Sobel', Prewitt', 'CentralDifference',
    % IntermediateDifference')

    [GRAD_X(:,:,k), GRAD_Y(:,:,k)] = imgradientxy(FILT_IMG(:,:,k),'Sobel');
    

    % Intermediate difference:
%     for m = 1:height-1
%         for n = 1:width-1
%             GRAD_X(m,n,k) = FILT_IMG(m,n+1,k) - FILT_IMG(m,n,k);
%             GRAD_Y(m,n,k) = FILT_IMG(m+1,n,k) - FILT_IMG(m,n,k);
%         end
%     end  

end



% clear FILT_IMG


%% Filter x and y gradient
for k = 1:N_det
    
%     if use_wbar == 1
%         waitbar(k/N_det,wbar,'Filtering x and y gradient components (step 5/9)')
%     end
%     
%     for m = 1:height
%         GRAD_X(m,:,k) = filtfilt(filt_b, filt_a,GRAD_X(m,:,k));
%         % GRAD_X(m,:,k) = filtfilt(butterfilt,GRAD_X(m,:,k));
%     end
%     for n = 1:width
%         GRAD_Y(:,n,k) = filtfilt(filt_b, filt_a,GRAD_Y(:,n,k));
%         % GRAD_Y(:,n,k) = filtfilt(butterfilt,GRAD_Y(:,n,k));
%     end
%    GRAD_X(:,:,k) = imfilter(GRAD_X(:,:,k), GAUSS, 'same');
%   GRAD_Y(:,:,k) = imfilter(GRAD_Y(:,:,k), GAUSS, 'same');


end


%% Apply threshold to gradient

% Calculate thresholds based on peaks
% Smaller thresholds detect shallower gradients (more out of focus cells)
% Denominator is typically 7-12

% thresh_px = max(max(max(GRAD_X)))/15;
% thresh_nx = -min(min(min(GRAD_X)))/15;
% thresh_py = max(max(max(GRAD_Y)))/15;
% thresh_ny = -min(min(min(GRAD_Y)))/15;

thresh_px = mean(mean(rms(GRAD_X,3)))*rms_fact;
thresh_nx = thresh_px;
thresh_py = mean(mean(rms(GRAD_Y,3)))*rms_fact;
thresh_ny = thresh_py;


% tmp = double(FRAMES);
% tmp2 = rms(tmp,3);
% tmp2 = imfilter(tmp2, GAUSS_C, 'same');

THRESH_X = mean(GRAD_X,3);
THRESH_X = imfilter(THRESH_X, GAUSS_C, 'same');
% THRESH_X = THRESH_X.*rms_fact;
THRESH_X = THRESH_X+rms_fact;

THRESH_Y = mean(GRAD_Y,3);
THRESH_Y = imfilter(THRESH_Y, GAUSS_C, 'same');
% THRESH_Y = THRESH_Y.*rms_fact;
THRESH_Y = THRESH_Y+rms_fact;







% % Use single threshold value (different varibles)
% 
% GRAD_XT = GRAD_X;
% GRAD_YT = GRAD_Y;
% 
% GRAD_XT(find(GRAD_XT<thresh_px & GRAD_XT > -thresh_nx)) = 0;
% GRAD_XT(GRAD_XT>= thresh_px) = 1;
% GRAD_XT(GRAD_XT<= -thresh_nx) = -1;
% 
% GRAD_YT(find(GRAD_YT<thresh_py & GRAD_YT > -thresh_ny)) = 0;
% GRAD_YT(GRAD_YT>= thresh_py) = 1;
% GRAD_YT(GRAD_YT<= -thresh_ny) = -1;


% Use single threshold value (replace GRAD_X GRAD_Y)

GRAD_X(find(GRAD_X<thresh_px & GRAD_X > -thresh_nx)) = 0;
GRAD_X(GRAD_X>= thresh_px) = 1;
GRAD_X(GRAD_X<= -thresh_nx) = -1;

GRAD_Y(find(GRAD_Y<thresh_py & GRAD_Y > -thresh_ny)) = 0;
GRAD_Y(GRAD_Y>= thresh_py) = 1;
GRAD_Y(GRAD_Y<= -thresh_ny) = -1;




% Use threshold matrix

% THRESH_X2 = repmat(THRESH_X,[1, 1, N_det]);
% THRESH_Y2 = repmat(THRESH_Y,[1, 1, N_det]);
% 
% 
% if use_wbar
%     waitbar(1/4,wbar,'Applying threshold to gradient')
% end
% TH_XP = GRAD_X(:,:,1:N_det) > THRESH_X2;
% if use_wbar
%     waitbar(2/4,wbar,'Applying threshold to gradient')
% end
% TH_XN = GRAD_X(:,:,1:N_det) < -THRESH_X2;
% if use_wbar
%     waitbar(3/4,wbar,'Applying threshold to gradient')
% end
% TH_YP = GRAD_Y(:,:,1:N_det) > THRESH_Y2;
% if use_wbar
%     waitbar(4/4,wbar,'Applying threshold to gradient')
% end
% TH_YN = GRAD_Y(:,:,1:N_det) < -THRESH_Y2;
% 
% 
% 
% GRAD_XT = double(TH_XP) - double(TH_XN);
% GRAD_YT = double(TH_YP) - double(TH_YN);








% for k = 1:N_det
%     if use_wbar
%         waitbar(k/N_det,wbar,'Applying threshold to gradient')
%     end
%     for m = 1:height
%         for n = 1:width
%              if (GRAD_X(m,n,k) >= THRESH_X(m,n)) 
%                 GRAD_XT(m,n,k) = 1;
%             end
%             if (GRAD_X(m,n,k)>-THRESH_X(m,n) && GRAD_X(m,n,k) < THRESH_X(m,n)) 
%                 GRAD_XT(m,n,k) = 0;
%             end
%             if (GRAD_X(m,n,k) <= -THRESH_X(m,n)) 
%                 GRAD_XT(m,n,k) = -1;    
%             end
% 
%             if (GRAD_Y(m,n,k) >= THRESH_Y(m,n)) 
%                 GRAD_YT(m,n,k) = 1;
%             end
%             if (GRAD_Y(m,n,k)>-THRESH_Y(m,n) && GRAD_Y(m,n,k) < THRESH_Y(m,n)) 
%                 GRAD_YT(m,n,k) = 0;
%             end
%             if (GRAD_Y(m,n,k) <= -THRESH_Y(m,n)) 
%                 GRAD_YT(m,n,k) = -1;    
%             end
%         end
%     end
% end







% clear threshold_* GRAD_X GRAD_Y THRESH_* TH_*




%% Match cell shape



for k = 1:N_det
     if use_wbar
         waitbar(k/N_det,wbar,'Matching gradient with cell shape')
     end



%     for m = 1:height
%         x_conv(m,:) = conv(GRAD_XT(m,:,k),cell_conv,'same');
%     end
%     x_conv(find(x_conv(:,:)< thresh_cell)) = 0;
%     CONV_X(:,:,k) = logical(x_conv);
% 
%     
%     for n = 1:width
%         y_conv(:,n) = conv(GRAD_YT(:,n,k),cell_conv,'same');
%     end
%     y_conv(find(y_conv(:,:)< thresh_cell)) = 0;
%     CONV_Y(:,:,k) = logical(y_conv);





    wait_neg = 0;

    for m = 1:height
        state = 0;
        cnt = 0;
        wait_neg = 0;
        for n = 1:width

            switch state
                case 0
                    if GRAD_X(m,n,k) == 1
                        state = 2;
                        cnt = 0;
                        wait_neg = 0;
                    end   
                case 2
                    if GRAD_X(m,n,k) == 1
                        state = 3;
                    else 
                        state = 0;
                    end
                case 3
                    cnt = cnt + 1;

                    if GRAD_X(m,n,k) == 0
                        wait_neg = 1;
                    end

                    if GRAD_X(m,n,k) == -1
                        wait_neg = 0;
                        state = 4; 
                    end

                    if (GRAD_X(m,n,k) == 1) && (wait_neg == 1)
                        state = 0;
                    end


                    if cnt > cell_size
                        state = 0;
                    end


                case 4
                    if GRAD_X(m,n,k) == -1
                        state = 6;
                    else 
                        state = 0;
                    end

                case 6 
                    % CONV_X(m,n-round(cnt/2),k) = 1;
                    cell_indx_min = max(1,n-round(cnt/2)-floor(cell_mask_size/2));
                    cell_indx_max = min(width,n-round(cnt/2)+floor(cell_mask_size/2));

                    CONV_X(m,cell_indx_min-1:cell_indx_max-1,k) = 1;
                    state = 0;
            end
        end
    end






     for m = 1:width
        state = 0;
        cnt = 0;
        wait_neg = 0;
        for n = 1:height

            switch state
                case 0
                    if GRAD_Y(n,m,k) == 1
                        state = 2;
                        cnt = 0;
                        wait_neg = 0;
                    end

                case 2
                    if GRAD_Y(n,m,k) == 1
                        state = 3;
                    else
                        state = 0;
                    end
                 case 3
                    cnt = cnt +1;


                    if GRAD_Y(n,m,k) == 0
                        wait_neg = 1;
                    end

                    if GRAD_Y(n,m,k) == -1
                        wait_neg = 0;
                        state = 4; 
                    end

                    if (GRAD_Y(n,m,k) == 1) && (wait_neg == 1)
                        state = 0;
                    end


                    if cnt > cell_size
                        state = 0;
                    end

                case 4
                    if GRAD_Y(n,m,k) == -1
                        state = 6;
                    else 
                        state = 0;
                    end

                case 6 
                    % CONV_Y(n-round(cnt/2),m,k) = 1;


                    cell_indx_min = max(1,n-round(cnt/2)-floor(cell_mask_size/2));
                    cell_indx_max = min(height,n-round(cnt/2)+floor(cell_mask_size/2));

                    CONV_Y(cell_indx_min-1:cell_indx_max-1,m,k) = 1;
                    state = 0;

            end

        end

    end





end


clear GRAD_X GRAD_Y


%% Match x-y locations and time persistence


CELLS2 = logical(zeros(height,width,N_det));
CELLS_DYN_SINGLE = logical(zeros(height,width,N_det));
curr_cell = logical(zeros(height, width));

for k = 1:N_det
%     if use_wbar == 1
%         waitbar(k/N_det,wbar,'Matching x-y cell locations')
%     end

    % [indx_y indx_x] = find(CONV_Y(:,:,k)>0);
    CELLS2(:,:,k) = CONV_X(:,:,k) & CONV_Y(:,:,k);
    % CELLS2(:,:,k) =  CONV_Y(:,:,k);


    % Accumulate cells masks
    if (k > T_cell)
        % CELLS(:,:,k) = CELLS(:,:,k-1) | CELLS(:,:,k);



        tmp_cell = logical(ones(height, width));

        for m = 1:T_cell
            curr_cell = CELLS2(:,:,k-m);
            tmp_cell = tmp_cell & curr_cell;
        end
        CELLS(:,:,k) = CELLS(:,:,k-1) | tmp_cell;
        CELLS_DYN_SINGLE(:,:,k) = tmp_cell;
    end
end

clear tmp_cell curr_cell CELLS2 CONV_X CONV_Y

% figure; imagesc(FRAMES(:,:,N_det)); colormap(gray); caxis([0,100]); hold on;
% plot(indx_x, indx_y, 'ro');



%% Identify cells and create cell mask


[identified_cells, CELL_AREAS, CENT_X, CENT_Y, CELL_REGIONS, CELLS_FILT] = find_cells(CELLS(:,:,N_det), CELL_MASK, area_min, min_cell_distance);

CELLS_DYN = logical(CELLS_DYN_SINGLE);
CELLS_DYN_SUM = logical(CELLS);
CELLS = logical(CELLS(:,:,N_det));







%% Visulaize results 

% subplot(2,2,3); imagesc(GRAD_X); colormap(gray); caxis([0,50]);
% subplot(2,2,4); imagesc(GRAD_Y); colormap(gray); caxis([0,50]);
% 
% figure; plot(FRAMES(195:205,:,10)'); hold on; % plot(FILT_IMG(200,:),'r');
% 
% 
% 
% 
% figure; plot(GRAD_Y(:,135));
% 
% figure; imagesc(FRAMES(:,:,curr_frame)); colormap(gray); caxis([0,500]);hold on;


% vis_frame = 328;
% vis_col = 223;
% vis_row = 171;
% 
% figure; imagesc(FRAMES(:,:,vis_frame)); colormap(gray); 
% hold on; line([1:400], vis_row.*ones(1,400),'Color','red');
% line(vis_col.*ones(1,400),[1:400],'Color','red');
% 
% 
% figure; plot(FRAMES(vis_row,:,vis_frame)); hold on; plot(FRAMES(:,vis_col,vis_frame),'r');
% xlabel('Pixel number'); ylabel(['Pixel value'])
% legend('Selected row', 'Selected column');
% 
% figure; plot(GRAD_X(vis_row,:,vis_frame)); hold on; plot(400*GRAD_XT(vis_row,:,vis_frame),'r'); plot(400*CONV_X(vis_row,:,vis_frame),'g');
% xlabel('Pixel number'); ylabel(['Gradient'])
% legend('X gradient', 'Filtered x gradient','Identified cell');
% 
% figure; plot(GRAD_Y(:,vis_col,vis_frame)); hold on; plot(400*GRAD_YT(:,vis_col,vis_frame),'r'); plot(400*CONV_Y(:,vis_col,vis_frame),'g');
% xlabel('Pixel number'); ylabel(['Gradient'])
% legend('Y gradient', 'Filtered y gradient','Identified cell');
% 
% 
% figure; imagesc(FRAMES(:,:,vis_frame)); colormap(gray); hold on;
% red = cat(3, ones(400, 400), zeros(400, 400), zeros(400, 400));
% green = cat(3, zeros(400, 400), ones(400, 400), zeros(400, 400));
% blue = cat(3, zeros(400, 400), zeros(400, 400), ones(400, 400));
% hh = imshow(green);
% set(hh, 'AlphaData', 0.85.*CELLS(:,:,end));
% hold on;plot(CENT_X,CENT_Y,'or')




% figure; subplot(1,2,1); imagesc(FRAMES(:,:,vis_frame)); colormap(gray); 
% subplot(1,2,2); imagesc(FILT_IMG(:,:,vis_frame)); colormap(gray);
% figure; plot(FILT_IMG(vis_row,:,vis_frame)); hold on; plot(GRAD_X(vis_row,:,vis_frame),'r');
% figure; plot(GRAD_X(vis_row,:,vis_frame)); hold on ; plot(100*GRAD_XT(vis_row,:,vis_frame),'r'); plot(100*CONV_X(vis_row,:,vis_frame),'g');
% 
% figure; plot(FILT_IMG(:,vis_col,vis_frame)); hold on; plot(GRAD_Y(:,vis_col,vis_frame),'r');
% figure; plot(GRAD_Y(:,vis_col,vis_frame)); hold on ; plot(100*GRAD_YT(:,vis_col,vis_frame),'r'); plot(100*CONV_Y(:,vis_col,vis_frame),'g');
% 





% green = cat(3, zeros(height, width), ones(height, width), zeros(height, width));
% blue = cat(3, zeros(height, width), zeros(height, width), ones(height, width));
% 
% 
% figure; imagesc(FILT_IMG(:,:,vis_frame)); colormap(gray); hold on; 
% h = imshow(green); 
% hh = imshow(blue); 
% set(h, 'AlphaData', 0.2*CONV_X(:,:,end));
% set(hh, 'AlphaData', 0.2*CONV_Y(:,:,end));
% 
% 
% figure; imagesc(CELLS_FILT); colormap(gray);



 % green = cat(3, zeros(height, width), ones(height, width), zeros(height, width));
 % hold on 
 % h = imshow(green); 
 % hold off
 % set(h, 'AlphaData', CONV_X)


% figure; imagesc(uint16(FRAMES(:,:,end))); colormap(gray(256)); 
% axis equal; xlim([1,width]); ylim([1, height]);
% caxis([0,col_depth]); hold on;

% green = cat(3, zeros(height, width), ones(height, width), zeros(height, width));
% hh = imshow(green);
% set(hh, 'AlphaData', CONV_X(:,:,end)); 

% plot(CENT_X,CENT_Y,'go');



% a = zeros(height, width,3);
% a(:,:,1) = CELLS(:,:,3000)-CELLS(:,:,2000);
% a(:,:,2) = CELLS(:,:,2000);
% a(:,:,3) = zeros(height, width);
% figure; imagesc(a)
% 
% 
% a = zeros(height, width,3);
% a(:,:,1) = CELLS(:,:,end);
% a(:,:,2) = CELLS_FILT(:,:);
% a(:,:,3) = zeros(height, width);
% figure; imagesc(a)
% 
% figure; imagesc(CELLS(:,:,end)); colormap(gray); hold on; plot(CENT_X, CENT_Y,'go');





%% Save current workspace for each session


if use_wbar
    close(wbar);    
end




CNT_BAD_FRAMES = cnt_bad_frames;
X_SHIFT = x_shift;
Y_SHIFT = y_shift;
IDENTIFIED_CELLS = identified_cells; 

% FRAMES = uint16(FRAMES);



%         clear  FRAMES2* CURR_CELL_AREAS CURR_CELLS_FILT CURR_CENT_X CURR_CENT_Y curr_dist curr_identified_cells CURR_REGIONS  FILT_IMG folder GRAD_X GRAD_XT GRAD_Y GRAD_YT k m  min_dist  n N_frames str_ind  x_conv y_conv CELL_MASK_*
% 
% 
%         workspace_var = whos('*');
% 
% 
%         for k = 1:numel(workspace_var)
%             var_name = workspace_var(k).name;
%             new_var_name = strcat(var_name, header);
%             % if (strcmp('exp_numb', var_name) == 0) & (strcmp('N_days', var_name) == 0) & (strcmp('header', var_name) == 0) & (strcmp('batch_proc', var_name) == 0)
%             %     eval([new_var_name, ' = ', var_name, '; clear  ', var_name])
%             % end
%             eval([new_var_name, ' = ', var_name, '; '])
%         end
% 
% 
% 
%         display('Current data saved to workspace.');











%         clear cell_indx_* CELLS CELLS_FILT CENT_X CENT_Y ff* file* GAUSS* k MAX_FRAME MIN_FRAME AVG_FRAME x_shift_1 y_shift_1 x_shift_f y_shift_f FRAMES col_depth_* day_numb_* dir_content_* height_* mouse_numb_* N_det_* N_ss_* play_mov* thresh_nx_* thresh_ny_* thresh_px_* thresh_py_* use_wbar* wbar* width_* exp_numb_* filt_a_* filt_b_* cell_type_* MacOS* cell_conv_* thresh_cell_* area_min* aw_co* min_cell_distance_* slsh*
%         clear x_shift* y_shift* wait_neg* thresh* T_cell_* sz_* state* remove_* register_* N_days_* header* cell_mask_size* cell_size* baseDir*
%         clear FRAMES2*
% 
% 
%         wspace = strcat('data_exp_', num2str(exp_numb), '.mat');
%         display('All done. Workspace saved to disk.');
%         save(wspace, '-v7.3');
%        
%         beep
%         pause(1)
%         beep
%         pause(0.5)
%         beep



 curr_time = round(toc);
 curr_min = floor(curr_time/60);
 curr_s = rem(curr_time,60);
 
 fprintf('Elapsed time is %d min and %d s.\n', curr_min, curr_s)




%% Play movie

if show_movie
    play_movie(FRAMES, CELLS, CELLS, CENT_X, CENT_Y, 20, 1, 0, 170,'false');
end





     