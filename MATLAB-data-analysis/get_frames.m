function [FRAMES] = get_frames(MacOS, use_wbar, remove_bad_frames, register_image, filePath, DARK_FRAME, N_det, width, height)
%GET_FRAMES Load raw images and apply image registration algorithm
%
% MacOS:                    (bool) 1 for OSX, 0 for Windows        
% use_wbar                  (bool) use status bar
% remove_bad_frames         (bool) remove and mark bad frames
% register_image_flag       (bool) perform image registration
% filePath                  (string array) full path
% N_det                     (scalar) numer of frames used for cell detection
% width                     (scalar) image width
% height                    (scalar) image height
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


FILT_IMG = zeros(height, width, N_det,'uint16');


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

% % Use most frequent value for very large datasets (mostly dark, N_frames>3000)
% MIN_FRAME = uint16(zeros(height,width));
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

MIN_FRAME = uint16(zeros(height,width));
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
if ~isempty(frame_idx([frame_idx<size(FRAMES,3)] & [frame_idx>min(10,size(FRAMES,3))]))
    fprintf('WARNING: large intensity changes detected, possible LED issue.');
    figure; plot(frame_mean); ylabel('Average pixel intensity'); xlabel('Frame number');
    hold on;
    frame_err = nan(1,size(FRAMES,3));
    frame_err(frame_idx([frame_idx<size(FRAMES,3)] & [frame_idx>min(10,size(FRAMES,3))])) = frame_mean(frame_idx([frame_idx<size(FRAMES,3)] & [frame_idx>min(10,size(FRAMES,3))]));
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
        MIN_FRAME(n,m) = nanmin(FRAMES_2(n,m,min(10,size(FRAMES,3)):max(1,size(FRAMES,3)-5)));
    end
end
MIN_FRAME = uint16(MIN_FRAME);


    
    


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

    b_width = 1.6*width;
    b_height = 1.6*height;
    blackman_filt_w = blackman(b_width,'symmetric');
    blackman_filt_h = blackman(b_height,'symmetric');
    ff=(blackman_filt_w*blackman_filt_h'); 

    ff = ff(floor(b_width-width)/2:floor(b_width-width)/2 + width-1, floor(b_height-height)/2:floor(b_height-height)/2 + height-1);

    [filt_b,filt_a] = butter(16,0.2,'low');


    NFFTX = width;%1024; %2^nextpow2(width); % Next power of 2 from width
    freqx = 1/2*linspace(0,1,NFFTX/2+1);
    NFFTY = height; %1024; %2^nextpow2(width); % Next power of 2 from width
    freqy = 1/2*linspace(0,1,NFFTY/2+1);

    corr = zeros(NFFTY, NFFTY, N_frames);

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

    corr = zeros(NFFTY, NFFTY,N_frames);

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
    F1 = fft2(img1,NFFTX,NFFTY);
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
        F2 = fft2(img2,NFFTX,NFFTY);

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

if ~register_image
   
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




if use_wbar == 1
    close(wbar)
end
fprintf('\n');



