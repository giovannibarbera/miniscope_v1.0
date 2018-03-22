function [SPIKES F ] = cell_activity3(identified_cells, CENT_X, CENT_Y, ROI, CIRCLE, FRAMES, MIN_FRAME)
%CELL_ACTIVITY3 Evaluate activity of each cell in the mask CELLS, by
% averaging the pixel values over a fixed square area around the cell. 
%
% identified cells: scalar, total number of identified (filtered) cells N
% CENT_X:           (identified_cells), x position of cell centroid
% CENT_Y:           (identified_cells), y position of cell centroid
% ROI:              (m,n), mask of the region of interest around the cell
% FRAMES:           (height,width,N_frames), image data (background removed)
% 
% SPIKES:           (N, identified_cells), average normalized (deltaF/F) pixel intensity per cell
% F:                (identified_cells), reference fluorescence value for each cell
% 
%   giobarbera@neuralmapper.com

size_frames = size(FRAMES);
height = size(FRAMES,1);
width = size(FRAMES,2);
SPIKES = zeros(size_frames(3),identified_cells);
N = size(FRAMES,3);
F = zeros(identified_cells,1);
a = size(ROI,1);
b = size(ROI,2);
aa = size(CIRCLE,1);
bb = size(CIRCLE,2);
c = sum(sum(ROI));
F = zeros(N,identified_cells);
F0 = zeros(1,identified_cells);
zeros_roi = reshape(ROI,[1,a*b]);
zeros_circle = reshape(CIRCLE,[1,aa*bb]);
sz_circle = length(find(zeros_circle==1)); 
mean_bkgnd = mean(mean(MIN_FRAME(round(height/2)-25:round(height/2)+25,round(height/2)-25:round(height/2)+25)));


% Avoids rounding errors when CENT_X(k) = XX.5
CENT_X = CENT_X - 0.000001;
CENT_Y = CENT_Y - 0.000001;



% cell_act = figure;
% cell_act2 = figure;
% cell_act3 = figure;
% for p = [127 128] %

for p = 1:identified_cells 

    
   
   
    
    curr_mask = false(size(FRAMES,1),size(FRAMES,2));
    curr_circ = false(size(FRAMES,1),size(FRAMES,2));
    FRAMES2 = zeros(a,b,N); % ROI
    FRAMES3 = zeros(aa,bb,N); % Ring
    if (CENT_X(p) > (bb/2+2)) & (CENT_X(p) < size(FRAMES,2)-(bb/2+2)) & (CENT_Y(p) > (aa/2+2)) & (CENT_Y(p) < size(FRAMES,1)-(aa/2+2))
        curr_mask(max(round(CENT_Y(p)-floor(a/2)),1):min(round(CENT_Y(p)+ceil(a/2)-1),size(FRAMES,1)), max(round(CENT_X(p)-floor(b/2)),1):min(round(CENT_X(p)+ceil(b/2)-1),size(FRAMES,2))) = ROI;
        curr_circ(max(round(CENT_Y(p)-floor(aa/2)),1):min(round(CENT_Y(p)+ceil(aa/2)-1),size(FRAMES,1)), max(round(CENT_X(p)-floor(bb/2)),1):min(round(CENT_X(p)+ceil(bb/2)-1),size(FRAMES,2))) = CIRCLE;
        FRAMES2(:,:,:) = FRAMES(max(round(CENT_Y(p)-floor(a/2)),1):min(round(CENT_Y(p)+ceil(a/2)-1),size(FRAMES,1)), max(round(CENT_X(p)-floor(b/2)),1):min(round(CENT_X(p)+ceil(b/2)-1),size(FRAMES,2)),:);
        FRAMES3(:,:,:) = FRAMES(max(round(CENT_Y(p)-floor(aa/2)),1):min(round(CENT_Y(p)+ceil(aa/2)-1),size(FRAMES,1)), max(round(CENT_X(p)-floor(bb/2)),1):min(round(CENT_X(p)+ceil(bb/2)-1),size(FRAMES,2)),:);

    else
        curr_mask = false(size(FRAMES,1),size(FRAMES,2));
        curr_circ = false(size(FRAMES,1),size(FRAMES,2));
        FRAMES2 = zeros(a,b,N);
        FRAMES3 = zeros(aa,bb,N);
    end
    
    
%        figure; imagesc(FRAMES2(:,:,10));
%        figure; imagesc(FRAMES3(:,:,10));

    FRAMES2 = reshape(FRAMES2,[a*b,N]);
    FRAMES2(zeros_roi==0,:) = nan;
    
    FRAMES3 = reshape(FRAMES3,[aa*bb,N]);
    FRAMES3(zeros_circle==0,:) = nan; 
    
%         figure; imagesc(FRAMES2(:,10));
%         figure; imagesc(FRAMES3(:,10));   
%         FRAMES2b = reshape(FRAMES2,[a,b,N]);
%         FRAMES3b = reshape(FRAMES3,[aa,bb,N]);
%         figure; imagesc(FRAMES2b(:,:,10));
%         figure; imagesc(FRAMES3b(:,:,10));
    
    % SPIKES(:,p) = nanmax(FRAMES2,[],1);
    SPIKES(:,p) = nanmean(FRAMES2,1);
    gf = nanmean(FRAMES2,1);
    
    % F(:,p) = sum(FRAMES3,1); F = F./sz_circle;
    
    
    
    
    F(:,p) = nanmean(FRAMES3,1); 
    % F(:,p) = nanmin(FRAMES3,[],1); 


   
    
    % Saturation
    % F(F>200)=0;

    aaa = 40; bbb = 40;
    F00 = MIN_FRAME(max(round(CENT_Y(p)-floor(aaa/2)),1):min(round(CENT_Y(p)+ceil(aaa/2)-1),size(FRAMES,1)), max(round(CENT_X(p)-floor(bbb/2)),1):min(round(CENT_X(p)+ceil(bbb/2)-1),size(FRAMES,2)));
    F00 = double(F00);
    F00(F00<1) = nan;
    F0(p) = nanmean(nanmean(F00));
    
   %  fprintf('%d\n',p);
    
    % SPIKES(:,p) = mean(mean(FRAMES2,2),1);
    % SPIKES(:,p) = max(max(FRAMES2));
    
    SPIKES(:,p) = 100.*((SPIKES(:,p))-F(:,p))./F0(p);
    % SPIKES(:,p) = 100.*((SPIKES(:,p))-F(:,p))./mean_bkgnd;



%         % curr_frames = 75;%
%         curr_frames = 850;%
%         % curr_frames = [735 771];% 797];
%         for curr_frame = curr_frames 
%             figure; imagesc(FRAMES(:,:,curr_frame)); colormap(gray); hold on; plot(CENT_X(p),CENT_Y(p),'or');
%             ff = double(FRAMES(:,:,curr_frame));
%             ff = ff+100.*double(curr_mask) + 100.*double(curr_circ); 
%             
% 
%             figure; imagesc(ff); colormap(gray); caxis([0 400]);
%         
%             
%             figure(cell_act); hold on;
%             plot(0.1.*[1:size(FRAMES,3)],100.*nanmean(FRAMES2,1)./F0(p));
%             line([0.1.*curr_frame   0.1.*curr_frame],get(gca,'YLim'),'Color','k','LineStyle','--')
%             xlim([0.1*min(curr_frames)-1, 0.1*max(curr_frames)+7]);
%             title('Average soma');
% 
% 
%             figure(cell_act2); hold on;
%             plot(0.1.*[1:size(FRAMES,3)],SPIKES(:,p));
%             line([0.1.*curr_frame   0.1.*curr_frame],get(gca,'YLim'),'Color','k','LineStyle','--')
%             xlim([0.1*min(curr_frames)-1, 0.1*max(curr_frames)+7]);
%             title('Extracted trace');
% 
%             figure(cell_act3); hold on;
%             plot(0.1.*[1:size(FRAMES,3)],nanmean(FRAMES3,1));
%             line([0.1.*curr_frame   0.1.*curr_frame],get(gca,'YLim'),'Color','k','LineStyle','--')
%             xlim([0.1*min(curr_frames)-1, 0.1*max(curr_frames)+7]);
%              title('Annular region');
% 
% 
%             %    % Plot ROI
%             FRAMES2 = reshape(FRAMES2,[a b N]);
%             FRAMES3 = reshape(FRAMES3,[aa bb N]);
%             figure; imagesc(FRAMES2(:,:,curr_frame)); colormap(gray); caxis([0 400]); axis equal;
%             figure; imagesc(FRAMES3(:,:,curr_frame)); colormap(gray); caxis([0 400]); axis equal;
%         end
%         figure(cell_act); hold on;




end


