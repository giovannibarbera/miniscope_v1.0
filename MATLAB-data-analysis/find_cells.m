function [identified_cells, CELL_AREAS, CENT_X, CENT_Y, REGIONS, CELLS_FILT] = find_cells(CELLS, CELL_MASK, area_min, dist_min)
%FIND_CELLS Locate cell position in a binary mask. It also filters and
%sorts cells by area.
%
% CELLS:        (height,width), cell binary mask
% CELL_MASK:    (m, m), binary mask to associate with each cell centroid
% area_min:     scalar, minimum cell area [pixel]
% area_min:     scalar, minimum cell distance [pixel]
% 
% identified cells: scalar, total number of identified (filtered) cells N
% CELL_AREAS:       (N), cell area
% CENT_X:           (N), x position of cell centroid
% CENT_Y:           (N), y position of cell centroid
% REGIONS:          (cell), pixel list of every cell
% CELLS_FILT:       (height, width), filtered cell mask
%
%   giobarbera@neuralmapper.com

size_cells = size(CELLS);
width = size(CELLS,2);
height = size(CELLS,1);
size_mask = size(CELL_MASK,1); % should be square
CELLS_FILT = zeros(size_cells(1), size_cells(2));
IMG_FILL = imfill(CELLS,'holes');
CellsLabel = bwlabel(IMG_FILL);
CENT = regionprops(CellsLabel,'centroid');
CELL_AREAS = regionprops(CellsLabel,'Area');
CELL_REGIONS = regionprops(CellsLabel,'PixelList');
AREAS = 0;

% Find cells with area < min_area
rm_indx = 0;
for p = 1:numel(CENT)
    AREAS(p) = CELL_AREAS(p).Area;
    if AREAS(p) < area_min
        rm_indx = [rm_indx; p];
    end
end

% [CELL_AREAS sort_ind] = sort(AREAS);
% indx_min = find(CELL_AREAS<area_min);
% if isempty(indx_min)
%     indx_min = [0];
% end
% sort_ind = sort_ind(indx_min(end)+1:end);






% Remove entries with area < area_min

identified_cells = numel(CENT) - length(rm_indx) +1;
CENT_X = zeros(1, identified_cells);
CENT_Y = zeros(1, identified_cells);
REGIONS = cell(identified_cells);


if length(rm_indx) > 1
    curr_rm_ind = 2;
    curr_ind = 1;
    
    for k = 1:numel(CENT)

        if k ~= rm_indx(curr_rm_ind)
            CENT_X(curr_ind) = CENT(k).Centroid(1);
            CENT_Y(curr_ind) = CENT(k).Centroid(2);
            REGIONS{curr_ind} = CELL_REGIONS(k).PixelList;
            curr_region = REGIONS{curr_ind};
            sz_curr_region = size(curr_region);

            for m = 1:sz_curr_region(1)
                CELLS_FILT(curr_region(m,2),curr_region(m,1)) = 1;
            end
            clear curr_region;

            curr_ind = curr_ind + 1;

        else 
            curr_rm_ind = min(curr_rm_ind + 1, length(rm_indx));
        end
    end

else
     for k = 1:numel(CENT)

            CENT_X(k) = CENT(k).Centroid(1);
            CENT_Y(k) = CENT(k).Centroid(2);
            REGIONS{k} = CELL_REGIONS(k).PixelList;
            curr_region = REGIONS{k};
            sz_curr_region = size(curr_region);

            for m = 1:sz_curr_region(1)
                CELLS_FILT(curr_region(m,2),curr_region(m,1)) = 1;
            end
            clear curr_region;


     end
end





% Remove cell too close to others
rm_ind = 0;
not_done = 1;
curr_dist = 0;
if length(CENT_X) >4
    while (length(rm_ind) > 1 | not_done == 1)
        for n = 1:length(CENT_X)-1
            rm_ind = n;
            for m = (n+1):length(CENT_X)
                curr_dist = sqrt(((CENT_X(n)-CENT_X(m))^2)  + ((CENT_Y(n)-CENT_Y(m))^2));
                if curr_dist < dist_min
                    rm_ind = [rm_ind; m];
                end
            end

            if length(rm_ind) > 1

                % Find centroid of cell cluster and replace first cell too close to
                % others
                CENT_X(rm_ind(1)) = round(mean(CENT_X(rm_ind)));
                CENT_Y(rm_ind(1)) = round(mean(CENT_Y(rm_ind)));

                % Remove cells too close to each other
                for k = 2:length(rm_ind)
                    CENT_X = [CENT_X(1:rm_ind(k)-k+1),CENT_X(rm_ind(k)-k+3:end)];
                    CENT_Y = [CENT_Y(1:rm_ind(k)-k+1),CENT_Y(rm_ind(k)-k+3:end)];
                end
                break;
            end
        end
        not_done = 0;
    end
end








identified_cells = length(CENT_X);



% % Adjust offset
% CENT_X = CENT_X + 1;
% CENT_Y = CENT_Y + 1;
% CENT_X(CENT_X>size(CELLS,2)) = size(CELLS,2);
% CENT_Y(CENT_Y>size(CELLS,1)) = size(CELLS,1);




% 
% CELL_SHAPE = logical([0 0 1 1 0 0;
%               0 1 1 1 1 0;
%               1 1 1 1 1 1;
%               1 1 1 1 1 1;
%               0 1 1 1 1 0;
%               0 0 1 1 0 0]);
% 
% CELLS_FILT = zeros(size_cells(1), size_cells(2));
% for k = 1:identified_cells
%     CELLS_FILT(floor(CENT_Y(k))-2:floor(CENT_Y(k))+3, floor(CENT_X(k))-2:floor(CENT_X(k))+3) = CELLS_FILT(floor(CENT_Y(k))-2:floor(CENT_Y(k))+3, floor(CENT_X(k))-2:floor(CENT_X(k))+3) | CELL_SHAPE;
% end


% CELL_SHAPE = logical([1 1;
%                       1 1]);
% CELLS_FILT = zeros(size_cells(1), size_cells(2));
% for k = 1:identified_cells
%     CELLS_FILT(floor(CENT_Y(k))-1:floor(CENT_Y(k)), floor(CENT_X(k))-1:floor(CENT_X(k))) = CELLS_FILT(floor(CENT_Y(k))-1:floor(CENT_Y(k)), floor(CENT_X(k))-1:floor(CENT_X(k))) | CELL_SHAPE;
% end


CELL_SHAPE = logical(CELL_MASK);
center = floor(size_mask/2);
offset = mod(size_mask,2);
if offset == 0
    offset = 1;
else 
    offset = 0;
end


CELLS_FILT = zeros(size_cells(1), size_cells(2));
for k = 1:identified_cells
    % Exclude cells close to the border
    if CENT_X(k) > center +1 && CENT_X(k) < width - center -1 && CENT_Y(k) > center +1 && CENT_Y(k) < height - center -1 
        CELLS_FILT(floor(CENT_Y(k))-center:floor(CENT_Y(k))+center-offset, floor(CENT_X(k))-center:floor(CENT_X(k))+center-offset) = CELLS_FILT(floor(CENT_Y(k))-center:floor(CENT_Y(k)+center-offset), floor(CENT_X(k))-center:floor(CENT_X(k))+center-offset) | CELL_SHAPE;

    end
end





% identified_cells =  length(sort_ind);
% 
% CENT_X = zeros(1, identified_cells);
% CENT_Y = zeros(1, identified_cells);
% % PERIMETER = zeros(1, identified_cells);
% REGIONS = cell(identified_cells);
% for k = 1:identified_cells
% %     if use_wbar == 1
% %         waitbar(k/identified_cells,wbar,'Filtering cells by area (step 9/10)')
% %     end
%     CENT_X(k) = CENT(sort_ind(k)).Centroid(1);
%     CENT_Y(k) = CENT(sort_ind(k)).Centroid(2);
%     REGIONS{k} = CELL_REGIONS(sort_ind(k)).PixelList;
%     curr_region = REGIONS{k};
%     for m = 1:length(curr_region)
%         CELLS_FILT(curr_region(m,2),curr_region(m,1)) = 1;
%     end
%     clear curr_region;
% end
% 
% 
