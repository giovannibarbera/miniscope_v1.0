function [N_days, N_ss, T, width, height, days, sessions, cell_mask_size, cell_size, min_cell_distance, rms_fact, area_min, T_cell, N_det, MIN_FRAME, MIN_FRAME_FILT, AVG_FRAME, X_SHIFT, Y_SHIFT, CENT_X_GLOBAL, CENT_Y_GLOBAL, IDX_MAP, CA_TRACES_RAW_G, CA_TRACES_FILT_G, CA_TRACES_BIN_G, CA_TRACES_ONSET_G, clust_idx, clust_idx_s, sort_idx,  CLUST_CONN_FULL, clust_conn_thresh, clust_iterations, clust_bin_size, clust_curr_days, clust_curr_sessions, clust_metric, clust_tot_kmeans_iter, clust_max_iter, global_clusters, eartag, ref_day, ref_angle] = load_mouse_light(MOUSE)
%LOAD_MOUSE_FULL Load full data from mouse
%
%   giobarbera@neuralmapper.com


%% Load single mouse data


METADATA = MOUSE{1};
N_days = length(MOUSE{1}{1});
width = METADATA{5}(1);
height = METADATA{5}(2);
days = METADATA{1};
sessions = METADATA{2};
T = METADATA{8};
eartag = METADATA{10};

N_ss = 0;           % Maximum number of session in a day
max_day = 0;            % Day with maximum number of sessions
for k = 1:length(days)
   N_ss = max(N_ss,length(sessions{k}));
   if length(sessions{k})==N_ss
       max_day = k;
   end
end


cell_mask_size = METADATA{3}(1);
cell_size = METADATA{3}(2);
min_cell_distance = METADATA{3}(3);
rms_fact = METADATA{3}(4);
area_min = METADATA{3}(5);
T_cell = METADATA{3}(6);
N_det = METADATA{4}(1);
% ref_reg = MOUSE{16}{4};
% if ~isempty(ref_reg)
%     ref_day = MOUSE{16}{4}(1);
%     ref_angle = MOUSE{16}{4}(2);
% else
    ref_day = [];
    ref_angle = [];
% end

% FRAMES = MOUSE{2};
MIN_FRAME = MOUSE{3};
MIN_FRAME_FILT = MOUSE{4};
AVG_FRAME = MOUSE{5};
X_SHIFT = MOUSE{6};
Y_SHIFT = MOUSE{7};
CELLS_DYN = MOUSE{8};
CELLS_DYN_SUM = MOUSE{9};
CELLS = MOUSE{10};
CELLS_FILT = MOUSE{11};
CENT_X = MOUSE{12};
CENT_Y = MOUSE{13};

% CENT_X_DAY = MOUSE{17}{3};
% CENT_Y_DAY = MOUSE{17}{4};
% CA_TRACES_RAW = MOUSE{18}{1};
% CA_TRACES_FILT = MOUSE{18}{2};
% CA_TRACES_BIN = MOUSE{18}{3};
% CA_TRACES_ONSET = MOUSE{18}{4};

if ~isempty(MOUSE{20}) && ~isempty(MOUSE{21})
    CENT_X_GLOBAL = MOUSE{20}{1};
    CENT_Y_GLOBAL = MOUSE{20}{2};
    IDX_MAP = MOUSE{20}{3};
    CA_TRACES_RAW_G = MOUSE{21}{1};
    CA_TRACES_FILT_G = MOUSE{21}{2};
    CA_TRACES_BIN_G = MOUSE{21}{3};
    CA_TRACES_ONSET_G = MOUSE{21}{4};
else
    CENT_X_GLOBAL = [];
    CENT_Y_GLOBAL = [];
    IDX_MAP = [];
    CA_TRACES_RAW_G = [];
    CA_TRACES_FILT_G = [];
    CA_TRACES_BIN_G = [];
    CA_TRACES_ONSET_G = [];
end


clust_idx = MOUSE{19}{1};
clust_idx_s = MOUSE{19}{2};
sort_idx = MOUSE{19}{3};
CLUST_CONN_FULL =  MOUSE{19}{4};
clust_conn_thresh = MOUSE{19}{5};
clust_iterations = MOUSE{19}{6};
clust_bin_size = MOUSE{19}{7};
clust_curr_days = MOUSE{19}{8};
clust_curr_sessions = MOUSE{19}{9};
clust_metric = MOUSE{19}{10};
clust_tot_kmeans_iter = MOUSE{19}{11};
clust_max_iter = MOUSE{19}{12};
global_clusters = MOUSE{19}{13};




