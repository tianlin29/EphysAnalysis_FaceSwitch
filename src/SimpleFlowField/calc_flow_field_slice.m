function data = calc_flow_field_slice(Tstamp, raster, ID, opt)
% function data = calc_flow_field_slice(Tstamp, raster, ID, opt)
% Tstamp = 1 x time
% raster = cell x time x trial
% ID = cell x trial
%
% Sliding window options:
% - window_size: proportion of trials in each window (e.g., 0.1 = 10% of trials)
% - window_step: proportion of trials to slide (e.g., 0.01 = 1% of trials step)

def.PC_coef = [];
def.detPSTH = [];

def.t_range = [0 600];
def.flow_field_dim = [1 2];
def.slice_dim = 3;
def.window_size = 0.1;  % proportion of trials in each window
def.window_step = 0.01; % proportion of trials to slide

def.conv_kernel = [];
def.bin_num = 20;

opt = safeStructAssign(def, opt);

if isempty(opt.PC_coef)
    error('opt.PC_coef should be filled');
end

ndim = 3;

%% smooth raster
raster = raster * 1e3; % cell x time x trial

% detrend
if ~isempty(opt.detPSTH)
    raster = raster - opt.detPSTH;
end

psth = nan(size(raster,1), size(raster, 2), nanmax(ID(:)));
for id = 1:nanmax(ID(:))
    psth(:, :, id) = nanmean(raster(:, :, ID(1,:) == id),3);
end

%% convert to PCA
T = opt.t_range(1):opt.t_range(2);
ind = Tstamp >= opt.t_range(1) & Tstamp <= opt.t_range(2);

nT = length(T);
ntri = size(raster, 3);
ncond = size(psth, 3);

score = nan(nT, ntri, ndim); % time x trial x dim
psth_score = nan(nT, ncond, ndim); 
for dim=1:ndim
    PSTH1 = raster(:, ind, :);
    sc = PSTH1(:,:)' * opt.PC_coef(:, dim);
    score(:,:,dim) = reshape(sc, [nT ntri]);
    
    PSTH1 = psth(:, ind, :);
    sc = PSTH1(:,:)' * opt.PC_coef(:, dim);
    psth_score(:,:,dim) = reshape(sc, [nT ncond]);
end

if ~isempty(opt.conv_kernel)
    for d=1:ndim
        for n=1:ntri
            score(:, n, d) = nanconv(score(:, n, d), opt.conv_kernel(:), 'same');
        end
        for n=1:ncond
            psth_score(:, n, d) = nanconv(psth_score(:, n, d), opt.conv_kernel(:), 'same');
        end
    end
end

%% compute flow field with sliding windows

% Compute time derivative: flow = T-1 x N x 2
flow = score(2:end,:,opt.flow_field_dim) - score(1:end-1,:,opt.flow_field_dim);  % (T-2) x N x 2

% For time points 1 to T-1
pos = score(1:end-1,:,opt.flow_field_dim);  % position at t
slice_pos = score(1:end-1,:,opt.slice_dim);
pos = reshape(pos, [], 2);
flow = reshape(flow, [], 2);
slice_pos = slice_pos(:);

% Bin positions
binSize = (max(pos(:,1)) - min(pos(:,1)) + max(pos(:,2)) - min(pos(:,2)))/2 / opt.bin_num;
x_edges = min(pos(:,1)):binSize:max(pos(:,1)) + binSize;
y_edges = min(pos(:,2)):binSize:max(pos(:,2)) + binSize;

% Create sliding windows based on slice dimension
% Sort trials by slice dimension value
[sorted_slice_pos, sort_idx] = sort(slice_pos);
sorted_pos = pos(sort_idx, :);
sorted_flow = flow(sort_idx, :);

% Calculate window parameters
n_trials = length(sorted_slice_pos);
window_size_trials = round(opt.window_size * n_trials);
window_step_trials = round(opt.window_step * n_trials);

% Calculate number of windows
n_windows = floor((n_trials - window_size_trials) / window_step_trials) + 1;

% Initialize storage
[Xgrid, Ygrid] = meshgrid( ...
    x_edges(1:end-1) + binSize/2, ...
    y_edges(1:end-1) + binSize/2 );
U = zeros([size(Xgrid), n_windows]);
V = zeros([size(Ygrid), n_windows]);
Cnt = zeros([size(Ygrid), n_windows]);
Z = zeros(n_windows, 1); % average slice position for each window

% Process each sliding window
for w = 1:n_windows
    % Calculate window boundaries
    start_idx = (w-1) * window_step_trials + 1;
    end_idx = start_idx + window_size_trials - 1;
    
    % Extract data for this window
    window_pos = sorted_pos(start_idx:end_idx, :);
    window_flow = sorted_flow(start_idx:end_idx, :);
    window_slice_pos = sorted_slice_pos(start_idx:end_idx);
    
    % Store average slice position for this window
    Z(w) = mean(window_slice_pos);
    
    % Average flow per bin for this window
    for x = 1:size(Xgrid, 1)
        for y = 1:size(Xgrid, 2)
            cx = Xgrid(x,y);
            cy = Ygrid(x,y);

            % Define bin boundaries
            xlo = cx - binSize/2;
            xhi = cx + binSize/2;
            ylo = cy - binSize/2;
            yhi = cy + binSize/2;

            % Find points in bin
            idx = window_pos(:,1) >= xlo & window_pos(:,1) < xhi & ...
                  window_pos(:,2) >= ylo & window_pos(:,2) < yhi;
            
            Cnt(x,y,w) = sum(idx);
            if any(idx)
                % Average flow
                avg_flow = mean(window_flow(idx,:),1);
                U(x,y,w) = avg_flow(1);
                V(x,y,w) = avg_flow(2);
            end
        end
    end
end

data.Xgrid = Xgrid;
data.Ygrid = Ygrid;
data.U = U;
data.V = V;
data.Z = Z;
data.Cnt = Cnt;
data.mean_trajectory = psth_score(:, :, [opt.flow_field_dim, opt.slice_dim]);
data.window_info = struct('window_size', opt.window_size, ...
                         'window_step', opt.window_step, ...
                         'n_windows', n_windows, ...
                         'window_size_trials', window_size_trials, ...
                         'window_step_trials', window_step_trials);

end
