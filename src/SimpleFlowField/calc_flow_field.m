function data = calc_flow_field(tstamp, raster, ID, opt)
% function data = calc_flow_field(Tstamp, raster, ID, opt)
% Tstamp = 1 x time
% raster = cell x time x trial
% ID = cell x trial

def.PC_coef = [];
def.detPSTH = [];

def.t_range = [0 600];
def.dim = [1 2];
def.conv_kernel = [];
def.bin_num = 20;

opt = safeStructAssign(def, opt);

if isempty(opt.PC_coef)
    error('opt.PC_coef should be filled');
end

ndim = max(opt.dim); % 2pcs or 3 pcs, choose pc1-2 or pc1-3 later

%% smooth raster
raster = raster * 1e3; % (unit, time, trial)

% detrend
if ~isempty(opt.detPSTH)
    raster = raster - opt.detPSTH;
end

% average raster across trials
ncond = nanmax(ID(:));
psth = nan(size(raster,1), size(raster,2), ncond);
for c = 1:ncond
    psth(:,:,c) = nanmean(raster(:,:,ID(1,:)==c),3);
end

%% conver to PCA
I_valid_trial = ~isnan(ID(1,:));
fprintf('remove %d trials with nan in ID.\n', sum(isnan(ID(1,:))))
ntrial = sum(I_valid_trial);

T = opt.t_range(1):opt.t_range(2);
I_valid_time = tstamp>=opt.t_range(1) & tstamp<=opt.t_range(2);
ntime = length(T);

ncond = size(psth, 3);

score = nan(ntime, ntrial, ndim); % raster projects to PC space (time, trial, dim)
psth_score = nan(ntime, ncond, ndim); 
for d = 1:ndim
    psth_ = raster(:, I_valid_time, I_valid_trial);
    score_ = psth_(:,:)' * opt.PC_coef(:,d);
    score(:,:,d) = reshape(score_, [ntime ntrial]);
    
    psth_ = psth(:, I_valid_time, :);
    score_ = psth_(:,:)' * opt.PC_coef(:,d);
    psth_score(:,:,d) = reshape(score_, [ntime ncond]);
end

if ~isempty(opt.conv_kernel)
    for d = 1:ndim
        for tr = 1:ntrial
            score(:,tr,d) = nanconv(score(:,tr,d), opt.conv_kernel(:), 'same'); % smooth along time dimension
        end
        for c = 1:ncond
            psth_score(:,c,d) = nanconv(psth_score(:,c,d), opt.conv_kernel(:), 'same');
        end
    end
end

%% compute flow field
% compute time derivative
flow = score(2:end,:,opt.dim) - score(1:end-1,:,opt.dim);  % velocity at t, (time-1, trial, dim)

% for time points 1 to t-1
pos = score(1:end-1,:,opt.dim);  % position at t
pos = reshape(pos, [], 2); % (time*trial, dim)
flow = reshape(flow, [], 2); % (time*trial, dim)

% bin positions
binSize = (max(pos(:,1)) - min(pos(:,1)) + max(pos(:,2)) - min(pos(:,2)))/2 / opt.bin_num;
x_edges = min(pos(:,1)) : binSize : (max(pos(:,1)) + binSize);
y_edges = min(pos(:,2)) : binSize : (max(pos(:,2)) + binSize);

% initialize storage
[Xgrid, Ygrid] = meshgrid( ...
    x_edges(1:end-1) + binSize/2, ...
    y_edges(1:end-1) + binSize/2);
U = zeros(size(Xgrid));
V = zeros(size(Ygrid));
Cnt = zeros(size(Ygrid));

% average flow per bin
for xi = 1:numel(Xgrid)  
    cx = Xgrid(xi);
    cy = Ygrid(xi);
    
    % Define bin boundaries
    xlo = cx - binSize/2;
    xhi = cx + binSize/2;
    ylo = cy - binSize/2;
    yhi = cy + binSize/2;
    
    % Find points in bin
    idx = pos(:,1) >= xlo & pos(:,1) < xhi & ...
          pos(:,2) >= ylo & pos(:,2) < yhi;
    Cnt(xi) = sum(idx);
    if any(idx)
        % Average flow
        avg_flow = mean(flow(idx,:),1); % average across trials
        U(xi) = avg_flow(1);
        V(xi) = avg_flow(2);
    end
end

data.Xgrid = Xgrid;
data.Ygrid = Ygrid;
data.U = U;
data.V = V;
data.Cnt = Cnt;
data.mean_trajectory = psth_score(:,:,opt.dim);

end


