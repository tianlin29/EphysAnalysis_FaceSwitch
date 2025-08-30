function data = calc_input_flow_field(Tstamp, raster, ID, input, opt)
% function data = calc_input_flow_field(Tstamp, raster, ID, input, opt)
% Tstamp = 1 x time
% raster = cell x time x trial
% ID = cell x trial
% input = time x trial x param

def.PC_coef = [];
def.detPSTH = [];

def.t_range = [0 600];
def.dim = [1 2];
def.conv_kernel = [];
def.bin_num = 20;
def.min_sample = 100;

def.xlim = [];
def.ylim = [];
def.zlim = [];
def.view_angle = [48, 28];
def.color = [1 0 0; 0 0 1];

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


%% conver to PCA

valid_tri = ~isnan(ID(1,:));

T = opt.t_range(1):opt.t_range(2);
ind = Tstamp >= opt.t_range(1) & Tstamp <= opt.t_range(2);

input = input(ind, :, :);

nT = length(T);

ntri = sum(valid_tri);
ncond = size(psth, 3);

score = nan(nT, ntri, ndim); % time x trial x dim
psth_score = nan(nT, ncond, ndim); 
for dim=1:ndim
    PSTH1 = raster(:, ind, valid_tri);
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


%% compute effective dynamics (input drive + internal dynamics)

% Compute time derivative: flow = T-1 x N x 2
flow = score(2:end, :, opt.dim) - score(1:end-1, :, opt.dim);  % (T-2) x N x 2

% For time points 1 to T-1
pos = score(1:end-1, :, opt.dim);  % position at t
pos = reshape(pos, [], 2);
flow = reshape(flow, [], 2);
input = reshape(input(1:end-1, valid_tri, :), [], 2);

% Bin positions
binSize = (max(pos(:,1)) - min(pos(:,1)) + max(pos(:,2)) - min(pos(:,2)))/2 / opt.bin_num;
x_edges = min(pos(:,1)):binSize:max(pos(:,1)) + binSize;
y_edges = min(pos(:,2)):binSize:max(pos(:,2)) + binSize;

% Initialize storage
[Xgrid, Ygrid] = meshgrid( ...
    x_edges(1:end-1) + binSize/2, ...
    y_edges(1:end-1) + binSize/2 );

pos_idx = nan(size(flow,1), 1);

% Average flow per bin
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
    if any(idx)
        pos_idx(idx) = xi;
    end
end

nind = any(isnan(flow), 2) | any(isnan(input), 2);

flow(nind,:) = [];
input(nind,:) = [];
pos_idx(nind) = [];

%% dissociate input drive and internal dynamics

% flow = dynamics + input * W

W_init = input \ flow;

options = optimset('Display', 'iter', 'MaxFunEvals', 1e3, 'MaxIter', 1e3, 'TolX', 0.001, 'TolFun', 0.001);

W = fminsearch(@dynamic_fun, W_init, options);

    function err = dynamic_fun(W1)
        D1 = get_dynamics(W1);
        flow_est = D1 + input * W1;
        err = sqrt(sum((flow(:) - flow_est(:)).^2));
    end

    function D = get_dynamics(W1)
        D = flow - input * W1;
        for n1 = 1:numel(Xgrid)
            D(pos_idx == n1, :) = ones(sum(pos_idx == n1), 1) * mean(D(pos_idx == n1, :), 1);
        end
    end

D = get_dynamics(W);
U = zeros(size(Xgrid));
V = zeros(size(Ygrid));
Cnt = zeros(size(Ygrid));
for xi = 1:numel(Xgrid)
    U(xi) = mean(D(pos_idx == xi, 1));
    V(xi) = mean(D(pos_idx == xi, 2));
    Cnt(xi) = sum(pos_idx == xi);
end

data.Xgrid = Xgrid;
data.Ygrid = Ygrid;
data.U = U;
data.V = V;
data.Cnt = Cnt;
data.W = W;
data.mean_trajectory = psth_score(:, :, opt.dim);

end







