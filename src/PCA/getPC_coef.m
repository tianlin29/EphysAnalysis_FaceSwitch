function [coef, detPSTH] = getPC_coef(tstamp, data, opt)

def.epoch = 1;
def.PC_range = [250 600];
def.t_range = [0 600];
def.anchor_t = 0:50:600;
def.interval = 10;
def.dim = [1 2 3];
def.interp_method = {'polyfit',3};

def.xlim = [];
def.ylim = [];
def.zlim = [];
def.view_angle = [48, 28];
def.plot = [];

opt = safeStructAssign(def, opt);

PSTH = data{opt.epoch};
Tstamp = tstamp{opt.epoch};


%% preprocess
% mean imputation to nan
if sum(isnan(PSTH(:)))>0
    fill_val = bsxfun(@times, isnan(PSTH), nanmean(PSTH(:,:),2));
    PSTH(isnan(PSTH(:))) = fill_val(isnan(PSTH(:)));
end

% detrend
detPSTH = nanmean(PSTH,3);
PSTH = bsxfun(@minus, PSTH, detPSTH);

% PC range
ind = Tstamp >= opt.PC_range(1) & Tstamp <= opt.PC_range(2);
pcPSTH = PSTH(:, ind, :); % cell x time x cond

%% run PCA
coef = pca(pcPSTH(:,:)');


