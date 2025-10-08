function data = gen_popresp_dPCA_axis(fnd, ID, opt)
def.epoch = 1;
def.target = {};
def.coefficient_data = [];
def.detrend = false;

opt = safeStructAssign(def, opt);

if isempty(opt.coefficient_data)
    error('you should put coefficient_data');
end

coef = opt.coefficient_data;


%% preprocess

psth = fnd.PSTH(ID, [], [], [], opt.epoch); % no smoothing
psth = psth{opt.epoch}; % cell x time x cond
raster = fnd.raster(opt.epoch); raster = raster{1}*1000; % (unit, time, trial)

if any(coef.exclude_unit)
    psth(coef.exclude_unit,:,:) = [];
    raster(coef.exclude_unit,:,:) = [];
end

if opt.detrend
    psth = psth - nanmean(psth, 3);
    raster = raster - nanmean(psth, 3);
end

psth = bsxfun(@rdivide, bsxfun(@minus, psth, coef.MU), coef.SIGMA);
raster = bsxfun(@rdivide, bsxfun(@minus, raster, coef.MU), coef.SIGMA);

%% apply dPCA weight

ntarg = length(opt.target);
tdim = nan(ntarg,1);
for n=1:ntarg
    I = find(strcmp(opt.target{n}{1}, coef.margNames));
    if isempty(I)
        error('No such dimension: %s', opt.target{n}{1});
    end
    tmpdim = find(coef.whichMarg == I, opt.target{n}{2});
    tdim(n) = tmpdim(end);
end

dpc = nan(ntarg, size(psth,2), size(psth,3));
dpc_trial = nan(ntarg, size(raster,2), size(raster,3));

for n = 1:ntarg
    dpc(n,:) = psth(:,:)' * coef.W(:, tdim(n));
    dpc_trial(n,:) = raster(:,:)' * coef.W(:, tdim(n));
end

%% output

data.dpc = dpc;
data.dpc_trial = dpc_trial;
data.target = opt.target;
data.tdim = tdim;
data.tstamp = fnd.tstamp{opt.epoch};
data.cutoff = fnd.cutoff{opt.epoch};















