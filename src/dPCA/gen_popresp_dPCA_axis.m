function data = gen_popresp_dPCA_axis(fnd, ID, opt)
def.epoch = 1;
def.target = {};
def.coefficient_data = [];

opt = safeStructAssign(def, opt);

if isempty(opt.coefficient_data)
    error('you should put coefficient_data');
end

coef = opt.coefficient_data;


%% preprocess

psth = fnd.PSTH(ID, [], [], [], opt.epoch); % no smoothing
psth = psth{opt.epoch}; % cell x time x cond

if any(coef.exclude_unit)
    psth(coef.exclude_unit,:,:) = [];
end

if isfield(coef, 'mTrend')
    trend = nanmean(psth, 3);
    psth = psth - trend;
end

psth = bsxfun(@rdivide, bsxfun(@minus, psth, coef.MU), coef.SIGMA);


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

for n=1:ntarg
    dpc(n,:) = psth(:,:)' * coef.W(:, tdim(n));
end

%% output

data.dpc = dpc;
data.target = opt.target;
data.tdim = tdim;
data.tstamp = fnd.tstamp{opt.epoch};
data.cutoff = fnd.cutoff{opt.epoch};















