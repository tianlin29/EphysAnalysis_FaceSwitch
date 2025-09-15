function [coef, detPSTH] = getPC_coef(tstamp, data, opt)

def.epoch = 1;
def.PC_range = [250 600];
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
detPSTH = nanmean(PSTH,3); % (unit, time)
PSTH = bsxfun(@minus, PSTH, detPSTH);

% PC range
ind = Tstamp>=opt.PC_range(1) & Tstamp<=opt.PC_range(2);
pcPSTH = PSTH(:, ind, :); % (unit, time, cond)

%% run PCA
coef = pca(pcPSTH(:,:)'); % coeff ..(feature, pc)

end


