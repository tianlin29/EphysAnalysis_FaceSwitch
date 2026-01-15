function [coeff, psth_mn] = getPC_coef(tstamp, psth, opt)

def.epoch = 1;
def.PC_range = [250 600];
opt = safeStructAssign(def, opt);

% extract time epoch
psth = psth{opt.epoch};
tstamp = tstamp{opt.epoch};

% mean imputation to nan
if sum(isnan(psth(:)))>0
    fill_val = bsxfun(@times, isnan(psth), nanmean(psth(:,:),2));
    psth(isnan(psth(:))) = fill_val(isnan(psth(:)));
end

% detrend
psth_mn = nanmean(psth, 3); % (unit, time)
psth = bsxfun(@minus, psth, psth_mn);

% PC time range
I = tstamp>=opt.PC_range(1) & tstamp<=opt.PC_range(2);
psth = psth(:, I, :); % (unit, time, cond)

% run PCA
coeff = pca(psth(:,:)'); % coeff ..(unit, pc)

end


