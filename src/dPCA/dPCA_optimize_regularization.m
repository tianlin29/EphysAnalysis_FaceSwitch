function [fh, output_data] = dPCA_optimize_regularization(raster, tstamp, opt)

def.param_combination = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
def.param_name = {'Difficulty', 'Choice', 'Condition-independent', 'D/C Interaction'};
def.t_range = []; % ms; the range of time that would be averaged
def.tbin = 10; % ms; average raster into FR within each time bin
def.detrend = false; % detrend for each unit and timepoint
def.zscoring = false; % perform zscoring or not during centering for each unit
def.numRep = 10; % how many cross-validation iterations to perform in the dPCA program
def.minTrial = 2; % remove units that do not have minimum # of trials in any of cond
opt = safeStructAssign(def, opt);

if any(raster(:)~=0 & raster(:)~=1 & ~isnan(raster(:)))
    error('Input raster is supposed to be a raster data.');
end

%% setup
% convert raster to firing rate
ntbin = length(opt.t_range(1):opt.tbin:opt.t_range(2));
firingRates = nan(size(raster,1), size(raster,2), ntbin, size(raster,4)); % (unit, trial, time, condition)
for t = 1:ntbin
    tind = (tstamp >= opt.t_range(1) + (t-1)*opt.tbin) & (tstamp < opt.t_range(1) + t*opt.tbin);
    firingRates(:,:,t,:) = mean(raster(:,:,tind,:),3) * 1e3; % not use nanmean so any nan will make the data point nan
end

% reorganize FR
firingRates = permute(firingRates, [1 4 3 2]); % (unit, trial, time, condition) -> (unit, condition, time, trial)
[nunit, ncond, ntime, ntrial] = size(firingRates);
firingRates = reshape(firingRates, [nunit, ncond/2, 2, ntime, ntrial]); % (unit, stimulus, 2 choices, time, trial)

% detrend for each unit and timepoint
if opt.detrend
    mTrend = nanmean(nanmean(nanmean(firingRates,5),3),2); % (unit, 1, 1, time, 1)
    firingRates = bsxfun(@minus, firingRates, mTrend);
end

%% validation of data
[nunit, nstim, ndec, ~, ~] = size(firingRates);

% count number of trials, remove trials with nan within the range
ntri = nan(nunit, nstim, ndec);
for n = 1:nunit
    for s = 1:nstim
        for c = 1:ndec
            nanTr = any(isnan(squeeze(firingRates(n,s,c,:,:))), 1);
            validTr = ~nanTr;
            firingRates(n,s,c,:,nanTr) = NaN;
            firingRates(n,s,c,:,:) = cat(5, firingRates(n,s,c,:,validTr), firingRates(n,s,c,:,nanTr)); % arrange dim 5 as valid trials + nan trials
            ntri(n,s,c) = sum(validTr);
        end
    end
end

% remove units that do not have minimum # of trials in any of cond
exclude_unit = any(reshape(ntri, [nunit, nstim * ndec]) < opt.minTrial, 2);
if any(exclude_unit)
    fprintf('%d/%d units removed due to lack of trials\n', sum(exclude_unit), length(exclude_unit));
    ntri(exclude_unit,:,:) = [];
    firingRates(exclude_unit,:,:,:,:) = [];
end

%% centering
firingRates = firingRates(:,:,:,:,1:max(ntri(:)));
firingRatesAverage = nanmean(firingRates, 5); % (unit, stimulus, 2 choices, time)

% centering for each unit
X = firingRatesAverage(:,:); % (unit, rest)
MU = mean(X, 2); % (unit, 1)
if opt.zscoring
    SIGMA = std(X, [], 2);
else
    SIGMA = ones(size(X,1), 1);
end
firingRates = bsxfun(@rdivide, bsxfun(@minus, firingRates, MU), SIGMA); % FR = (FR-mu)/sigma
firingRatesAverage = bsxfun(@rdivide, bsxfun(@minus, firingRatesAverage, MU), SIGMA);

%% find optimal lambda
optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, ntri, ...
    'combinedParams', opt.param_combination, ...
    'numRep', opt.numRep);  % increase this number to ~10 for better accuracy
fh = gcf;

fprintf('optimal lambda = %g\n', optimalLambda);

%% save data
output_data.optimalLambda = optimalLambda;
output_data.exclude_unit = exclude_unit;

end




