function [fh, output_data] = dPCA_stim_choice(raster, tstamp, opt)

def.param_combination = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
def.param_name = {'Difficulty', 'Choice', 'Condition-independent', 'D/C Interaction'};
def.t_range = []; % ms; the range of time for dPCA
def.tbin = 10; % ms
def.detrend = false;
def.zscoring = false;
def.regularization = [];
def.coefficient_data = [];
def.max_dpc = 20; % number of dPCA components
def.show_figure = true;
opt = safeStructAssign(def, opt);

if any(raster(:)~=0 & raster(:)~=1 & ~isnan(raster(:)))
    error('Input raster is supposed to be a raster data.');
end

%% setup
if ~isempty(opt.regularization) && any(opt.regularization.exclude_unit) % if regularizing
    raster(opt.regularization.exclude_unit,:,:,:) = [];
    exclude_unit = opt.regularization.exclude_unit;
else
    exclude_unit = false(size(raster,1),1);
end

% convert raster to firing rate
dtstamp = (opt.t_range(1):opt.tbin:opt.t_range(2)) + opt.tbin/2;
ntbin = length(dtstamp);
firingRates = nan(size(raster,1), size(raster, 2), ntbin, size(raster, 4)); % (unit, trial, time bin, condition)
for t = 1:ntbin
    tind = (tstamp >= opt.t_range(1) + (t-1)*opt.tbin) & (tstamp < opt.t_range(1) + t*opt.tbin);
    firingRates(:,:,t,:) = mean(raster(:,:,tind,:),3) * 1e3; % not use nanmean so any nan will make the data point nan
end

% reorganize FR
firingRates = permute(firingRates, [1 4 3 2]); % (unit, trial, time, condition) -> (unit, condition, time, trial)
[nunit, ncond, ntime, ntrial] = size(firingRates);
firingRates = reshape(firingRates, [nunit, ncond/2, 2, ntime, ntrial]); % (unit, stimulus, 2 choices, time, trial)

% detrend
if opt.detrend
    mTrend = nanmean(nanmean(nanmean(firingRates,5),3),2); % (unit, 1, 1, time)
    firingRates = bsxfun(@minus, firingRates, mTrend);
end

firingRatesAverage = nanmean(firingRates,5);
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of time-points (note that all the trials should have the
% same length in time!)

%% centering
% centering
X = firingRatesAverage(:,:);
MU = mean(X,2); % (unit, 1)
if opt.zscoring
    SIGMA = std(X,[],2);
else
    SIGMA = ones(size(X,1),1);
end
firingRatesAverage = bsxfun(@rdivide, bsxfun(@minus, firingRatesAverage, MU), SIGMA);

%% run dPCA
% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

if isempty(opt.coefficient_data)
    % run dpca
    if ~isempty(opt.regularization)
        lambda = opt.regularization.optimalLambda;
    else
        lambda = 0;
    end

    [W,V,whichMarg] = dpca(firingRatesAverage, opt.max_dpc, ...
        'combinedParams', opt.param_combination, 'lambda', lambda);

    explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
        'combinedParams', opt.param_combination);
else
    W = opt.coefficient_data.W;
    V = opt.coefficient_data.V;
    whichMarg = opt.coefficient_data.whichMarg;
    explVar = opt.coefficient_data.explVar;
end

fh = nan;
if opt.show_figure
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
    dpca_plot(firingRatesAverage, W, V, @dpca_plot_custom, ...
        'explainedVar', explVar, ...
        'marginalizationNames', opt.param_name, ...
        'marginalizationColours', margColours, ...
        'whichMarg', whichMarg,                 ...
        'time', dtstamp,                        ...
        'timeEvents', 0,               ...
        'timeMarginalization', 3, ...
        'legendSubplot', 16);
    fh = gcf;
end

%% save data
if nargout > 1
    output_data.MU = MU;
    output_data.SIGMA = SIGMA;
    output_data.W = W;
    output_data.V = V;
    output_data.whichMarg = whichMarg;
    output_data.explVar = explVar;
    output_data.combinedParams = opt.param_combination;
    output_data.max_dpc = opt.max_dpc;
    output_data.margNames = opt.param_name;

    output_data.exclude_unit = exclude_unit;
    output_data.tstamp = dtstamp;
    if opt.detrend
        output_data.mTrend = mTrend;
    end

    si = size(firingRatesAverage);
    score = (firingRatesAverage(:,:)' * W)';
    output_data.score = reshape(score, size(W,2), si(2), si(3), si(4));
    % dPCA dim x stim x choice x time
end

end




