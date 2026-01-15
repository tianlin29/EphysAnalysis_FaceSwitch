function data = PopulationClassifier(fnd, ID, opt)

def.epoch = [1 2];
def.tstamp = {0:50:450, -250:50:50};
def.t_win = 100;
def.min_trial = 10; % min trial for each category
def.max_nan_rate = .5;
def.Kfold = 5; % to build training and testing set
def.lasso_kFold = 5; % to find the best lambda
def.glm_type = 'ordinary'; % ordinary, lasso
def.lambda_set = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 1, 10, 100];
def.sq_rt_spike = false;
def.refBeta = [];
def.standardize_spike = false;
def.balance_trial = true; % balance number of trials for two categories
def.save_trial_param = true;
opt = safeStructAssign(def, opt);

%% check data
if any(ID(:)~= 1 & ID(:)~= 2 & ~isnan(ID(:)))
    error('condition should be either 1, 2, or NaN');
end

%% preprocess spike data
fprintf('Preprocessing...\n');
nepoch = length(opt.epoch);
raster = cell(nepoch,1);

for e = 1:nepoch
    raster(e) = fnd.raster(opt.epoch(e)); % (unit, time, trial)
    % sq root transform if requested
    if opt.sq_rt_spike
        raster{e} = sqrt(raster{e});
    end
end

% standardize each unit if requested
if opt.standardize_spike
    for n = 1:size(raster{1},1) % loop through units
        v = cell2mat(cellfun(@(x) x(n,:)', raster, 'uni', 0));
        mv = nanmean(v);
        msd = nanstd(v);
        for e = 1:nepoch
            raster{e}(n,:,:) = (raster{e}(n,:,:) - mv)/msd;
        end
    end
end

%% main calculation
warning('off', 'stats:glmfit:PerfectSeparation');
warning('off', 'stats:glmfit:IterationLimit');
warning('off', 'stats:glmfit:IllConditioned');

fprintf('Running glm...\n');
[data_ind, LogOdds, Correct, Beta] = deal(cell(1, nepoch));

target = ID(1,:)'; % (trial, 1)
I_invalid_trial = isnan(target);
target(I_invalid_trial) = []; % remove NaN trials

ntrial = length(target);
nunit = size(fnd.data{1}, 1);

if opt.save_trial_param
    try
        morph_level = fnd.getp('morph_level');
    catch
        try
            morph_level = fnd.getp('morph');
        catch
            error('check morph');
        end
    end
    task_set = fnd.getp('task_set');
    targ_cor = fnd.getp('targ_cor');
    targ_cho = fnd.getp('targ_cho');

    param.morph_level = morph_level(1, ~I_invalid_trial)';
    param.task_set = task_set(1, ~I_invalid_trial)';
    param.targ_cor = targ_cor(1, ~I_invalid_trial)';
    param.targ_cho = targ_cho(1, ~I_invalid_trial)';
end

for e = 1:nepoch
    r = raster{e}; % (unit, time, trial)
    r(:,:,I_invalid_trial) = [];

    ntime_bin = length(opt.tstamp{e});
    data_ind{e} = cell(1, ntime_bin);
    LogOdds{e} = nan(ntrial, ntime_bin);
    Correct{e} = nan(ntrial, ntime_bin);
    Beta{e} = nan(nunit+1, ntime_bin);

    for t = 1:ntime_bin
        I_time = fnd.tstamp{e}>=(opt.tstamp{e}(t)-opt.t_win/2) & fnd.tstamp{e}<(opt.tstamp{e}(t)+opt.t_win/2);
        r_ = squeeze(nanmean(r(:,I_time,:),2))'; % (trial, unit)
        nan_percent = squeeze(mean(isnan(r(:,I_time,:)),2))';
        r_(nan_percent(:) > opt.max_nan_rate) = NaN;
        if isempty(opt.refBeta)
            data_ind{e}{t} = compute_neural_performance(r_, target, [], opt);
        else
            data_ind{e}{t} = compute_neural_performance(r_, target, opt.refBeta{e}{t}.ncv.bBeta, opt);
        end
        if isempty(data_ind{e}{t})
            continue;
        end
        LogOdds{e}(:,t) = data_ind{e}{t}.LogOdds; % (trial, 1); range=[-Inf, Inf]
        Correct{e}(:,t) = data_ind{e}{t}.correct; % (trial, 1); classification accuracy
        Beta{e}(1:end-1,t) = mean(cell2mat(data_ind{e}{t}.bBeta(:,1)'),2); % (unit+1, 1); weight and intercept
        Beta{e}(end,t) = mean(cell2mat(data_ind{e}{t}.bBeta(:,2)));
    end
end

warning('on', 'stats:glmfit:PerfectSeparation');
warning('on', 'stats:glmfit:IterationLimit');
warning('on', 'stats:glmfit:IllConditioned');

%% data
data.tstamp = opt.tstamp;
data.data_ind = data_ind;
data.LogOdds = LogOdds;
data.Correct = Correct;
data.targ_cor = target;
data.Beta = Beta;
data.param = param;

end



function data = compute_neural_performance(raster, target, refBeta, opt)

% remove trial with NaN in raster
target(any(isnan(raster),2)) = NaN;

% balance trials
ntrial = [sum(target==1), sum(target==2)];
if min(ntrial) < opt.min_trial
    data = [];
    return;
end

if opt.balance_trial && ntrial(1) ~= ntrial(2)
    tr = min(ntrial(1), ntrial(2));
    if ntrial(1) < ntrial(2)
        ind = crandsample(find(target==2), ntrial(2) - tr);
    else
        ind = crandsample(find(target==1), ntrial(1) - tr);
    end
    target(ind) = NaN;
    ntrial = [sum(target==1), sum(target==2)];
end

% get cross-validation index
idx_CV = nan(length(target),1);
for n = 1:2
    if exist('crossvalind', 'file') % in bioinformatics toolbox
        idx_CV_sub = crossvalind('Kfold', ntrial(n), opt.Kfold);
    else
        idx_CV_sub = repmat(1:opt.Kfold, ceil(ntrial(n)/opt.Kfold), 1)';
        idx_CV_sub = idx_CV_sub(1:ntrial(n));
        idx_CV_sub = idx_CV_sub(randperm(length(idx_CV_sub)));
    end
    idx_CV(target==n) = idx_CV_sub; % assign fold number to 2 targets seperately
end

%% classification
if isempty(refBeta)
    % run cross validation
    Beta = cell(opt.Kfold,1);
    bBeta = cell(opt.Kfold,2);
    FitInfo = cell(opt.Kfold,1);
    LogOdds = nan(length(target),1);

    for f = 1:opt.Kfold
        if opt.Kfold==1 % no CV
            I_train = ~isnan(idx_CV);
            I_test = I_train;
        else
            I_train = idx_CV~=f & ~isnan(idx_CV);
            I_test = idx_CV==f;
        end

        switch opt.glm_type
            case 'lasso'
                [Beta{f}, FitInfo{f}] = lassoglm(raster(I_train, :), target(I_train)==2, 'binomial', ...
                    'link', 'logit', 'Alpha', 1, 'Lambda', opt.lambda_set, 'CV', opt.lasso_kFold);
                minB = Beta{f}(:, FitInfo{f}.IndexMinDeviance); % (unit, 1); weight
                minB0 = FitInfo{f}.Intercept(FitInfo{f}.IndexMinDeviance); % (1, 1); intercept
            case 'ordinary'
                Beta{f} = glmfit(raster(I_train, :), target(I_train)==2, 'binomial', ...
                    'link', 'logit');
                minB = Beta{f}(2:end);
                minB0 = Beta{f}(1);
        end

        bBeta{f,1} = minB;
        bBeta{f,2} = minB0;
        LogOdds(I_test) = minB0 + raster(I_test,:) * minB; % (trial, 1); range=[-Inf, Inf], project to [0, 1] by sigmoid
    end

    data.Beta = Beta; % {fold, 1}; (unit, lambda)
    data.bBeta = bBeta; % {fold, 2}; (unit, 1), (1, 1); weight and intercept
    data.FitInfo = FitInfo;
    data.LogOdds = LogOdds; % (trial, 1); range=[-Inf, Inf]
    data.Prob1 = exp(data.LogOdds)./(1+exp(data.LogOdds)); % (trial, 1); range=[0, 1]; classification result
    data.correct = double((data.LogOdds>0)==(target==2)); % (trial, 1); classification accuracy
    data.correct(isnan(target)) = NaN;
    data.target = target; % (trial, 1); true category

    % don't run cross validation
    I_valid_trial = ~isnan(target);
    switch opt.glm_type
        case 'lasso'
            [BetaAll, FitInfoAll] = lassoglm(raster(I_valid_trial,:), target(I_valid_trial)==2, 'binomial','link','logit','Alpha',1, 'Lambda',opt.lambda_set,'CV', opt.lasso_kFold);
            minB = BetaAll(:,FitInfoAll.IndexMinDeviance);
            minB0 = FitInfoAll.Intercept(FitInfoAll.IndexMinDeviance);
        case 'ordinary'
            BetaAll = glmfit(raster(I_valid_trial,:), target(I_valid_trial)==2, 'binomial','link','logit');
            FitInfoAll = [];
            minB = BetaAll{k}(2:end);
            minB0 = BetaAll{k}(1);
    end

    data.ncv.LogOdds = minB0 + raster * minB;
    data.ncv.bBeta = {minB, minB0};
    data.ncv.Beta = BetaAll;
    data.ncv.FitInfo = FitInfoAll;
    data.ncv.Prob1 = exp(data.ncv.LogOdds)./(1+exp(data.ncv.LogOdds));
    data.ncv.correct = double((data.ncv.LogOdds>0)==(target==2));
    data.ncv.correct(isnan(target)) = NaN;
else
    LogOdds = refBeta{2} + raster * refBeta{1}; % (trial, 1); range=[-Inf, Inf]
    LogOdds(isnan(idx_CV),:) = NaN;
    data.LogOdds = LogOdds;
    data.Prob1 = exp(data.LogOdds)./(1+exp(data.LogOdds)); % (trial, 1); range=[0, 1]; classification result
    data.correct = double((data.LogOdds>0)==(target==2)); % (trial, 1); classification accuracy
    data.correct(isnan(target)) = NaN;
    data.bBeta = refBeta;
end

end

