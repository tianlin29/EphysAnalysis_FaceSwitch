run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, threeExemplar, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);

%% classifier
% n_files = 2;
n_cond = 3; n_repeat = 10;
for r = 1:n_repeat
    stat = cell(n_files, n_cond);
    for n = 15
        fprintf('repeat%d session%d\n', r, n)
        fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
        fnd = fnd.extract_epoch(2);
        raster = fnd.raster(1); raster = raster{1}; % (unit, time, trial)

        targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)';
        task_set = fnd.getp('task_set'); task_set = task_set(1,:)';

        clear opt;
        opt.t_win = [100, 200, 300, 400, 500, 600, 700];

        % 1 to 1
        opt.Kfold = 5;
        target = targ_cor; target(task_set~=1) = NaN;
        stat{n,1} = get_classifier(fnd, raster, target, opt);

        % 2 to 2
        opt.Kfold = 5;
        target = targ_cor; target(task_set~=2) = NaN;
        stat{n,2} = get_classifier(fnd, raster, target, opt);

        % 2 to 1
        opt.Kfold = 1;
        target = targ_cor; target(task_set~=1) = NaN;
        stat_tmp = get_classifier(fnd, raster, target, opt);

        opt.coef_pretrained = stat_tmp.coef;
        target = targ_cor; target(task_set~=2) = NaN;
        stat{n,3} = get_classifier(fnd, raster, target, opt);
    end
    save(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r)), 'stat')
end


%% merge data
n_files = 2; n_cond = 3; n_repeat = 2;
[session, repeat, cond, correct, morph_level] = deal([]);
for r = 1:n_repeat
    stat = load(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r))).stat;
    for n = 1:n_files
        for c = 1:n_cond
            correct_tmp = stat{n,c}.correct;
            morph_level_tmp = stat{n,c}.param.morph_level;
            n_trial = length(morph_level_tmp);

            session = cat(1, session, n*ones(n_trial, 1));
            repeat = cat(1, repeat, r*ones(n_trial, 1));
            cond = cat(1, cond, c*ones(n_trial, 1));
            correct = cat(1, correct, correct_tmp);
            morph_level = cat(1, morph_level, morph_level_tmp);
        end
    end
end

t_win = stat{1,1}.opt.t_win;
n_win = length(t_win)-1;

[data.lucoh, data.p, data.pse, data.fcoh, data.fresp, data.thres] = deal(cell(n_files, n_cond, n_win, n_repeat));
for r = 1:n_repeat
    for n = 1:n_files
        for c = 1:n_cond
            for t = 1:n_win
                I = repeat==r & session==n & cond==c;
                correct_tmp = correct(I,t);
                morph_level_tmp = morph_level(I);

                abs_coh = abs(morph_level_tmp);
                ucoh = unique(abs_coh);

                lucoh = log(ucoh);
                mincoh = ucoh(1);
                if mincoh == 0
                    mincoh = ucoh(2)/ucoh(3) * ucoh(2);
                    lucoh(1) = log(mincoh);
                end

                fcoh = linspace(0.03, max(abs_coh), 100)';

                [p, pse] = calcGroupMean(correct_tmp, abs_coh, ucoh, 'binary');

                [alpha, ~, stat_glmfit] = glmfit(abs_coh, correct_tmp, 'binomial', 'link', 'logit', 'Constant', 'off'); % logit(P) = a0 + a1*coh
                fresp = glmval(alpha, fcoh, 'logit', 'Constant', 'off');
                thres = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit', 'Constant', 'off') - 0.816), .2);

                data.lucoh{n,c,t,r} = lucoh;
                data.p{n,c,t,r} = p;
                data.pse{n,c,t,r} = pse;
                data.fcoh{n,c,t,r} = fcoh;
                data.fresp{n,c,t,r} = fresp;
                data.thres{n,c,t,r} = thres;
            end
        end
    end
end

%% plot threshold of all sessions
opt.color = [0 0 0; 0 0 1; 1 0 0];

thres = cell2mat(data.thres);
thres_mn = mean(thres, 4);
thres_se = std(thres, [], 4) ./ sqrt(n_repeat);

figure('Position', [50 50 800 800]);
for n = 1:n_files
    subplot(5,5,n); hold on
    for c = 1:n_cond
        plot(1:n_win, squeeze(thres_mn(n,c,:)), '.-', 'markers', 7, 'Color', opt.color(c,:));
        cerrorbar(1:n_win, squeeze(thres_mn(n,c,:)), squeeze(thres_se(n,c,:)), 'Color', opt.color(c,:));
    end
    xlabel('#Time window')
    ylabel('Threshold')
    title(sprintf('session %d', n))
end

%% plot threshold of all time windows
opt.color = [0 0 0; 0 0 1; 1 0 0];

thres = cell2mat(data.thres);
thres_mn = mean(thres, 4);
thres_se = std(thres, [], 4) ./ sqrt(n_repeat);

figure('Position', [50 50 800 800]);
for t = 1:n_win
    subplot(5,5,t); hold on
    for c = 1:n_cond
        plot(1:n_files, thres_mn(:,c,t), '.-', 'markers', 7, 'Color', opt.color(c,:));
        cerrorbar(1:n_files, thres_mn(:,c,t), thres_se(:,c,t), 'Color', opt.color(c,:));
    end
    xlabel('#Session')
    ylabel('Threshold')
    title(sprintf('%d-%d ms', t_win(t), t_win(t+1)))
end

%% plot figures of all sessions
opt.color = [0 0 0; 0 0 1; 1 0 0];

t = 3; r = 1;

figure('Position', [50 50 800 800]);
for n = 1:n_files
    subplot(5,5,n); hold on
    for c = 1:n_cond
        lucoh = data.lucoh{n,c,t,r};
        p = data.p{n,c,t,r};
        pse = data.pse{n,c,t,r};
        fcoh = data.fcoh{n,c,t,r};
        fresp = data.fresp{n,c,t,r};

        plot(log(fcoh), fresp, 'Color', opt.color(c,:));
        plot(lucoh, p, '.', 'markers', 7, 'Color', opt.color(c,:));
        cerrorbar(lucoh, p, pse, 'Color', opt.color(c,:));
    end

    opt.xlim = log([2.8 100]/100);
    opt.xtick = log([3 6 12 24 48 96]'/100);
    opt.xticklabel = {'0', '6', '', '24', '', '96'};
    format_panel(gca, ...
        'ylim', [.5 1], 'ytick', .5:.125:1, ...
        'xlim', opt.xlim, 'xtick', opt.xtick, 'xticklabel', opt.xticklabel, ...
        'xlabel', {'Stimulus strength','(% Morph)'}, 'ylabel', 'P(correct)');
    title(sprintf('session %d', n))
end

%% plot figures of all time windows
opt.color = [0 0 0; 0 0 1; 1 0 0];

n = 1; r = 1;

figure('Position', [50 50 800 800]);
for t = 1:n_win
    subplot(5,5,t); hold on
    for c = 1:n_cond
        lucoh = data.lucoh{n,c,t,r};
        p = data.p{n,c,t,r};
        pse = data.pse{n,c,t,r};
        fcoh = data.fcoh{n,c,t,r};
        fresp = data.fresp{n,c,t,r};

        plot(log(fcoh), fresp, 'Color', opt.color(c,:));
        plot(lucoh, p, '.', 'markers', 7, 'Color', opt.color(c,:));
        cerrorbar(lucoh, p, pse, 'Color', opt.color(c,:));
    end

    opt.xlim = log([2.8 100]/100);
    opt.xtick = log([3 6 12 24 48 96]'/100);
    opt.xticklabel = {'0', '6', '', '24', '', '96'};
    format_panel(gca, ...
        'ylim', [.5 1], 'ytick', .5:.125:1, ...
        'xlim', opt.xlim, 'xtick', opt.xtick, 'xticklabel', opt.xticklabel, ...
        'xlabel', {'Stimulus strength','(% Morph)'}, 'ylabel', 'P(correct)');
    title(sprintf('%d-%d ms', t_win(t), t_win(t+1)))
end

%% function
function stat = get_classifier(fnd, raster, target, opt)
def.Kfold = 5;
def.lambda_set = [0.005, 0.01, 0.015, 0.02, 0.025];
def.lasso_kFold = 5;
def.t_win = [100, 200, 300, 400, 500, 600, 700];
def.coef_pretrained = [];
opt = safeStructAssign(def, opt);

% get data
tstamp = fnd.tstamp{1};

% balance number of trials
ntrial = [sum(target==1), sum(target==2)];
if ntrial(1) ~= ntrial(2)
    tr = min(ntrial(1), ntrial(2));
    if ntrial(1) < ntrial(2)
        idx = crandsample(find(target==2), ntrial(2)-tr);
    else
        idx = crandsample(find(target==1), ntrial(1)-tr);
    end
else
    idx = [];
end
valid = ~ismember((1:length(target))', idx) & ~isnan(target);
raster = raster(:,:,valid);
target = target(valid);
ntrial = [sum(target==1), sum(target==2)];

% save parameters
targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,valid)';
targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,valid)';
task_set = fnd.getp('task_set'); task_set = task_set(1,valid)';
morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,valid)';
param = struct('targ_cor', targ_cor, 'targ_cho', targ_cho, 'task_set', task_set, 'morph_level', morph_level);

% CV
idx_CV = nan(length(target),1);
for n = 1:2
    idx_CV_sub = crossvalind('Kfold', ntrial(n), opt.Kfold);
    idx_CV(target==n) = idx_CV_sub; % assign fold number to 2 IDs seperately
end

% lassoglm
if isempty(opt.coef_pretrained) % perform lassoglm
    n_win = length(opt.t_win)-1;
    coef = nan(size(raster,1)+1, n_win, opt.Kfold);
    correct = nan(length(target), n_win);
    for t = 1:n_win
        r_tmp = squeeze(mean(raster(:, tstamp>=opt.t_win(t) & tstamp<=opt.t_win(t+1), :), 2))'; % (trial, unit)
        for f = 1:opt.Kfold
            if opt.Kfold==1
                I_train = true(size(target));
                I_test = true(size(target));
            else
                I_train = idx_CV~=f;
                I_test = idx_CV==f;
            end

            [B, FitInfo] = lassoglm(r_tmp(I_train,:), target(I_train)==2, 'binomial', ...
                'link', 'logit', 'Alpha', 1, 'Lambda',opt.lambda_set, 'CV', opt.lasso_kFold);

            B0 = FitInfo.Intercept(FitInfo.IndexMinDeviance);
            coef(:,t,f) = [B0; B(:,FitInfo.IndexMinDeviance)];

            pred = ((glmval(coef(:,t,f), r_tmp(I_test,:), 'logit'))>=0.5)+1;
            correct(I_test,t) = double(pred==target(I_test));
        end
    end
else % prediction only
    n_win = length(opt.t_win)-1;
    coef = opt.coef_pretrained;
    correct = nan(length(target), n_win);
    for t = 1:n_win
        r_tmp = squeeze(mean(raster(:, tstamp>=opt.t_win(t) & tstamp<=opt.t_win(t+1), :), 2))'; % (trial, unit)
        pred = ((glmval(coef(:,t), r_tmp(:,:), 'logit'))>=0.5)+1;
        correct(:,t) = double(pred==target(:));
    end
end

% save
stat.param = param;
stat.coef = coef;
stat.correct = correct;
stat.opt = opt;

end















