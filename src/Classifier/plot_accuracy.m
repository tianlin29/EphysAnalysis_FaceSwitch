function [fh_session, fh_timebin] = plot_accuracy(coh, cor, cond, session, repeat, opt_classifier, experiment)

% get threshold
n_timebins = length(opt_classifier.tstamp{1});
n_repeats = max(repeat);
n_files = max(session);
n_conds = max(cond);

[facc] = deal(nan(n_repeats, n_files, n_timebins, n_conds)); % (repeat, session, time bin, 3 conditions)
for r = 1:n_repeats
    fprintf('repeat %d\n', r)
    for n = 1:n_files
        for t = 1:n_timebins
            I = repeat==r & session==n & ~isnan(cor(:,t));
            stat = get_threshold(cond(I), coh(I), cor(I,t));
            facc(r,n,t,:) = stat.facc;
        end
    end
end

facc_mn = squeeze(mean(facc, 1));
facc_se = squeeze(std(facc, [], 1) ./ sqrt(size(facc,1)));

ylim_range = [.5 1];

% plotting options
switch experiment
    case {'learnTask2', 'faceColor'}
        opt.color = [0 0 0; 44 145 224; 255 0 0]/255;
    case 'learnTask3'
        opt.color = [0 0 0; 58 191 153; 255 0 0]/255;
    case 'learnTask4'
        opt.color = [0 0 0; 240 169 58; 255 0 0]/255;
    otherwise
        opt.color = [0 0 0; 44 145 224; 255 0 0]/255;
end

% plot how threshold changes with session
fh_session = figure('Position', [50 100 600 600]);
for n = 1:n_files
    subplot(5,5,n); hold on
    title(sprintf('Session %d', n))
    for c = 1:n_conds
        plot(facc_mn(n,:,c), '.-', 'markers', 5, 'Color', opt.color(c,:));
        cerrorbar(1:size(facc_mn,2), facc_mn(n,:,c), facc_se(n,:,c), 'Color', opt.color(c,:));
    end
    format_panel(gca, 'xlabel', '#Time bin', 'ylabel', 'Accuracy', ...
        'xlim', [1-0.5 n_timebins+0.5], 'xtick', 1:2:n_timebins, 'ylim', ylim_range)
    xtickangle(0)
end

% plot how threshold changes with time in a trial
fh_timebin = figure('Position', [50 100 600 600]);
for t = 1:n_timebins
    subplot(5,5,t); hold on
    title(sprintf('%d-%d ms', opt_classifier.tstamp{1}(t)-opt_classifier.t_win/2, opt_classifier.tstamp{1}(t)+opt_classifier.t_win/2))
    for c = 1:n_conds
        plot(facc_mn(:,t,c), '.-', 'markers', 5, 'Color', opt.color(c,:));
        cerrorbar(1:size(facc_mn,1), facc_mn(:,t,c), facc_se(:,t,c), 'Color', opt.color(c,:));
    end
    format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Accuracy', ...
        'xlim', [1-0.5 n_files+0.5], 'xtick', 1:2:n_files, 'ylim', ylim_range)
    xtickangle(0)
end

end


%% get_threshold
function stat = get_threshold(cond, coh, cor, opt)

def.log = true;
def.constant = false;
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

% process coh
abs_coh = abs(coh); % absolute coh
ucoh = unique(abs_coh); % unique coh
ncoh = length(ucoh);

if opt.log
    lucoh = log(ucoh); % log coh
    mincoh = ucoh(1);
    if mincoh == 0
        mincoh = ucoh(2)/ucoh(3) * ucoh(2);
        lucoh(1) = log(mincoh);
    end
else
    lucoh = ucoh;
    mincoh = 0;
end

fcoh = linspace(mincoh, max(abs_coh), 100); fcoh = fcoh(:); % coh for function

% get choice accuracy
ncond = length(unique(cond));
[p, pse] = deal(nan(ncoh, ncond));
for c = 1:ncond
    [p(:,c), pse(:,c)] = calcGroupMean(cor(cond==c), abs_coh(cond==c), ucoh, 'binary');
end

% fit to psychometric function, get threshold
% logit(P) = a0 + a1*cond + (a2 + cond*a3)*coh
%          = a0 + a1*cond + a2*coh + a3*coh*cond
% no constant term:
% logit(P) =                a2*coh + a3*coh*cond
fresp = nan(length(fcoh), ncond);
thres = nan(1, ncond);
thres_se = nan(1, ncond);
for c = 1:ncond
    if opt.constant
        [alpha, ~, stat_glmfit] = glmfit([abs_coh(cond==c)], cor(cond==c), 'binomial', 'link', 'logit'); % logit(P) = a0 + a2*coh
        fresp(:,c) = glmval(alpha, fcoh, 'logit');
        thres(c) = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit') - 0.816), .2);
    else
        [alpha, ~, stat_glmfit] = glmfit([abs_coh(cond==c)], cor(cond==c), 'binomial', 'link', 'logit', 'Constant', 'off'); % logit(P) = a2*coh
        fresp(:,c) = glmval(alpha, fcoh, 'logit', 'Constant', 'off');
        thres(c) = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit', 'Constant', 'off') - 0.816), .2);

        logit_target = log(0.816 / (1 - 0.816)); 
        thres_se(c) = abs(logit_target / (alpha^2)) * sqrt(stat_glmfit.covb);
    end
end

facc = nan(1, ncond);
for c = 1:ncond
    facc(c) = mean(fresp(:,c));
end

stat = struct('ucoh', ucoh, 'lucoh', lucoh, 'ncond', ncond, 'p', p, ... % plot data
    'fcoh', fcoh, 'fresp', fresp, 'facc', facc, ... % plot fit result
    'thres', thres, 'thres_se', thres_se, ...
    'ntrial1', sum(cond==1), 'ntrial2', sum(cond==2), 'ntrial3', sum(cond==3), ... % trial number for quick summary
    'acc1', mean(cor(cond==1)), 'acc2', mean(cor(cond==2)), 'acc3', mean(cor(cond==3)), ... % accuracy for quick summary
    'lapse1', 1-p(end,1), 'lapse2', 1-p(end,2), 'lapse3', 1-p(end,3)); % lapse rate for quick summary

end