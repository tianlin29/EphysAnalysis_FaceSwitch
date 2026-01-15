function [fh, fh_idv, fh_summary, stat] = run_unsigned_choice_3cond(cond, coh, cor, session, opt)

% select data
def.session_list = [];
def.session_idx = []; % will be defined later
% process data
def.log = true;
def.constant = false;
def.verbose = true;
% plot
if length(opt.session_list)<=10; def.tile_fh_idv = true; else; def.tile_fh_idv = false; end
def.color = [0 0 0; 1 0 0; 0 1 0];
def.linewidth = [0.5, 0.5, 0.5];
def.average = true;
def.legend = {'Cond 1', 'Cond 2'};
% quick summary plot
def.normalize_threshold = false;
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

% plot individual subject
nsubj = length(opt.session_list);
[stat, fh_idv] = deal(cell(1, nsubj));
fprintf('remove %d trials with NaN in cor.\n', sum(isnan(cor)))
for s = 1:nsubj
    opt.session_idx = s;
    I = session==opt.session_list(s) & ~isnan(cor);
    [stat{s}, fh_idv{s}] = show_unsigned_choice_3cond(cond(I), coh(I), cor(I), opt);
end

% average across subjects
fh = show_unsigned_choice_3cond_average(stat, opt); % TBU

% show quick summary across sessions
fh_summary = show_quick_summary(stat, opt);

end


%% show_quick_summary
function fh_summary = show_quick_summary(stat, opt)

def.color = [0 0 0; 1 0 0; 0 1 0];
def.linewidth = [0.5, 0.5, 0.5];
def.normalize_threshold = false;
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

ntrial1 = cellfun(@(x) x.ntrial1, stat(:));
ntrial2 = cellfun(@(x) x.ntrial2, stat(:));
ntrial3 = cellfun(@(x) x.ntrial3, stat(:));
acc1 = cellfun(@(x) round(x.acc1*100), stat(:));
acc2 = cellfun(@(x) round(x.acc2*100), stat(:));
acc3 = cellfun(@(x) round(x.acc3*100), stat(:));
thres = cellfun(@(x) round(x.thres*100), stat(:), 'uni', 0); thres = cell2mat(thres);
tbl = table(ntrial1, ntrial2, ntrial3, acc1, acc2, acc3);
disp(tbl)

fh_summary = figure('pos', [50+240 300 300 120]);
ncond = stat{1}.ncond;

subplot(1,2,1); hold on
plot(acc1, '.-', 'markers', 7, 'Color', opt.color(1,:), 'LineWidth', opt.linewidth(1));
plot(acc2, '.-', 'markers', 7, 'Color', opt.color(2,:), 'LineWidth', opt.linewidth(2));
plot(acc3, '.-', 'markers', 7, 'Color', opt.color(3,:), 'LineWidth', opt.linewidth(3));
format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Accuracy (%)', 'xlim', [0.5 length(acc1)+0.5], 'xtick', 1:length(stat))
subplot(1,2,2); hold on
if opt.normalize_threshold; thres_for_plot = thres./(thres(:,1)); else; thres_for_plot = thres; end
for c = 1:ncond
    plot(thres_for_plot(:,c), '.-', 'markers', 7, 'Color', opt.color(c,:), 'LineWidth', opt.linewidth(c));
end
format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Threshold (% Morph)', 'xlim', [0.5 length(acc1)+0.5], 'xtick', 1:2:length(stat))

end


%% show_unsigned_choice_2cond_average
function fh = show_unsigned_choice_3cond_average(stat, opt)

def.average = true;
def.log = true;
def.color = [0 0 0; 1 0 0; 0 1 0];
def.linewidth = [0.5, 0.5, 0.5];
def.legend = {'Cond 1', 'Cond 2', 'Cond 3'};
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

% coh level check
if ~opt.average; fh = []; return; end

ncoh = cellfun(@(x) length(x.ucoh), stat);
if length(unique(ncoh))~=1
    warning('Coherence levels are not consistant across sessions. Cannot plot averaged result.')
    fh = [];
    return
end

% extract data
nsubj = length(stat);
ucoh = stat{1}.ucoh;
lucoh = stat{1}.lucoh;
fcoh = stat{1}.fcoh;

% data averaging
p1 = cellfun(@(x) x.p1, stat, 'uni', 0); p1 = cell2mat(p1); p1_mean = mean(p1, 2); p1_sem = std(p1, 0, 2) / sqrt(nsubj);
p2 = cellfun(@(x) x.p2, stat, 'uni', 0); p2 = cell2mat(p2); p2_mean = mean(p2, 2); p2_sem = std(p2, 0, 2) / sqrt(nsubj);

% fit curve averaging
fresp1 = cellfun(@(x) x.fresp1, stat, 'uni', 0); fresp1 = cell2mat(fresp1); fresp1 = mean(fresp1, 2);
fresp2 = cellfun(@(x) x.fresp2, stat, 'uni', 0); fresp2 = cell2mat(fresp2); fresp2 = mean(fresp2, 2);

% plot
fh = figure('color', 'w', 'Position', [50 300 120*2 120]);
hold on
if opt.log
    plot(log(fcoh), fresp1, 'Color', opt.color(1,:));
    plot(log(fcoh), fresp2, 'Color', opt.color(2,:));
else
    plot(fcoh, fresp1, 'Color', opt.color(1,:));
    plot(fcoh, fresp2, 'Color', opt.color(2,:));
end
plot(lucoh, p1_mean, '.', 'markers', 7, 'Color', opt.color(1,:));
plot(lucoh, p2_mean, '.', 'markers', 7, 'Color', opt.color(2,:));
cerrorbar(lucoh, p1_mean, p1_sem, 'Color', opt.color(1,:));
cerrorbar(lucoh, p2_mean, p2_sem, 'Color', opt.color(2,:));
if opt.log
    opt.xlim = log([3 100]/100);
    opt.xtick = log([3 6 12 24 48 96]'/100);
    opt.xticklabel = {'0', '6', '', '24', '', '96'};
else % default
    opt.xlim = [0 1];
    opt.xtick = lucoh;
    ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
    ucohl(2:2:end) = {''};
    opt.xticklabel = ucohl;
end
format_panel(gca, ...
    'ylim', [.5 1], 'ytick', .5:.125:1, ...
    'xlim', opt.xlim, 'xtick', opt.xtick, 'xticklabel', opt.xticklabel, ...
    'xlabel', {'Stimulus strength','(% Morph)'}, 'ylabel', 'P(correct)');
legend(opt.legend, 'location', 'bestoutside'); legend boxoff

end


%% show_unsigned_choice_2cond
function [stat, fh_idv] = show_unsigned_choice_3cond(cond, coh, cor, opt)

def.session_list = [];
def.session_idx = [];
def.log = true;
def.constant = false;
def.verbose = true;
def.tile_fh_idv = [];
def.color = [0 0 0; 1 0 0; 0 1 0];
def.linewidth = [0.5, 0.5, 0.5];
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
for c = 1:ncond
    if opt.constant
        [alpha, ~, stat_glmfit] = glmfit([abs_coh(cond==c)], cor(cond==c), 'binomial', 'link', 'logit'); % logit(P) = a0 + a2*coh
        fresp(:,c) = glmval(alpha, fcoh, 'logit');
        thres(c) = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit') - 0.816), .2);
    else
        [alpha, ~, stat_glmfit] = glmfit([abs_coh(cond==c)], cor(cond==c), 'binomial', 'link', 'logit', 'Constant', 'off'); % logit(P) = a2*coh
        fresp(:,c) = glmval(alpha, fcoh, 'logit', 'Constant', 'off');
        thres(c) = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit', 'Constant', 'off') - 0.816), .2);
    end
end

% plot
if opt.tile_fh_idv
    fh_idv = figure('color', 'w', 'position', [50+120*(opt.session_idx-1) 100 120 120]);
else
    fh_idv = figure('color', 'w', 'position', [50 100 120 120]);
end
hold on;
if opt.log
    for c = 1:ncond
        plot(log(fcoh), fresp(:,c), 'Color', opt.color(c,:), 'LineWidth', opt.linewidth(c));
    end
else
    for c = 1:ncond
        plot(fcoh, fresp(:,c), 'Color', opt.color(c,:), 'LineWidth', opt.linewidth(c));
    end
end
for c = 1:ncond
    plot(lucoh, p(:,c), '.', 'markers', 7, 'Color', opt.color(c,:), 'LineWidth', opt.linewidth(c));
    cerrorbar(lucoh, p(:,c), pse(:,c), 'Color', opt.color(c,:), 'LineWidth', opt.linewidth(c));
end
if opt.log
    opt.xlim = log([3 100]/100);
    opt.xtick = log([3 6 12 24 48 96]'/100);
    opt.xticklabel = {'0', '6', '', '24', '', '96'};
else % default
    opt.xlim = [0 1];
    opt.xtick = lucoh;
    ucohl = arrayfun(@(x)num2str(x, '%g'), ucoh*100, 'uni', 0);
    ucohl(2:2:end) = {''};
    opt.xticklabel = ucohl;
end
format_panel(gca, ...
    'ylim', [.5 1], 'ytick', .5:.125:1, ...
    'xlim', opt.xlim, 'xtick', opt.xtick, 'xticklabel', opt.xticklabel, ...
    'xlabel', {'Stimulus strength','(% Morph)'}, 'ylabel', 'P(correct)');
if ~isempty(opt.session_list); title(['session', num2str(opt.session_list(opt.session_idx))]); end
plot(xlim, [0.816 0.816], ':', 'Color', [.5 .5 .5])

% stat
stat = struct('ucoh', ucoh, 'lucoh', lucoh, 'ncond', ncond, 'p', p, ... % plot data
    'fcoh', fcoh, 'fresp', fresp, ... % plot fit result
    'thres', thres, ...
    'ntrial1', sum(cond==1), 'ntrial2', sum(cond==2), 'ntrial3', sum(cond==3), ... % trial number for quick summary
    'acc1', mean(cor(cond==1)), 'acc2', mean(cor(cond==2)), 'acc3', mean(cor(cond==3)), ... % accuracy for quick summary
    'lapse', 1-p(end,:)); % lapse rate for quick summary

end



