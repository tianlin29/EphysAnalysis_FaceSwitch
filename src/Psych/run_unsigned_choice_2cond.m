function [fh, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond, coh, cor, session, opt)

% select data
def.session_list = [];
def.session_idx = []; % will be defined later
% process data
def.log = true;
def.constant = false;
def.verbose = true;
% plot
if length(opt.session_list)<=10; def.tile_fh_idv = true; else; def.tile_fh_idv = false; end
def.color = [0 0 0; 1 0 0];
def.average = true;
def.legend = {};
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

% plot individual subject
nsubj = length(opt.session_list);
[stat, fh_idv] = deal(cell(1, nsubj));
for s = 1:nsubj
    opt.session_idx = s;
    I = session==opt.session_list(s);
    [stat{s}, fh_idv{s}] = show_unsigned_choice_2cond(cond(I), coh(I), cor(I), opt);
end

% average across subjects
fh = show_unsigned_choice_2cond_average(stat, opt);

% show quick summary across sessions
fh_summary = show_quick_summary(stat, opt);

end


%% show_quick_summary
function fh_summary = show_quick_summary(stat, opt)

def.color = [0 0 0; 1 0 0];
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

ntrial1 = cellfun(@(x) x.ntrial1, stat(:));
ntrial2 = cellfun(@(x) x.ntrial2, stat(:));
acc1 = cellfun(@(x) round(x.acc1*100), stat(:));
acc2 = cellfun(@(x) round(x.acc2*100), stat(:));
thres1 = cellfun(@(x) round(x.thres1*100), stat(:)); 
% thres1(thres1>=100 | thres1<0) = NaN;
thres2 = cellfun(@(x) round(x.thres2*100), stat(:)); 
% thres2(thres2>=100 | thres2<0) = NaN;
tbl = table(ntrial1, ntrial2, acc1, acc2, thres1, thres2);
disp(tbl)

fh_summary = figure('pos', [50+240 300 300 120]);
subplot(1,2,1); hold on
plot(acc1, '.-', 'markers', 7, 'Color', opt.color(1,:));
plot(acc2, '.-', 'markers', 7, 'Color', opt.color(2,:));
format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Accuracy (%)', 'xlim', [0.5 length(acc1)+0.5], 'xtick', 1:length(stat))
subplot(1,2,2); hold on
plot(thres1, '.-', 'markers', 7, 'Color', opt.color(1,:));
plot(thres2, '.-', 'markers', 7, 'Color', opt.color(2,:));
format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Threshold (% Morph)', 'xlim', [0.5 length(acc1)+0.5], 'xtick', 1:length(stat))

end


%% show_unsigned_choice_2cond_average
function fh = show_unsigned_choice_2cond_average(stat, opt)

def.average = true;
def.log = true;
def.color = [0 0 0; 1 0 0];
def.legend = {};
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

% coh level check
if ~opt.average; fh = []; return; end

ncoh = cellfun(@(x) length(x.ucoh), stat);
if length(unique(ncoh))~=1
    warning('Coherence levels are not consistant across sessions. Delete some coherences.')
    for n = 1:length(stat)
        I = ismember(stat{n}.ucoh, [0 0.06 0.12 0.18 0.24 0.32 0.4 0.6 0.8 0.96]);
        stat{n}.ucoh = stat{n}.ucoh(I);
        stat{n}.lucoh = stat{n}.lucoh(I);
        stat{n}.p1 = stat{n}.p1(I);
        stat{n}.p2 = stat{n}.p2(I);
    end
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
if ~isempty(opt.legend); legend(opt.legend, 'location', 'bestoutside'); legend boxoff; end

end


%% show_unsigned_choice_2cond
function [stat, fh_idv] = show_unsigned_choice_2cond(cond, coh, cor, opt)

def.session_list = [];
def.session_idx = [];
def.log = true;
def.constant = false;
def.verbose = true;
def.tile_fh_idv = [];
def.color = [0 0 0; 1 0 0];
if exist('opt', 'var')
    opt = safeStructAssign(def, opt);
else
    opt = def;
end

% process coh
abs_coh = abs(coh); % absolute coh
ucoh = unique(abs_coh); % unique coh

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
I1 = cond==1; I2 = cond==2;
[p1, pse1] = calcGroupMean(cor(I1), abs_coh(I1), ucoh, 'binary');
[p2, pse2] = calcGroupMean(cor(I2), abs_coh(I2), ucoh, 'binary');

% fit to psychometric function, get threshold
% logit(P) = a0 + a1*cond + (a2 + cond*a3)*coh
%          = a0 + a1*cond + a2*coh + a3*coh*cond
% no constant term:
% logit(P) =                a2*coh + a3*coh*cond
if opt.constant
    [alpha, ~, stat_glmfit] = glmfit([cond==2, abs_coh, abs_coh.*(cond==2)], cor, 'binomial', 'link', 'logit'); % logit(P) = a0 + a1*cond + a2*coh + a3*coh*cond
    fresp1 = glmval(alpha, [zeros(size(fcoh)), fcoh, fcoh.*zeros(size(fcoh))], 'logit');
    fresp2 = glmval(alpha, [ones(size(fcoh)), fcoh, fcoh.*ones(size(fcoh))], 'logit');
    thres1 = fminsearch(@(abs_coh) abs(glmval(alpha, [0, abs_coh, 0], 'logit') - 0.816), .2);
    thres2 = fminsearch(@(abs_coh) abs(glmval(alpha, [1, abs_coh, abs_coh], 'logit') - 0.816), .2);
else
    [alpha, ~, stat_glmfit] = glmfit([abs_coh, abs_coh.*(cond==2)], cor, 'binomial', 'link', 'logit', 'Constant', 'off'); % logit(P) = a2*coh + a3*coh*cond
    fresp1 = glmval(alpha, [fcoh, fcoh.*zeros(size(fcoh))], 'logit', 'Constant', 'off');
    fresp2 = glmval(alpha, [fcoh, fcoh.*ones(size(fcoh))], 'logit', 'Constant', 'off');
    thres1 = fminsearch(@(abs_coh) abs(glmval(alpha, [abs_coh, 0], 'logit', 'Constant', 'off') - 0.816), .2);
    thres2 = fminsearch(@(abs_coh) abs(glmval(alpha, [abs_coh, abs_coh], 'logit', 'Constant', 'off') - 0.816), .2);
end

if opt.verbose
    fprintf('session %s\n', opt.session_list(opt.session_idx))
    if length(alpha)==4; fprintf('alpha0 (bias) = %1.3f (p = %s)\n', alpha(end-3), p2str(stat_glmfit.p(end-3))); end
    if length(alpha)==4; fprintf('alpha1 (conditional bias) = %1.3f (p = %s)\n', alpha(end-2), p2str(stat_glmfit.p(end-2))); end
    fprintf('alpha2 (slope) = %1.3f (p = %s)\n', alpha(end-1), p2str(stat_glmfit.p(end-1)));
    fprintf('alpha3 (conditional slope) = %1.3f (p = %s)\n', alpha(end), p2str(stat_glmfit.p(end)));
    fprintf('threshold of condition 1 (81.6%% correct): %1.1f%% coh\n', thres1*1e2);
    fprintf('threshold of condition 2 (81.6%% correct): %1.1f%% coh\n\n', thres2*1e2);
end

% plot
if opt.tile_fh_idv
    fh_idv = figure('color', 'w', 'position', [50+120*(opt.session_idx-1) 100 120 120]);
else
    fh_idv = figure('color', 'w', 'position', [50 100 120 120]);
end
hold on;
if opt.log
    plot(log(fcoh), fresp1, 'Color', opt.color(1,:));
    plot(log(fcoh), fresp2, 'Color', opt.color(2,:));
else
    plot(fcoh, fresp1, 'Color', opt.color(1,:));
    plot(fcoh, fresp2, 'Color', opt.color(2,:));
end
plot(lucoh, p1, '.', 'markers', 7, 'Color', opt.color(1,:));
plot(lucoh, p2, '.', 'markers', 7, 'Color', opt.color(2,:));
% cerrorbar(lucoh, p1, pse1, 'Color', opt.color(1,:));
% cerrorbar(lucoh, p2, pse2, 'Color', opt.color(2,:));

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

% stat
stat = struct('ucoh', ucoh, 'lucoh', lucoh, 'p1', p1, 'p2', p2, ... % plot data
    'fcoh', fcoh, 'fresp1', fresp1, 'fresp2', fresp2, ... % plot fit result
    'alpha', alpha(:), 'stat_glmfit', stat_glmfit, ...
    'thres1', thres1, 'thres2', thres2, ...
    'ntrial1', sum(cond==1), 'ntrial2', sum(cond==2), ... % trial number for quick summary
    'acc1', mean(cor(cond==1)), 'acc2', mean(cor(cond==2)), ... % accuracy for quick summary
    'lapse1', 1-p1(end), 'lapse2', 1-p2(end)); % lapse rate for quick summary

end



