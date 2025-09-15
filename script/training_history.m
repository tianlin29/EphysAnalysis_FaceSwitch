run('../Initialize.m');
monkey = 'Woody';
FigDir = fullfile(MainFigDir, 'training_history'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'training_history'); mkdir(InterimDir);

%% load data
D_all = load(fullfile(PreprocDir, [monkey, '_training.mat'])).trial_data;

%% calculate threshold and lapse rate
% define options
clear opt
opt.constant = false;

% get task (1-4)
task_list = {'human_monkey2_generated', 'human_monkey1_generated', 'human_monkey3_generated', 'human_monkey4_generated'};
task_list_backup = {'human_monkey2', 'human_monkey1'};
for n = 1:length(D_all)
    task = find(strcmp(task_list, D_all{n}.stim_set));
    if isempty(task); task = find(strcmp(task_list_backup, D_all{n}.stim_set)); end
    D_all{n}.task = task;
end

% for each session, get threshold and lapse rate
nses = D_all{end}.session_id;
[threshold, lapse] = deal(nan(nses, 4)); % 5 tasks

for s = 1:nses % s for session
    fprintf('%d/%d\n', s, nses)

    % load data of one session
    D = D_all(cellfun(@(x) x.session_id, D_all)==s);

    % get parameters
    task = cellfun(@(x) x.task, D);
    coh = cellfun(@(x) x.coh, D);
    resp = cellfun(@(x) x.resp, D); % there won't be NaN in resp
    if sum(isnan(resp))~=0; error('NaN in resp.'); end
    targ_cor = cellfun(@(x) x.targ_cor, D);
    cor = double(resp==targ_cor);

    for t = 1:4 % t for task
        % one session has two tasks, get one task
        I = task==t;
        if sum(I)==0; continue; end

        % coherence
        abs_coh = abs(coh(I));
        u_abs_coh = unique(abs_coh);
        fcoh = (0.011:0.01:1)';

        % correct rate
        [resp, respse] = calcGroupMean(cor(I), abs_coh, u_abs_coh, 'binary');

        % fit to get threshold
        if opt.constant
            alpha = glmfit(abs_coh, cor(I), 'binomial', 'link', 'logit'); % logit(P) = a0 + a1*coh
            fresp = glmval(alpha, fcoh, 'logit');
            threshold(s, t) = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit') - 0.816), .2);
        else
            alpha = glmfit(abs_coh, cor(I), 'binomial', 'link', 'logit', 'Constant', 'off'); % logit(P) = a1*coh
            fresp = glmval(alpha, fcoh, 'logit', 'Constant', 'off');
            threshold(s, t) = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit', 'Constant', 'off') - 0.816), .2);
        end

        % lapse rate
        lapse(s, t) = 1 - resp(end);
    end
end

% remove inaccurate threshold
threshold(threshold>1) = NaN;
threshold(threshold<0.05) = NaN;

% save
save(fullfile(InterimDir, sprintf('lapse_%s.mat', monkey)), 'lapse')
save(fullfile(InterimDir, sprintf('threshold_%s.mat', monkey)), 'threshold')

%% [poster] plot threshold of Nick and Woody
% load
threshold_Nick = load(fullfile(InterimDir, sprintf('threshold_%s.mat', 'Nick'))).threshold;
threshold_Woody = load(fullfile(InterimDir, sprintf('threshold_%s.mat', 'Woody'))).threshold;

% define color of 4 tasks
clear opt
opt.color = [0 0 0;
    44 145 224;
    58 191 153;
    240 169 58]/255;

% plot
fig_size = [280 240]; 
fh = figure('Position', [50 100 fig_size]);

% threshold (Nick)
ax1 = subplot(2,2,1); hold on
for t = 1:4
    h(t) = plot(threshold_Nick(:,t)*100, '-', 'Color', opt.color(t,:));
    scatter(1:size(threshold_Nick,1), threshold_Nick(:,t)*100, 7, 'filled', 'o', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', opt.color(t,:))
end
format_panel(ax1, 'axis', 'normal', 'ylim', [10 100], 'ytick', 10:15:100, 'xtick', 0:10:91, 'xlim', [0 91], ...
    'xlabel', '#Session', 'ylabel', {'Threshold', '(% morph)'})
ax1.Position = [0.15, 0.6, 0.8, 0.3];
axPosition = ax1.Position;
ax1.TickLength = [0.01/axPosition(3), 0.01/axPosition(4)];
title('Monkey N')

% threshold (Woody)
ax3 = subplot(2,2,3); hold on
ntr = 25;
plot(1:ntr, threshold_Woody(1:ntr,1)*100, '-', 'Color', opt.color(1,:));
scatter(1:ntr, threshold_Woody(1:ntr,1)*100, 7, 'filled', 'o', ...
    'MarkerEdgeColor', 'none', 'MarkerFaceColor', opt.color(1,:))
format_panel(ax3, 'axis', 'normal', 'ylim', [10 100], 'ytick', 10:15:100, 'xtick', 0:20:ntr, 'xlim', [0 ntr], ...
    'xlabel', '#Session', 'ylabel', {'Threshold', '(% morph)'})
ax3.Position = [0.15, 0.15, ntr/size(threshold_Woody,1)*1.12, 0.3];
axPosition = ax3.Position;
ax3.TickLength = [0.01/axPosition(4), 0.01/axPosition(4)];
xtickangle(0)
title('Monkey W')

ax4 = subplot(2,2,4); hold on
for t = 1:4
    plot(77:169, threshold_Woody(77:169,t)*100, '-', 'Color', opt.color(t,:));
    scatter(77:169, threshold_Woody(77:169,t)*100, 7, 'filled', 'o', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', opt.color(t,:))
end
format_panel(ax4, 'axis', 'normal', 'ylim', [10 100], 'ytick', 10:15:100, 'xtick', 80:20:170, 'xlim', [77 169], ...
    'xlabel', '#Session', 'ylabel', {'Threshold', '(% morph)'})
ax4.Position = [0.15+ntr/size(threshold_Woody,1)*1.12+0.02, 0.15, 93/size(threshold_Woody,1)*1.12, 0.3];
set(ax4, 'YColor', 'none');
axPosition = ax4.Position;
ax4.TickLength = [0.01/axPosition(3), 0.01/axPosition(4)];
xtickangle(0)

% print(fh, '-dpdf', fullfile(FigDir, sprintf('training_history_threshold.pdf')));

%% plot threshold and lapse rate
% load
lapse = load(fullfile(InterimDir, sprintf('lapse_%s.mat', monkey))).lapse;
threshold = load(fullfile(InterimDir, sprintf('threshold_%s.mat', monkey))).threshold;

% define color of 4 tasks
clear opt
opt.color = [0 0 0;
    44 145 224;
    58 191 153;
    240 169 58]/255;

fh = figure('Position', [50 100 280 115]); hold on
% plot([1,size(threshold,1)], [20 20], ':', 'Color', [.5 .5 .5])
h = nan(1,4);
for t = 1:4
    h(t) = plot(threshold(:,t)*100, '-', 'Color', opt.color(t,:));
    scatter(1:size(threshold,1), threshold(:,t)*100, 7, 'filled', 'o', ...
        'MarkerEdgeColor', 'white', 'MarkerFaceColor', opt.color(t,:))
end
if strcmp(monkey, 'Nick')
    format_panel(fh, 'xlabel', '#Session', 'ylabel', {'Threshold', '(% morph)'}, 'ylim', [10 100], 'axis', 'normal', ...
        'ytick', 10:15:100, 'xtick', 0:15:size(threshold,1), 'xlim', [0 95])
    title('Monkey N');
else
    format_panel(fh, 'xlabel', '#Session', 'ylabel', {'Threshold', '(% morph)'}, 'ylim', [10 100], 'axis', 'normal', ...
        'ytick', 10:15:100, 'xtick', 0:20:size(threshold,1), 'xlim', [0 170])
    title('Monkey W');
    xtickangle(0)
end
legend(h, {'Pair 1', 'Pair 2', 'Pair 3', 'Pair 4'}, 'location', 'bestoutside'); legend boxoff
print(fh, '-dpdf', fullfile(FigDir, sprintf('training_history_threshold_%s.pdf', monkey)));

% lapse rate
fh = figure('Position', [50 100 300 130]); hold on
for t = 1:4
    plot(lapse(:,t)*100, '.-', 'Color', opt.color(t,:))
end
format_panel(fh, 'xlabel', '#Session', 'ylabel', 'Lapse rate (%)', 'ylim', [0 40], 'axis', 'normal')
legend({'Pair 1', 'Pair 2', 'Pair 3', 'Pair 4'}, 'location', 'bestoutside'); legend boxoff
print(fh, '-dpdf', fullfile(FigDir, sprintf('training_history_lapse_%s.pdf', monkey)));

